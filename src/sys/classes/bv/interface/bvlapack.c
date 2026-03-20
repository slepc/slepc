/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV private kernels that use the LAPACK
*/

#include <slepc/private/bvimpl.h>
#include <slepcblaslapack.h>

/*
    Reduction operation to compute sqrt(x**2+y**2) when normalizing vectors
*/
SLEPC_EXTERN void MPIAPI SlepcPythag(void *in,void *inout,PetscMPIInt *len,MPI_Datatype *datatype)
{
  PetscBLASInt i,n=*len;
  PetscReal    *x = (PetscReal*)in,*y = (PetscReal*)inout;

  PetscFunctionBegin;
  if (PetscUnlikely(*datatype!=MPIU_REAL)) {
    (void)(*PetscErrorPrintf)("Only implemented for MPIU_REAL data type");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  for (i=0;i<n;i++) y[i] = SlepcAbs(x[i],y[i]);
  PetscFunctionReturnVoid();
}

/*
    Compute ||A|| for an mxn matrix
*/
PetscErrorCode BVNorm_LAPACK_Private(BV bv,PetscInt m_,PetscInt n_,const PetscScalar *A,PetscInt lda_,NormType type,PetscReal *nrm,PetscBool mpi)
{
  PetscBLASInt   m,n,lda,i,j;
  PetscMPIInt    len;
  PetscReal      lnrm,*rwork=NULL,*rwork2=NULL;

  PetscFunctionBegin;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  if (type==NORM_FROBENIUS || type==NORM_2) {
    lnrm = LAPACKlange_("F",&m,&n,(PetscScalar*)A,&lda,rwork);
    if (mpi) PetscCallMPI(MPIU_Allreduce(&lnrm,nrm,1,MPIU_REAL,MPIU_LAPY2,PetscObjectComm((PetscObject)bv)));
    else *nrm = lnrm;
    PetscCall(PetscLogFlops(2.0*m*n));
  } else if (type==NORM_1) {
    if (mpi) {
      PetscCall(BVAllocateWork_Private(bv,2*n_));
      rwork = (PetscReal*)bv->work;
      rwork2 = rwork+n_;
      PetscCall(PetscArrayzero(rwork,n_));
      PetscCall(PetscArrayzero(rwork2,n_));
      for (j=0;j<n_;j++) {
        for (i=0;i<m_;i++) {
          rwork[j] += PetscAbsScalar(A[i+j*lda_]);
        }
      }
      PetscCall(PetscMPIIntCast(n_,&len));
      PetscCallMPI(MPIU_Allreduce(rwork,rwork2,len,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)bv)));
      *nrm = 0.0;
      for (j=0;j<n_;j++) if (rwork2[j] > *nrm) *nrm = rwork2[j];
    } else {
      *nrm = LAPACKlange_("O",&m,&n,(PetscScalar*)A,&lda,rwork);
    }
    PetscCall(PetscLogFlops(1.0*m*n));
  } else if (type==NORM_INFINITY) {
    PetscCall(BVAllocateWork_Private(bv,m_));
    rwork = (PetscReal*)bv->work;
    lnrm = LAPACKlange_("I",&m,&n,(PetscScalar*)A,&lda,rwork);
    if (mpi) PetscCallMPI(MPIU_Allreduce(&lnrm,nrm,1,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)bv)));
    else *nrm = lnrm;
    PetscCall(PetscLogFlops(1.0*m*n));
  }
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Normalize the columns of an mxn matrix A
*/
PetscErrorCode BVNormalize_LAPACK_Private(BV bv,PetscInt m_,PetscInt n_,const PetscScalar *A,PetscInt lda_,PetscScalar *eigi,PetscBool mpi)
{
  PetscBLASInt   m,lda,j,k,info,zero=0;
  PetscMPIInt    len;
  PetscReal      *norms,*rwork=NULL,*rwork2=NULL,done=1.0;

  PetscFunctionBegin;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(lda_,&lda));
  PetscCall(BVAllocateWork_Private(bv,2*n_));
  rwork = (PetscReal*)bv->work;
  rwork2 = rwork+n_;
  /* compute local norms */
  for (j=0;j<n_;j++) {
    k = 1;
    if (!PetscDefined(USE_COMPLEX) && eigi && eigi[j] != 0.0) k = 2;
    rwork[j] = LAPACKlange_("F",&m,&k,(PetscScalar*)(A+j*lda_),&lda,rwork2);
    if (k==2) { rwork[j+1] = rwork[j]; j++; }
  }
  /* reduction to get global norms */
  if (mpi) {
    PetscCall(PetscMPIIntCast(n_,&len));
    PetscCall(PetscArrayzero(rwork2,n_));
    PetscCallMPI(MPIU_Allreduce(rwork,rwork2,len,MPIU_REAL,MPIU_LAPY2,PetscObjectComm((PetscObject)bv)));
    norms = rwork2;
  } else norms = rwork;
  /* scale columns */
  for (j=0;j<n_;j++) {
    k = 1;
#if !defined(PETSC_USE_COMPLEX)
    if (eigi && eigi[j] != 0.0) k = 2;
#endif
    PetscCallBLAS("LAPACKlascl",LAPACKlascl_("G",&zero,&zero,norms+j,&done,&m,&k,(PetscScalar*)(A+j*lda_),&lda,&info));
    SlepcCheckLapackInfo("lascl",info);
    if (k==2) j++;
  }
  PetscCall(PetscLogFlops(3.0*m*n_));
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute the upper Cholesky factor in R and its inverse in S.
   If S == R then the inverse overwrites the Cholesky factor.
 */
PetscErrorCode BVMatCholInv_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscInt       i,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS;
  PetscBLASInt   info,n_ = 0,m_ = 0,ld_,lds_;

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  PetscCall(MatGetSize(R,&m,NULL));
  n = k-l;
  PetscCall(PetscBLASIntCast(m,&m_));
  PetscCall(PetscBLASIntCast(n,&n_));
  ld  = m;
  ld_ = m_;
  PetscCall(MatDenseGetArray(R,&pR));

  if (S==R) {
    PetscCall(BVAllocateWork_Private(bv,m*k));
    pS = bv->work;
    lds = ld;
    lds_ = ld_;
  } else {
    PetscCall(MatDenseGetArray(S,&pS));
    PetscCall(MatGetSize(S,&lds,NULL));
    PetscCall(PetscBLASIntCast(lds,&lds_));
  }

  /* save a copy of matrix in S */
  for (i=l;i<k;i++) PetscCall(PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n));

  /* compute upper Cholesky factor in R */
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n_,pR+l*ld+l,&ld_,&info));
  PetscCall(PetscLogFlops((1.0*n*n*n)/3.0));

  if (info) {  /* LAPACKpotrf failed, retry on diagonally perturbed matrix */
    for (i=l;i<k;i++) {
      PetscCall(PetscArraycpy(pR+i*ld+l,pS+i*lds+l,n));
      pR[i+i*ld] += 50.0*PETSC_MACHINE_EPSILON;
    }
    PetscCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n_,pR+l*ld+l,&ld_,&info));
    SlepcCheckLapackInfo("potrf",info);
    PetscCall(PetscLogFlops((1.0*n*n*n)/3.0));
  }

  /* compute S = inv(R) */
  if (S==R) {
    PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pR+l*ld+l,&ld_,&info));
  } else {
    PetscCall(PetscArrayzero(pS+l*lds,(k-l)*k));
    for (i=l;i<k;i++) PetscCall(PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n));
    PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pS+l*lds+l,&lds_,&info));
  }
  SlepcCheckLapackInfo("trtri",info);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogFlops(0.33*n*n*n));

  /* Zero out entries below the diagonal */
  for (i=l;i<k-1;i++) {
    PetscCall(PetscArrayzero(pR+i*ld+i+1,(k-i-1)));
    if (S!=R) PetscCall(PetscArrayzero(pS+i*lds+i+1,(k-i-1)));
  }
  PetscCall(MatDenseRestoreArray(R,&pR));
  if (S!=R) PetscCall(MatDenseRestoreArray(S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute the inverse of an upper triangular matrix R, store it in S.
   If S == R then the inverse overwrites R.
 */
PetscErrorCode BVMatTriInv_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscInt       i,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS;
  PetscBLASInt   info,n_,m_ = 0,ld_,lds_;

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  PetscCall(MatGetSize(R,&m,NULL));
  n = k-l;
  PetscCall(PetscBLASIntCast(m,&m_));
  PetscCall(PetscBLASIntCast(n,&n_));
  ld  = m;
  ld_ = m_;
  PetscCall(MatDenseGetArray(R,&pR));

  if (S==R) {
    PetscCall(BVAllocateWork_Private(bv,m*k));
    pS = bv->work;
    lds = ld;
    lds_ = ld_;
  } else {
    PetscCall(MatDenseGetArray(S,&pS));
    PetscCall(MatGetSize(S,&lds,NULL));
    PetscCall(PetscBLASIntCast(lds,&lds_));
  }

  /* compute S = inv(R) */
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  if (S==R) {
    PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pR+l*ld+l,&ld_,&info));
  } else {
    PetscCall(PetscArrayzero(pS+l*lds,(k-l)*k));
    for (i=l;i<k;i++) PetscCall(PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n));
    PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pS+l*lds+l,&lds_,&info));
  }
  SlepcCheckLapackInfo("trtri",info);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscLogFlops(0.33*n*n*n));

  PetscCall(MatDenseRestoreArray(R,&pR));
  if (S!=R) PetscCall(MatDenseRestoreArray(S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Compute the matrix to be used for post-multiplying the basis in the SVQB
   block orthogonalization method.
   On input R = V'*V, on output S = D*U*Lambda^{-1/2} where (U,Lambda) is
   the eigendecomposition of D*R*D with D=diag(R)^{-1/2}.
   If S == R then the result overwrites R.
 */
PetscErrorCode BVMatSVQB_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscInt       i,j,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS,*D,*work,a;
  PetscReal      *eig,dummy;
  PetscBLASInt   info,lwork,n_,m_ = 0,ld_,lds_;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,rdummy;
#endif

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  PetscCall(MatGetSize(R,&m,NULL));
  PetscCall(MatDenseGetLDA(R,&ld));
  n = k-l;
  PetscCall(PetscBLASIntCast(m,&m_));
  PetscCall(PetscBLASIntCast(n,&n_));
  ld_ = m_;
  PetscCall(MatDenseGetArray(R,&pR));

  if (S==R) {
    pS = pR;
    lds = ld;
    lds_ = ld_;
  } else {
    PetscCall(MatDenseGetArray(S,&pS));
    PetscCall(MatDenseGetLDA(S,&lds));
    PetscCall(PetscBLASIntCast(lds,&lds_));
  }
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS,&lds_,&dummy,&a,&lwork,&rdummy,&info));
  PetscCall(PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork));
  PetscCall(PetscMalloc4(n,&eig,n,&D,lwork,&work,PetscMax(1,3*n-2),&rwork));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS,&lds_,&dummy,&a,&lwork,&info));
  PetscCall(PetscBLASIntCast((PetscInt)a,&lwork));
  PetscCall(PetscMalloc3(n,&eig,n,&D,lwork,&work));
#endif

  /* copy and scale matrix */
  for (i=l;i<k;i++) D[i-l] = 1.0/PetscSqrtReal(PetscRealPart(pR[i+i*ld]));
  for (i=l;i<k;i++) for (j=l;j<k;j++) pS[i+j*lds] = pR[i+j*ld]*D[i-l];
  for (j=l;j<k;j++) for (i=l;i<k;i++) pS[i+j*lds] *= D[j-l];

  /* compute eigendecomposition */
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS+l*lds+l,&lds_,eig,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS+l*lds+l,&lds_,eig,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("syev",info);

  if (S!=R) {   /* R = U' */
    for (i=l;i<k;i++) for (j=l;j<k;j++) pR[i+j*ld] = pS[j+i*lds];
  }

  /* compute S = D*U*Lambda^{-1/2} */
  for (i=l;i<k;i++) for (j=l;j<k;j++) pS[i+j*lds] *= D[i-l];
  for (j=l;j<k;j++) for (i=l;i<k;i++) pS[i+j*lds] /= PetscSqrtReal(eig[j-l]);

  if (S!=R) {   /* compute R = inv(S) = Lambda^{1/2}*U'/D */
    for (i=l;i<k;i++) for (j=l;j<k;j++) pR[i+j*ld] *= PetscSqrtReal(eig[i-l]);
    for (j=l;j<k;j++) for (i=l;i<k;i++) pR[i+j*ld] /= D[j-l];
  }

#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree4(eig,D,work,rwork));
#else
  PetscCall(PetscFree3(eig,D,work));
#endif
  PetscCall(PetscLogFlops(9.0*n*n*n));
  PetscCall(PetscFPTrapPop());

  PetscCall(MatDenseRestoreArray(R,&pR));
  if (S!=R) PetscCall(MatDenseRestoreArray(S,&pS));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    QR factorization of an mxn matrix via parallel TSQR
*/
PetscErrorCode BVOrthogonalize_LAPACK_TSQR(BV bv,PetscInt m_,PetscInt n_,PetscScalar *Q,PetscInt ldq_,PetscScalar *R,PetscInt ldr)
{
  PetscInt       level,plevel,nlevels,lda,worklen;
  PetscBLASInt   m,n,ldq,i,j,k,l,nb,sz,lwork,info;
  PetscScalar    *tau,*work,*A=NULL,*QQ=NULL,*Qhalf,*C=NULL,one=1.0,zero=0.0;
  PetscMPIInt    rank,size,count,stride,powtwo,s = 0;
  MPI_Datatype   tmat;

  PetscFunctionBegin;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(ldq_,&ldq));
  k  = PetscMin(m,n);
  nb = 16;
  lda = 2*n;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)bv),&rank));
  nlevels = (PetscInt)PetscCeilReal(PetscLog2Real((PetscReal)size));
  PetscCall(PetscMPIIntCast(PetscPowInt(2,(PetscInt)PetscFloorReal(PetscLog2Real((PetscReal)size))),&powtwo));
  worklen = n+n*nb;
  if (nlevels) worklen += n*lda+n*lda*nlevels+n*lda;
  PetscCall(BVAllocateWork_Private(bv,worklen));
  tau  = bv->work;
  work = bv->work+n;
  PetscCall(PetscBLASIntCast(n*nb,&lwork));
  if (nlevels) {
    A  = bv->work+n+n*nb;
    QQ = bv->work+n+n*nb+n*lda;
    C  = bv->work+n+n*nb+n*lda+n*lda*nlevels;
  }

  /* Compute QR */
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&m,&n,Q,&ldq,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("geqrf",info);

  /* Extract R */
  if (R || nlevels) {
    for (j=0;j<n;j++) {
      for (i=0;i<=PetscMin(j,m-1);i++) {
        if (nlevels) A[i+j*lda] = Q[i+j*ldq];
        else R[i+j*ldr] = Q[i+j*ldq];
      }
      for (i=PetscMin(j,m-1)+1;i<n;i++) {
        if (nlevels) A[i+j*lda] = 0.0;
        else R[i+j*ldr] = 0.0;
      }
    }
  }

  /* Compute orthogonal matrix in Q */
  PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&m,&k,&k,Q,&ldq,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("orgqr",info);

  if (nlevels) {

    PetscCall(PetscMPIIntCast(n,&count));
    PetscCall(PetscMPIIntCast(lda,&stride));
    PetscCall(PetscBLASIntCast(lda,&l));
    PetscCallMPI(MPI_Type_vector(count,count,stride,MPIU_SCALAR,&tmat));
    PetscCallMPI(MPI_Type_commit(&tmat));

    for (level=nlevels;level>=1;level--) {

      plevel = PetscPowInt(2,level);
      PetscCall(PetscMPIIntCast(plevel*PetscFloorReal(rank/(PetscReal)plevel)+(rank+PetscPowInt(2,level-1))%plevel,&s));

      /* Stack triangular matrices */
      if (rank<s && s<size) {  /* send top part, receive bottom part */
        PetscCallMPI(MPI_Sendrecv(A,1,tmat,s,111,A+n,1,tmat,s,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE));
      } else if (s<size) {  /* copy top to bottom, receive top part */
        PetscCallMPI(MPI_Sendrecv(A,1,tmat,rank,111,A+n,1,tmat,rank,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE));
        PetscCallMPI(MPI_Sendrecv(A+n,1,tmat,s,111,A,1,tmat,s,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE));
      }
      if (level<nlevels && size!=powtwo) {  /* for cases when size is not a power of 2 */
        if (rank<size-powtwo) {  /* send bottom part */
          PetscCallMPI(MPI_Send(A+n,1,tmat,rank+powtwo,111,PetscObjectComm((PetscObject)bv)));
        } else if (rank>=powtwo) {  /* receive bottom part */
          PetscCallMPI(MPI_Recv(A+n,1,tmat,rank-powtwo,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE));
        }
      }
      /* Compute QR and build orthogonal matrix */
      if (level<nlevels || (level==nlevels && s<size)) {
        PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&l,&n,A,&l,tau,work,&lwork,&info));
        SlepcCheckLapackInfo("geqrf",info);
        PetscCall(PetscArraycpy(QQ+(level-1)*n*lda,A,n*lda));
        PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&l,&n,&n,QQ+(level-1)*n*lda,&l,tau,work,&lwork,&info));
        SlepcCheckLapackInfo("orgqr",info);
        for (j=0;j<n;j++) {
          for (i=j+1;i<n;i++) A[i+j*lda] = 0.0;
        }
      } else if (level==nlevels) {  /* only one triangular matrix, set Q=I */
        PetscCall(PetscArrayzero(QQ+(level-1)*n*lda,n*lda));
        for (j=0;j<n;j++) QQ[j+j*lda+(level-1)*n*lda] = 1.0;
      }
    }

    /* Extract R */
    if (R) {
      for (j=0;j<n;j++) {
        for (i=0;i<=j;i++) R[i+j*ldr] = A[i+j*lda];
        for (i=j+1;i<n;i++) R[i+j*ldr] = 0.0;
      }
    }

    /* Accumulate orthogonal matrices */
    for (level=1;level<=nlevels;level++) {
      plevel = PetscPowInt(2,level);
      PetscCall(PetscMPIIntCast(plevel*PetscFloorReal(rank/(PetscReal)plevel)+(rank+PetscPowInt(2,level-1))%plevel,&s));
      Qhalf = (rank<s)? QQ+(level-1)*n*lda: QQ+(level-1)*n*lda+n;
      if (level<nlevels) {
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&n,&one,QQ+level*n*lda,&l,Qhalf,&l,&zero,C,&l));
        PetscCall(PetscArraycpy(QQ+level*n*lda,C,n*lda));
      } else {
        for (i=0;i<m/l;i++) {
          PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&n,&one,Q+i*l,&ldq,Qhalf,&l,&zero,C,&l));
          for (j=0;j<n;j++) PetscCall(PetscArraycpy(Q+i*l+j*ldq,C+j*l,l));
        }
        sz = m%l;
        if (sz) {
          PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&sz,&n,&n,&one,Q+(m/l)*l,&ldq,Qhalf,&l,&zero,C,&l));
          for (j=0;j<n;j++) PetscCall(PetscArraycpy(Q+(m/l)*l+j*ldq,C+j*l,sz));
        }
      }
    }

    PetscCallMPI(MPI_Type_free(&tmat));
  }

  PetscCall(PetscLogFlops(3.0*m*n*n));
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Reduction operation to compute [~,Rout]=qr([Rin1;Rin2]) in the TSQR algorithm;
    all matrices are upper triangular stored in packed format
*/
SLEPC_EXTERN void MPIAPI SlepcGivensPacked(void *in,void *inout,PetscMPIInt *len,MPI_Datatype *datatype)
{
  PetscBLASInt   n,i,j,k,one=1;
  PetscMPIInt    tsize;
  PetscScalar    v,s,*R2=(PetscScalar*)in,*R1=(PetscScalar*)inout;
  PetscReal      c;

  PetscFunctionBegin;
  PetscCallMPIAbort(PETSC_COMM_SELF,MPI_Type_size(*datatype,&tsize));  /* we assume len=1 */
  tsize /= sizeof(PetscScalar);
  n = (-1+(PetscBLASInt)PetscSqrtReal(1+8*tsize))/2;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      LAPACKlartg_(R1+(2*n-j-1)*j/2+j,R2+(2*n-i-1)*i/2+j,&c,&s,&v);
      R1[(2*n-j-1)*j/2+j] = v;
      k = n-j-1;
      if (k) BLASrot_(&k,R1+(2*n-j-1)*j/2+j+1,&one,R2+(2*n-i-1)*i/2+j+1,&one,&c,&s);
    }
  }
  PetscFunctionReturnVoid();
}

/*
    Computes the R factor of the QR factorization of an mxn matrix via parallel TSQR
*/
PetscErrorCode BVOrthogonalize_LAPACK_TSQR_OnlyR(BV bv,PetscInt m_,PetscInt n_,PetscScalar *Q,PetscInt ldq_,PetscScalar *R,PetscInt ldr)
{
  PetscInt       worklen;
  PetscBLASInt   m,n,ldq,i,j,s,nb,lwork,info;
  PetscScalar    *tau,*work,*A=NULL,*R1=NULL,*R2=NULL;
  PetscMPIInt    size,count;
  MPI_Datatype   tmat;

  PetscFunctionBegin;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCall(PetscBLASIntCast(m_,&m));
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(ldq_,&ldq));
  nb = 16;
  s  = n+n*(n-1)/2;  /* length of packed triangular matrix */
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size));
  worklen = n+n*nb+2*s+ldq*n;
  PetscCall(BVAllocateWork_Private(bv,worklen));
  tau  = bv->work;
  work = bv->work+n;
  R1   = bv->work+n+n*nb;
  R2   = bv->work+n+n*nb+s;
  A    = bv->work+n+n*nb+2*s;
  PetscCall(PetscBLASIntCast(n*nb,&lwork));
  PetscCall(PetscArraycpy(A,Q,ldq*n));

  /* Compute QR */
  PetscCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&m,&n,A,&ldq,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("geqrf",info);

  if (size==1) {
    /* Extract R */
    for (j=0;j<n;j++) {
      for (i=0;i<=PetscMin(j,m-1);i++) R[i+j*ldr] = A[i+j*ldq];
      for (i=PetscMin(j,m-1)+1;i<n;i++) R[i+j*ldr] = 0.0;
    }
  } else {
    /* Use MPI reduction operation to obtain global R */
    PetscCall(PetscMPIIntCast(s,&count));
    PetscCallMPI(MPI_Type_contiguous(count,MPIU_SCALAR,&tmat));
    PetscCallMPI(MPI_Type_commit(&tmat));
    for (i=0;i<n;i++) {
      for (j=i;j<n;j++) R1[(2*n-i-1)*i/2+j] = (i<m)?A[i+j*ldq]:0.0;
    }
    PetscCallMPI(MPIU_Allreduce(R1,R2,1,tmat,MPIU_TSQR,PetscObjectComm((PetscObject)bv)));
    for (i=0;i<n;i++) {
      for (j=0;j<i;j++) R[i+j*ldr] = 0.0;
      for (j=i;j<n;j++) R[i+j*ldr] = R2[(2*n-i-1)*i/2+j];
    }
    PetscCallMPI(MPI_Type_free(&tmat));
  }

  PetscCall(PetscLogFlops(3.0*m*n*n));
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}
