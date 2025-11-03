/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Routines for solving dense matrix equations, in some cases calling SLICOT
*/

#include <slepc/private/lmeimpl.h>     /*I "slepclme.h" I*/
#include <slepcblaslapack.h>

/*
   LMEDenseRankSVD - given a square matrix A, compute its SVD U*S*V', and determine the
   numerical rank. On exit, U contains U*S and A is overwritten with V'
*/
PetscErrorCode LMEDenseRankSVD(LME lme,PetscInt n,PetscScalar *A,PetscInt lda,PetscScalar *U,PetscInt ldu,PetscInt *rank)
{
  PetscInt       i,j,rk=0;
  PetscScalar    *work;
  PetscReal      tol,*sg,*rwork;
  PetscBLASInt   n_,lda_,ldu_,info,lw_;

  PetscFunctionBegin;
  PetscCall(PetscCalloc3(n,&sg,10*n,&work,5*n,&rwork));
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(lda,&lda_));
  PetscCall(PetscBLASIntCast(ldu,&ldu_));
  lw_ = 10*n_;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined (PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&n_,&n_,A,&lda_,sg,U,&ldu_,NULL,&n_,work,&lw_,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&n_,&n_,A,&lda_,sg,U,&ldu_,NULL,&n_,work,&lw_,rwork,&info));
#endif
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("gesvd",info);
  tol = 10*PETSC_MACHINE_EPSILON*n*sg[0];
  for (j=0;j<n;j++) {
    if (sg[j]>tol) {
      for (i=0;i<n;i++) U[i+j*n] *= sg[j];
      rk++;
    } else break;
  }
  *rank = rk;
  PetscCall(PetscFree3(sg,work,rwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_INFO)
/*
   LyapunovCholResidual - compute the residual norm ||A*U'*U+U'*U*A'+B*B'||
*/
static PetscErrorCode LyapunovCholResidual(PetscInt m,PetscScalar *A,PetscInt lda,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscBLASInt   n,kk,la,lb,lu;
  PetscScalar    *M,*R,zero=0.0,done=1.0;

  PetscFunctionBegin;
  *res = 0;
  PetscCall(PetscBLASIntCast(lda,&la));
  PetscCall(PetscBLASIntCast(ldb,&lb));
  PetscCall(PetscBLASIntCast(ldu,&lu));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(k,&kk));
  PetscCall(PetscMalloc2(m*m,&M,m*m,&R));

  /* R = B*B' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&kk,&done,B,&lb,B,&lb,&zero,R,&n));
  /* M = A*U' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,A,&la,U,&lu,&zero,M,&n));
  /* R = R+M*U */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,M,&n,U,&lu,&done,R,&n));
  /* R = R+U'*M' */
  PetscCallBLAS("BLASgemm",BLASgemm_("C","C",&n,&n,&n,&done,U,&lu,M,&n,&done,R,&n));

  *res = LAPACKlange_("F",&n,&n,R,&n,NULL);
  PetscCall(PetscFree2(M,R));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   LyapunovResidual - compute the residual norm ||A*X+X*A'+B||
*/
static PetscErrorCode LyapunovResidual(PetscInt m,PetscScalar *A,PetscInt lda,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx,PetscReal *res)
{
  PetscInt       i;
  PetscBLASInt   n,la,lb,lx;
  PetscScalar    *R,done=1.0;

  PetscFunctionBegin;
  *res = 0;
  PetscCall(PetscBLASIntCast(lda,&la));
  PetscCall(PetscBLASIntCast(ldb,&lb));
  PetscCall(PetscBLASIntCast(ldx,&lx));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscMalloc1(m*m,&R));

  /* R = B+A*X */
  for (i=0;i<m;i++) PetscCall(PetscArraycpy(R+i*m,B+i*ldb,m));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,A,&la,X,&lx,&done,R,&n));
  /* R = R+X*A' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,X,&lx,A,&la,&done,R,&n));

  *res = LAPACKlange_("F",&n,&n,R,&n,NULL);
  PetscCall(PetscFree(R));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

#if defined(SLEPC_HAVE_SLICOT)
/*
   HessLyapunovChol_SLICOT - implementation used when SLICOT is available
*/
static PetscErrorCode HessLyapunovChol_SLICOT(PetscInt m,PetscScalar *H,PetscInt ldh,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscBLASInt   lwork,info,n,kk,lu,ione=1,sdim;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*W,*wr,*wi,*work;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ldu,&lu));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(k,&kk));
  PetscCall(PetscBLASIntCast(6*m,&lwork));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

  /* transpose W = H' */
  PetscCall(PetscMalloc5(m*m,&W,m*m,&Q,m,&wr,m,&wi,lwork,&work));
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[j+i*ldh];
  }

  /* compute the real Schur form of W */
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
  SlepcCheckLapackInfo("gees",info);
#if defined(PETSC_USE_DEBUG)
  for (i=0;i<m;i++) PetscCheck(PetscRealPart(wr[i])<0.0,PETSC_COMM_SELF,PETSC_ERR_USER_INPUT,"Eigenvalue with non-negative real part, the coefficient matrix is not stable");
#endif

  /* copy B' into first rows of U */
  for (i=0;i<k;i++) {
    for (j=0;j<m;j++) U[i+j*ldu] = B[j+i*ldb];
  }

  /* solve Lyapunov equation (Hammarling) */
  PetscCallBLAS("SLICOTsb03od",SLICOTsb03od_("C","F","N",&n,&kk,W,&n,Q,&n,U,&lu,&scal,wr,wi,work,&lwork,&info));
  PetscCheck(!info,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in SLICOT subroutine SB03OD: info=%" PetscBLASInt_FMT,info);
  PetscCheck(scal==1.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Current implementation cannot handle scale factor %g",scal);

  /* resnorm = norm(H(m+1,:)*U'*U), use Q(:,1) = U'*U(:,m) */
  if (res) {
    for (j=0;j<m;j++) Q[j] = U[j+(m-1)*ldu];
    PetscCallBLAS("BLAStrmv",BLAStrmv_("U","C","N",&n,U,&lu,Q,&ione));
    PetscCheck(k==1,PETSC_COMM_SELF,PETSC_ERR_LIB,"Residual error is intended for k=1 only, but you set k=%" PetscInt_FMT,k);
    *res *= BLASnrm2_(&n,Q,&ione);
  }

  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree5(W,Q,wr,wi,work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#else

/*
   Compute the upper Cholesky factor of A
 */
static PetscErrorCode CholeskyFactor(PetscInt m,PetscScalar *A,PetscInt lda)
{
  PetscInt       i;
  PetscScalar    *S;
  PetscBLASInt   info,n,ld;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(lda,&ld));
  PetscCall(PetscMalloc1(m*m,&S));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

  /* save a copy of matrix in S */
  for (i=0;i<m;i++) PetscCall(PetscArraycpy(S+i*m,A+i*lda,m));

  /* compute upper Cholesky factor in R */
  PetscCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n,A,&ld,&info));
  PetscCall(PetscLogFlops((1.0*n*n*n)/3.0));

  if (info) {
    for (i=0;i<m;i++) {
      PetscCall(PetscArraycpy(A+i*lda,S+i*m,m));
      A[i+i*lda] += 50.0*PETSC_MACHINE_EPSILON;
    }
    PetscCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n,A,&ld,&info));
    SlepcCheckLapackInfo("potrf",info);
    PetscCall(PetscLogFlops((1.0*n*n*n)/3.0));
  }

  /* Zero out entries below the diagonal */
  for (i=0;i<m-1;i++) PetscCall(PetscArrayzero(A+i*lda+i+1,m-i-1));
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree(S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   HessLyapunovChol_LAPACK - alternative implementation when SLICOT is not available
*/
static PetscErrorCode HessLyapunovChol_LAPACK(PetscInt m,PetscScalar *H,PetscInt ldh,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscBLASInt   ilo=1,lwork,info,n,kk,lu,lb,ione=1;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*C,*W,*Z,*wr,*work,zero=0.0,done=1.0,dmone=-1.0;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *wi;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ldb,&lb));
  PetscCall(PetscBLASIntCast(ldu,&lu));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(k,&kk));
  PetscCall(PetscBLASIntCast(6*m,&lwork));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  C = U;

#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc6(m*m,&Q,m*m,&W,m*k,&Z,m,&wr,m,&wi,lwork,&work));
#else
  PetscCall(PetscMalloc5(m*m,&Q,m*m,&W,m*k,&Z,m,&wr,lwork,&work));
#endif

  /* save a copy W = H */
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[i+j*ldh];
  }

  /* compute the (real) Schur form of W */
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","I",&n,&ilo,&n,W,&n,wr,wi,Q,&n,work,&lwork,&info));
#else
  PetscCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","I",&n,&ilo,&n,W,&n,wr,Q,&n,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("hseqr",info);
#if defined(PETSC_USE_DEBUG)
  for (i=0;i<m;i++) PetscCheck(PetscRealPart(wr[i])<0.0,PETSC_COMM_SELF,PETSC_ERR_USER_INPUT,"Eigenvalue with non-negative real part %g, the coefficient matrix is not stable",(double)PetscRealPart(wr[i]));
#endif

  /* C = -Z*Z', Z = Q'*B */
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n,&kk,&n,&done,Q,&n,B,&lb,&zero,Z,&n));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&kk,&dmone,Z,&n,Z,&n,&zero,C,&lu));

  /* solve triangular Sylvester equation */
  PetscCallBLAS("LAPACKtrsyl",LAPACKtrsyl_("N","C",&ione,&n,&n,W,&n,W,&n,C,&lu,&scal,&info));
  SlepcCheckLapackInfo("trsyl",info);
  PetscCheck(scal==1.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Current implementation cannot handle scale factor %g",(double)scal);

  /* back-transform C = Q*C*Q' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,Q,&n,C,&n,&zero,W,&n));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,W,&n,Q,&n,&zero,C,&lu));

  /* resnorm = norm(H(m+1,:)*Y) */
  if (res) {
    PetscCheck(k==1,PETSC_COMM_SELF,PETSC_ERR_LIB,"Residual error is intended for k=1 only, but you set k=%" PetscInt_FMT,k);
    *res *= BLASnrm2_(&n,C+m-1,&n);
  }

  /* U = chol(C) */
  PetscCall(CholeskyFactor(m,C,ldu));

  PetscCall(PetscFPTrapPop());
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree6(Q,W,Z,wr,wi,work));
#else
  PetscCall(PetscFree5(Q,W,Z,wr,work));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif /* SLEPC_HAVE_SLICOT */

/*@C
   LMEDenseHessLyapunovChol - Computes the Cholesky factor of the solution of a
   dense Lyapunov equation with an upper Hessenberg coefficient matrix.

   Logically Collective

   Input Parameters:
+  lme - the linear matrix equation solver context
.  m   - number of rows and columns of `H`
.  H   - coefficient matrix
.  ldh - leading dimension of `H`
.  k   - number of columns of `B`
.  B   - right-hand side matrix
.  ldb - leading dimension of `B`
-  ldu - leading dimension of `U`

   Output Parameters:
+  U   - Cholesky factor of the solution
-  res - (optional) residual norm, on input it should contain $h_{m+1,m}$

   Note:
   The Lyapunov equation has the form $HX + XH^* = -BB^*$, where $H$ is an $m\times m$
   upper Hessenberg matrix, $B$ is an $m\times k$ matrix and the solution is expressed
   as $X = U^*U$, where $U$ is upper triangular. $H$ is supposed to be stable.

   When `k=1` and the `res` argument is provided, the last row of `X` is used to
   compute the residual norm of a Lyapunov equation projected via Arnoldi.

   Level: developer

.seealso: [](ch:lme), `LMEDenseLyapunov()`, `LMESolve()`
@*/
PetscErrorCode LMEDenseHessLyapunovChol(LME lme,PetscInt m,PetscScalar H[],PetscInt ldh,PetscInt k,PetscScalar B[],PetscInt ldb,PetscScalar U[],PetscInt ldu,PetscReal *res)
{
#if defined(PETSC_USE_INFO)
  PetscReal      error;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,m,2);
  PetscAssertPointer(H,3);
  PetscValidLogicalCollectiveInt(lme,ldh,4);
  PetscValidLogicalCollectiveInt(lme,k,5);
  PetscAssertPointer(B,6);
  PetscValidLogicalCollectiveInt(lme,ldb,7);
  PetscAssertPointer(U,8);
  PetscValidLogicalCollectiveInt(lme,ldu,9);
  if (res) PetscValidLogicalCollectiveReal(lme,*res,10);

#if defined(SLEPC_HAVE_SLICOT)
  PetscCall(HessLyapunovChol_SLICOT(m,H,ldh,k,B,ldb,U,ldu,res));
#else
  PetscCall(HessLyapunovChol_LAPACK(m,H,ldh,k,B,ldb,U,ldu,res));
#endif

#if defined(PETSC_USE_INFO)
  if (PetscLogPrintInfo) {
    PetscCall(LyapunovCholResidual(m,H,ldh,k,B,ldb,U,ldu,&error));
    PetscCall(PetscInfo(lme,"Residual norm of dense Lyapunov equation = %g\n",(double)error));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(SLEPC_HAVE_SLICOT)
/*
   Lyapunov_SLICOT - implementation used when SLICOT is available
*/
static PetscErrorCode Lyapunov_SLICOT(PetscInt m,PetscScalar *H,PetscInt ldh,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx)
{
  PetscBLASInt   sdim,lwork,info,n,lx,*iwork;
  PetscInt       i,j;
  PetscReal      scal,sep,ferr,*work;
  PetscScalar    *Q,*W,*wr,*wi;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ldx,&lx));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(PetscMax(20,m*m),&lwork));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

  /* transpose W = H' */
  PetscCall(PetscMalloc6(m*m,&W,m*m,&Q,m,&wr,m,&wi,m*m,&iwork,lwork,&work));
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[j+i*ldh];
  }

  /* compute the real Schur form of W */
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
  SlepcCheckLapackInfo("gees",info);

  /* copy -B into X */
  for (i=0;i<m;i++) {
    for (j=0;j<m;j++) X[i+j*ldx] = -B[i+j*ldb];
  }

  /* solve Lyapunov equation (Hammarling) */
  PetscCallBLAS("SLICOTsb03md",SLICOTsb03md_("C","X","F","N",&n,W,&n,Q,&n,X,&lx,&scal,&sep,&ferr,wr,wi,iwork,work,&lwork,&info));
  PetscCheck(!info,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in SLICOT subroutine SB03OD: info=%" PetscBLASInt_FMT,info);
  PetscCheck(scal==1.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Current implementation cannot handle scale factor %g",scal);

  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree6(W,Q,wr,wi,iwork,work));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#else

/*
   Lyapunov_LAPACK - alternative implementation when SLICOT is not available
*/
static PetscErrorCode Lyapunov_LAPACK(PetscInt m,PetscScalar *A,PetscInt lda,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx)
{
  PetscBLASInt   sdim,lwork,info,n,lx,lb,ione=1;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*W,*Z,*wr,*work,zero=0.0,done=1.0,dmone=-1.0;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#else
  PetscScalar    *wi;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ldb,&lb));
  PetscCall(PetscBLASIntCast(ldx,&lx));
  PetscCall(PetscBLASIntCast(m,&n));
  PetscCall(PetscBLASIntCast(6*m,&lwork));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc6(m*m,&Q,m*m,&W,m*m,&Z,m,&wr,m,&wi,lwork,&work));
#else
  PetscCall(PetscMalloc6(m*m,&Q,m*m,&W,m*m,&Z,m,&wr,lwork,&work,m,&rwork));
#endif

  /* save a copy W = A */
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = A[i+j*lda];
  }

  /* compute the (real) Schur form of W */
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
#else
  PetscCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,Q,&n,work,&lwork,rwork,NULL,&info));
#endif
  SlepcCheckLapackInfo("gees",info);

  /* X = -Q'*B*Q */
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n,&n,&n,&done,Q,&n,B,&lb,&zero,Z,&n));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&dmone,Z,&n,Q,&n,&zero,X,&lx));

  /* solve triangular Sylvester equation */
  PetscCallBLAS("LAPACKtrsyl",LAPACKtrsyl_("N","C",&ione,&n,&n,W,&n,W,&n,X,&lx,&scal,&info));
  SlepcCheckLapackInfo("trsyl",info);
  PetscCheck(scal==1.0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Current implementation cannot handle scale factor %g",(double)scal);

  /* back-transform X = Q*X*Q' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,Q,&n,X,&n,&zero,W,&n));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,W,&n,Q,&n,&zero,X,&lx));

  PetscCall(PetscFPTrapPop());
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree6(Q,W,Z,wr,wi,work));
#else
  PetscCall(PetscFree6(Q,W,Z,wr,work,rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif /* SLEPC_HAVE_SLICOT */

/*@C
   LMEDenseLyapunov - Computes the solution of a dense continuous-time Lyapunov
   equation.

   Logically Collective

   Input Parameters:
+  lme - the linear matrix equation solver context
.  m   - number of rows and columns of `A`
.  A   - coefficient matrix
.  lda - leading dimension of `A`
.  B   - right-hand side matrix
.  ldb - leading dimension of `B`
-  ldx - leading dimension of `X`

   Output Parameter:
.  X   - the solution

   Note:
   The Lyapunov equation has the form $AX + XA^* = -B$, where all are $m\times m$
   matrices, and $B$ is symmetric.

   Level: developer

.seealso: [](ch:lme), `LMEDenseHessLyapunovChol()`, `LMESolve()`
@*/
PetscErrorCode LMEDenseLyapunov(LME lme,PetscInt m,PetscScalar A[],PetscInt lda,PetscScalar B[],PetscInt ldb,PetscScalar X[],PetscInt ldx)
{
#if defined(PETSC_USE_INFO)
  PetscReal      error;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,m,2);
  PetscAssertPointer(A,3);
  PetscValidLogicalCollectiveInt(lme,lda,4);
  PetscAssertPointer(B,5);
  PetscValidLogicalCollectiveInt(lme,ldb,6);
  PetscAssertPointer(X,7);
  PetscValidLogicalCollectiveInt(lme,ldx,8);

#if defined(SLEPC_HAVE_SLICOT)
  PetscCall(Lyapunov_SLICOT(m,A,lda,B,ldb,X,ldx));
#else
  PetscCall(Lyapunov_LAPACK(m,A,lda,B,ldb,X,ldx));
#endif

#if defined(PETSC_USE_INFO)
  if (PetscLogPrintInfo) {
    PetscCall(LyapunovResidual(m,A,lda,B,ldb,X,ldx,&error));
    PetscCall(PetscInfo(lme,"Residual norm of dense Lyapunov equation = %g\n",(double)error));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}
