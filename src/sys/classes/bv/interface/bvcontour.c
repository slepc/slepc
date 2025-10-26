/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV developer functions needed in contour integral methods
*/

#include <slepc/private/bvimpl.h>            /*I "slepcbv.h" I*/
#include <slepcblaslapack.h>

#define p_id(i) (i*subcomm->n + subcomm->color)

/*@
   BVScatter - Scatters the columns of a `BV` to another `BV` created in a
   subcommunicator.

   Collective

   Input Parameters:
+  Vin  - input basis vectors (defined on the whole communicator)
.  scat - `VecScatter` object that contains the info for the communication
-  xdup - an auxiliary vector

   Output Parameter:
.  Vout - output basis vectors (defined on the subcommunicator)

   Notes:
   Currently implemented as a loop for each active column, where each
   column is scattered independently. The vector `xdup` is defined on the
   contiguous parent communicator and should have enough space to store one
   duplicate of the original vector per each subcommunicator.

   Level: developer

.seealso: [](sec:bv), `BVGetColumn()`
@*/
PetscErrorCode BVScatter(BV Vin,BV Vout,VecScatter scat,Vec xdup)
{
  PetscInt          i;
  Vec               v,z;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vin,BV_CLASSID,1);
  PetscValidHeaderSpecific(Vout,BV_CLASSID,2);
  PetscValidHeaderSpecific(scat,PETSCSF_CLASSID,3);
  PetscValidHeaderSpecific(xdup,VEC_CLASSID,4);
  PetscCall(BVCreateVec(Vout,&z));
  for (i=Vin->l;i<Vin->k;i++) {
    PetscCall(BVGetColumn(Vin,i,&v));
    PetscCall(VecScatterBegin(scat,v,xdup,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scat,v,xdup,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(BVRestoreColumn(Vin,i,&v));
    PetscCall(VecGetArrayRead(xdup,&array));
    PetscCall(VecPlaceArray(z,array));
    PetscCall(BVInsertVec(Vout,i,z));
    PetscCall(VecResetArray(z));
    PetscCall(VecRestoreArrayRead(xdup,&array));
  }
  PetscCall(VecDestroy(&z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSumQuadrature - Computes the sum of terms required in the quadrature
   rule to approximate the contour integral.

   Collective

   Input Parameters:
+  Y       - input basis vectors
.  M       - number of moments
.  L       - block size
.  L_max   - maximum block size
.  w       - quadrature weights
.  zn      - normalized quadrature points
.  scat    - (optional) `VecScatter` object to communicate between subcommunicators
.  subcomm - subcommunicator layout
.  npoints - number of points to process by the subcommunicator
-  useconj - whether conjugate points can be used or not

   Output Parameter:
.  S       - output basis vectors

   Notes:
   This is a generalization of `BVMult()`. The resulting matrix `S` consists of `M`
   panels of `L` columns, and the following formula is computed for each panel
   $S_k = \sum_j w_j\zeta_j^kY_j$, where $Y_j$ is the $j$-th panel of $Y$ containing
   the result of solving $T(z_j)^{-1}X$ for each integration point $j$. `L_max` is
   the width of the panels in $Y$.

   When using subcommunicators, `Y` is stored in the subcommunicators for a subset
   of integration points. In that case, the computation is done in the `subcomm`
   and then scattered to the whole communicator in `S` using the `VecScatter` `scat`.
   The value `npoints` is the number of points to be processed in this `subcomm`
   and the flag `useconj` indicates whether symmetric points can be reused.

   Level: developer

.seealso: [](sec:bv), `BVMult()`, `BVScatter()`, `BVDotQuadrature()`, `RGComputeQuadrature()`, `RGCanUseConjugates()`
@*/
PetscErrorCode BVSumQuadrature(BV S,BV Y,PetscInt M,PetscInt L,PetscInt L_max,PetscScalar w[],PetscScalar zn[],VecScatter scat,PetscSubcomm subcomm,PetscInt npoints,PetscBool useconj)
{
  PetscInt       i,j,k,nloc;
  Vec            v,sj;
  PetscScalar    *ppk,*pv,one=1.0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,BV_CLASSID,1);
  PetscValidHeaderSpecific(Y,BV_CLASSID,2);
  if (scat) PetscValidHeaderSpecific(scat,PETSCSF_CLASSID,8);

  PetscCall(BVGetSizes(Y,&nloc,NULL,NULL));
  PetscCall(PetscMalloc1(npoints,&ppk));
  for (i=0;i<npoints;i++) ppk[i] = 1.0;
  PetscCall(BVCreateVec(Y,&v));
  for (k=0;k<M;k++) {
    for (j=0;j<L;j++) {
      PetscCall(VecSet(v,0.0));
      for (i=0;i<npoints;i++) {
        PetscCall(BVSetActiveColumns(Y,i*L_max+j,i*L_max+j+1));
        PetscCall(BVMultVec(Y,ppk[i]*w[p_id(i)],1.0,v,&one));
      }
      if (PetscUnlikely(useconj)) {
        PetscCall(VecGetArray(v,&pv));
        for (i=0;i<nloc;i++) pv[i] = 2.0*PetscRealPart(pv[i]);
        PetscCall(VecRestoreArray(v,&pv));
      }
      PetscCall(BVGetColumn(S,k*L+j,&sj));
      if (PetscUnlikely(scat)) {
        PetscCall(VecScatterBegin(scat,v,sj,ADD_VALUES,SCATTER_REVERSE));
        PetscCall(VecScatterEnd(scat,v,sj,ADD_VALUES,SCATTER_REVERSE));
      } else PetscCall(VecCopy(v,sj));
      PetscCall(BVRestoreColumn(S,k*L+j,&sj));
    }
    for (i=0;i<npoints;i++) ppk[i] *= zn[p_id(i)];
  }
  PetscCall(PetscFree(ppk));
  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDotQuadrature - Computes the projection terms required in the quadrature
   rule to approximate the contour integral.

   Collective

   Input Parameters:
+  Y       - first basis vectors
.  V       - second basis vectors
.  M       - number of moments
.  L       - block size
.  L_max   - maximum block size
.  w       - quadrature weights
.  zn      - normalized quadrature points
.  subcomm - subcommunicator layout
.  npoints - number of points to process by the subcommunicator
-  useconj - whether conjugate points can be used or not

   Output Parameter:
.  Mu      - computed result

   Notes:
   This is a generalization of `BVDot()`. The resulting matrix `Mu` consists of `M`
   blocks of size `L*L` (placed horizontally), each of them computed as
   $Mu_k = \sum_j w_j\zeta_j^k V^*Y_j$, where $Y_j$ is the $j$-th panel of $Y$ containing
   the result of solving $T(z_j)^{-1}X$ for each integration point $j$. `L_max` is
   the width of the panels in $Y$.

   When using subcommunicators, $Y$ is stored in the subcommunicators for a subset
   of integration points. In that case, the computation is done in the `subcomm`
   and then the final result is combined via reduction.
   The value `npoints` is the number of points to be processed in this `subcomm`
   and the flag `useconj` indicates whether symmetric points can be reused.

   Level: developer

.seealso: [](sec:bv), `BVDot()`, `BVScatter()`, `BVSumQuadrature()`, `RGComputeQuadrature()`, `RGCanUseConjugates()`
@*/
PetscErrorCode BVDotQuadrature(BV Y,BV V,PetscScalar Mu[],PetscInt M,PetscInt L,PetscInt L_max,PetscScalar w[],PetscScalar zn[],PetscSubcomm subcomm,PetscInt npoints,PetscBool useconj)
{
  PetscMPIInt       sub_size,count;
  PetscInt          i,j,k,s;
  PetscScalar       *temp,*temp2,*ppk,alp;
  Mat               H;
  const PetscScalar *pH;
  MPI_Comm          child,parent;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Y,BV_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);

  PetscCall(PetscSubcommGetChild(subcomm,&child));
  PetscCallMPI(MPI_Comm_size(child,&sub_size));
  PetscCall(PetscMalloc3(npoints*L*(L+1),&temp,2*M*L*L,&temp2,npoints,&ppk));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,L,L_max*npoints,NULL,&H));
  PetscCall(PetscArrayzero(temp2,2*M*L*L));
  PetscCall(BVSetActiveColumns(Y,0,L_max*npoints));
  PetscCall(BVSetActiveColumns(V,0,L));
  PetscCall(BVDot(Y,V,H));
  PetscCall(MatDenseGetArrayRead(H,&pH));
  for (i=0;i<npoints;i++) {
    for (j=0;j<L;j++) {
      for (k=0;k<L;k++) {
        temp[k+j*L+i*L*L] = pH[k+j*L+i*L*L_max];
      }
    }
  }
  PetscCall(MatDenseRestoreArrayRead(H,&pH));
  for (i=0;i<npoints;i++) ppk[i] = 1;
  for (k=0;k<2*M;k++) {
    for (j=0;j<L;j++) {
      for (i=0;i<npoints;i++) {
        alp = ppk[i]*w[p_id(i)];
        for (s=0;s<L;s++) {
          if (!useconj) temp2[s+(j+k*L)*L] += alp*temp[s+(j+i*L)*L];
          else temp2[s+(j+k*L)*L] += 2.0*PetscRealPart(alp*temp[s+(j+i*L)*L]);
        }
      }
    }
    for (i=0;i<npoints;i++) ppk[i] *= zn[p_id(i)];
  }
  for (i=0;i<2*M*L*L;i++) temp2[i] /= sub_size;
  PetscCall(PetscMPIIntCast(2*M*L*L,&count));
  PetscCall(PetscSubcommGetParent(subcomm,&parent));
  PetscCallMPI(MPIU_Allreduce(temp2,Mu,count,MPIU_SCALAR,MPIU_SUM,parent));
  PetscCall(PetscFree3(temp,temp2,ppk));
  PetscCall(MatDestroy(&H));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVTraceQuadrature - Computes an estimate of the number of eigenvalues
   inside a region via quantities computed in the quadrature rule of
   contour integral methods.

   Collective

   Input Parameters:
+  Y       - first basis vectors
.  V       - second basis vectors
.  L       - block size
.  L_max   - maximum block size
.  w       - quadrature weights
.  scat    - (optional) `VecScatter` object to communicate between subcommunicators
.  subcomm - subcommunicator layout
.  npoints - number of points to process by the subcommunicator
-  useconj - whether conjugate points can be used or not

   Output Parameter:
.  est_eig - estimated eigenvalue count

   Notes:
   This function returns an estimation of the number of eigenvalues in the
   region, computed as $\operatorname{tr}(V^*S_0)$, where $S_0$ is the first panel of $S$
   computed by `BVSumQuadrature()`.

   When using subcommunicators, `Y` is stored in the subcommunicators for a subset
   of integration points. In that case, the computation is done in the `subcomm`
   and then scattered to the whole communicator in `V` using the `VecScatter` `scat`.
   The value `npoints` is the number of points to be processed in this `subcomm`
   and the flag `useconj` indicates whether symmetric points can be reused.

   Level: developer

.seealso: [](sec:bv), `BVScatter()`, `BVDotQuadrature()`, `BVSumQuadrature()`, `RGComputeQuadrature()`, `RGCanUseConjugates()`
@*/
PetscErrorCode BVTraceQuadrature(BV Y,BV V,PetscInt L,PetscInt L_max,PetscScalar w[],VecScatter scat,PetscSubcomm subcomm,PetscInt npoints,PetscBool useconj,PetscReal *est_eig)
{
  PetscInt       i,j;
  Vec            y,yall,vj;
  PetscScalar    dot,sum=0.0,one=1.0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Y,BV_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);
  if (scat) PetscValidHeaderSpecific(scat,PETSCSF_CLASSID,6);

  PetscCall(BVCreateVec(Y,&y));
  PetscCall(BVCreateVec(V,&yall));
  for (j=0;j<L;j++) {
    PetscCall(VecSet(y,0.0));
    for (i=0;i<npoints;i++) {
      PetscCall(BVSetActiveColumns(Y,i*L_max+j,i*L_max+j+1));
      PetscCall(BVMultVec(Y,w[p_id(i)],1.0,y,&one));
    }
    PetscCall(BVGetColumn(V,j,&vj));
    if (scat) {
      PetscCall(VecScatterBegin(scat,y,yall,ADD_VALUES,SCATTER_REVERSE));
      PetscCall(VecScatterEnd(scat,y,yall,ADD_VALUES,SCATTER_REVERSE));
      PetscCall(VecDot(vj,yall,&dot));
    } else PetscCall(VecDot(vj,y,&dot));
    PetscCall(BVRestoreColumn(V,j,&vj));
    if (useconj) sum += 2.0*PetscRealPart(dot);
    else sum += dot;
  }
  *est_eig = PetscAbsScalar(sum)/(PetscReal)L;
  PetscCall(VecDestroy(&y));
  PetscCall(VecDestroy(&yall));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVSVDAndRank_Refine(BV S,PetscReal delta,PetscScalar *pA,PetscReal *sigma,PetscInt *rank)
{
  PetscInt       i,j,k,ml=S->k;
  PetscMPIInt    len;
  PetscScalar    *work,*B,*tempB,*sarray,*Q1,*Q2,*temp2,alpha=1.0,beta=0.0;
  PetscBLASInt   l,m,n,lda,ldu,ldvt,lwork,info,ldb,ldc,lds;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(S->ld,&lds));
  PetscCall(BVGetArray(S,&sarray));
  PetscCall(PetscMalloc6(ml*ml,&temp2,S->n*ml,&Q1,S->n*ml,&Q2,ml*ml,&B,ml*ml,&tempB,5*ml,&work));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc1(5*ml,&rwork));
#endif
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));

  PetscCall(PetscArrayzero(B,ml*ml));
  for (i=0;i<ml;i++) B[i*ml+i]=1;

  for (k=0;k<2;k++) {
    PetscCall(PetscBLASIntCast(S->n,&m));
    PetscCall(PetscBLASIntCast(ml,&l));
    n = l; lda = m; ldb = m; ldc = l;
    if (!k) PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&l,&n,&m,&alpha,sarray,&lds,sarray,&lds,&beta,pA,&ldc));
    else PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&l,&n,&m,&alpha,Q1,&lda,Q1,&ldb,&beta,pA,&ldc));
    PetscCall(PetscArrayzero(temp2,ml*ml));
    PetscCall(PetscMPIIntCast(ml*ml,&len));
    PetscCallMPI(MPIU_Allreduce(pA,temp2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)S)));

    PetscCall(PetscBLASIntCast(ml,&m));
    n = m; lda = m; lwork = 5*m, ldu = 1; ldvt = 1;
#if defined(PETSC_USE_COMPLEX)
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("O","N",&m,&n,temp2,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,rwork,&info));
#else
    PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("O","N",&m,&n,temp2,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,&info));
#endif
    SlepcCheckLapackInfo("gesvd",info);

    PetscCall(PetscBLASIntCast(S->n,&l));
    PetscCall(PetscBLASIntCast(ml,&n));
    m = n; lda = l; ldb = m; ldc = l;
    if (!k) PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&m,&alpha,sarray,&lds,temp2,&ldb,&beta,Q1,&ldc));
    else PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&m,&alpha,Q1,&lda,temp2,&ldb,&beta,Q2,&ldc));

    PetscCall(PetscBLASIntCast(ml,&l));
    m = l; n = l; lda = l; ldb = m; ldc = l;
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&m,&alpha,B,&lda,temp2,&ldb,&beta,tempB,&ldc));
    for (i=0;i<ml;i++) {
      sigma[i] = PetscSqrtReal(sigma[i]);
      for (j=0;j<S->n;j++) {
        if (k%2) Q2[j+i*S->n] /= sigma[i];
        else Q1[j+i*S->n] /= sigma[i];
      }
      for (j=0;j<ml;j++) B[j+i*ml] = tempB[j+i*ml]*sigma[i];
    }
  }

  PetscCall(PetscBLASIntCast(ml,&m));
  n = m; lda = m; ldu=1; ldvt=1;
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","O",&m,&n,B,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","O",&m,&n,B,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);

  PetscCall(PetscBLASIntCast(S->n,&l));
  PetscCall(PetscBLASIntCast(ml,&n));
  m = n; lda = l; ldb = m;
  if (k%2) PetscCallBLAS("BLASgemm",BLASgemm_("N","T",&l,&n,&m,&alpha,Q1,&lda,B,&ldb,&beta,sarray,&lds));
  else PetscCallBLAS("BLASgemm",BLASgemm_("N","T",&l,&n,&m,&alpha,Q2,&lda,B,&ldb,&beta,sarray,&lds));

  PetscCall(PetscFPTrapPop());
  PetscCall(BVRestoreArray(S,&sarray));

  if (rank) {
    (*rank) = 0;
    for (i=0;i<ml;i++) {
      if (sigma[i]/PetscMax(sigma[0],1.0)>delta) (*rank)++;
    }
  }
  PetscCall(PetscFree6(temp2,Q1,Q2,B,tempB,work));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVSVDAndRank_QR(BV S,PetscReal delta,PetscScalar *pA,PetscReal *sigma,PetscInt *rankout)
{
  PetscInt       i,n,ml=S->k,rank;
  PetscBLASInt   m,lda,lwork,info;
  PetscScalar    *work;
  PetscReal      *rwork;
  Mat            A;
  Vec            v;

  PetscFunctionBegin;
  /* Compute QR factorizaton of S */
  PetscCall(BVGetSizes(S,NULL,&n,NULL));
  n = PetscMin(n,ml);
  PetscCall(BVSetActiveColumns(S,0,n));
  PetscCall(PetscArrayzero(pA,ml*n));
  PetscCall(MatCreateDense(PETSC_COMM_SELF,n,n,PETSC_DECIDE,PETSC_DECIDE,pA,&A));
  PetscCall(BVOrthogonalize(S,A));
  if (n<ml) {
    /* the rest of the factorization */
    for (i=n;i<ml;i++) {
      PetscCall(BVGetColumn(S,i,&v));
      PetscCall(BVOrthogonalizeVec(S,v,pA+i*n,NULL,NULL));
      PetscCall(BVRestoreColumn(S,i,&v));
    }
  }
  PetscCall(PetscBLASIntCast(n,&lda));
  PetscCall(PetscBLASIntCast(ml,&m));
  PetscCall(PetscMalloc2(5*ml,&work,5*ml,&rwork));
  lwork = 5*m;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined (PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("O","N",&lda,&m,pA,&lda,sigma,NULL,&lda,NULL,&lda,work,&lwork,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("O","N",&lda,&m,pA,&lda,sigma,NULL,&lda,NULL,&lda,work,&lwork,rwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  PetscCall(PetscFPTrapPop());
  rank = 0;
  for (i=0;i<n;i++) {
    if (sigma[i]/PetscMax(sigma[0],1)>delta) rank++;
  }
  if (rankout) *rankout = rank;
  /* n first columns of A have the left singular vectors */
  PetscCall(BVMultInPlace(S,A,0,rank));
  PetscCall(PetscFree2(work,rwork));
  PetscCall(MatDestroy(&A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BVSVDAndRank_QR_CAA(BV S,PetscInt M,PetscInt L,PetscReal delta,PetscScalar *pA,PetscReal *sigma,PetscInt *rankout)
{
  PetscInt       i,j,n,ml=S->k,rank;
  PetscBLASInt   m,k_,lda,lwork,info;
  PetscScalar    *work,*T,*U,*R,sone=1.0,zero=0.0;
  PetscReal      *rwork;
  Mat            A;

  PetscFunctionBegin;
  /* Compute QR factorizaton of S */
  PetscCall(BVGetSizes(S,NULL,&n,NULL));
  PetscCheck(n>=ml,PetscObjectComm((PetscObject)S),PETSC_ERR_SUP,"The QR_CAA method does not support problem size n < m*L");
  PetscCall(BVSetActiveColumns(S,0,ml));
  PetscCall(PetscArrayzero(pA,ml*ml));
  PetscCall(MatCreateDense(PETSC_COMM_SELF,ml,ml,PETSC_DECIDE,PETSC_DECIDE,pA,&A));
  PetscCall(BVOrthogonalize(S,A));
  PetscCall(MatDestroy(&A));

  /* SVD of first (M-1)*L diagonal block */
  PetscCall(PetscBLASIntCast((M-1)*L,&m));
  PetscCall(PetscMalloc5(m*m,&T,m*m,&R,m*m,&U,5*ml,&work,5*ml,&rwork));
  for (j=0;j<m;j++) PetscCall(PetscArraycpy(R+j*m,pA+j*ml,m));
  lwork = 5*m;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined (PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&m,&m,R,&m,sigma,U,&m,NULL,&m,work,&lwork,&info));
#else
  PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&m,&m,R,&m,sigma,U,&m,NULL,&m,work,&lwork,rwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  PetscCall(PetscFPTrapPop());
  rank = 0;
  for (i=0;i<m;i++) {
    if (sigma[i]/PetscMax(sigma[0],1)>delta) rank++;
  }
  if (rankout) *rankout = rank;
  PetscCall(MatCreateDense(PETSC_COMM_SELF,m,m,PETSC_DECIDE,PETSC_DECIDE,U,&A));
  PetscCall(BVSetActiveColumns(S,0,m));
  PetscCall(BVMultInPlace(S,A,0,rank));
  PetscCall(MatDestroy(&A));
  /* Projected linear system */
  /* m first columns of A have the right singular vectors */
  PetscCall(PetscBLASIntCast(rank,&k_));
  PetscCall(PetscBLASIntCast(ml,&lda));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","C",&m,&k_,&m,&sone,pA+L*lda,&lda,R,&m,&zero,T,&m));
  PetscCall(PetscArrayzero(pA,ml*ml));
  PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&m,&sone,U,&m,T,&m,&zero,pA,&k_));
  for (j=0;j<k_;j++) for (i=0;i<k_;i++) pA[j*k_+i] /= sigma[j];
  PetscCall(PetscFree5(T,R,U,work,rwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSVDAndRank - Compute the SVD (left singular vectors only, and singular
   values) and determine the numerical rank according to a tolerance.

   Collective

   Input Parameters:
+  S     - the basis vectors object
.  m     - the moment degree
.  l     - the block size
.  delta - the tolerance used to determine the rank
-  meth  - the method to be used

   Output Parameters:
+  A     - workspace, on output contains relevant values in the CAA method
.  sigma - computed singular values
-  rank  - estimated rank (optional, pass `NULL` if not needed)

   Notes:
   This function computes the SVD $S=U\Sigma V^*$ and replaces $S$ with $U$.
   The current implementation computes this via $S^*S$, and it may include
   some kind of iterative refinement to improve accuracy in some cases.

   The parameters `m` and `l` refer to the moment and block size of contour
   integral methods. All columns up to `m*l` are modified, and the active
   columns are set to `0..m*l`.

   See `BVSVDMethod` for available methods.

   The `A` workspace should be `m*l*m*l` in size.

   Once the decomposition is computed, the numerical rank is estimated
   by counting the number of singular values that are larger than the
   tolerance `delta`, relative to the first singular value.

   Level: developer

.seealso: [](sec:bv), `BVSetActiveColumns()`
@*/
PetscErrorCode BVSVDAndRank(BV S,PetscInt m,PetscInt l,PetscReal delta,BVSVDMethod meth,PetscScalar A[],PetscReal sigma[],PetscInt *rank)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(S,m,2);
  PetscValidLogicalCollectiveInt(S,l,3);
  PetscValidLogicalCollectiveReal(S,delta,4);
  PetscValidLogicalCollectiveEnum(S,meth,5);
  PetscAssertPointer(A,6);
  PetscAssertPointer(sigma,7);

  PetscCall(PetscLogEventBegin(BV_SVDAndRank,S,0,0,0));
  PetscCall(BVSetActiveColumns(S,0,m*l));
  switch (meth) {
    case BV_SVD_METHOD_REFINE:
      PetscCall(BVSVDAndRank_Refine(S,delta,A,sigma,rank));
      break;
    case BV_SVD_METHOD_QR:
      PetscCall(BVSVDAndRank_QR(S,delta,A,sigma,rank));
      break;
    case BV_SVD_METHOD_QR_CAA:
      PetscCall(BVSVDAndRank_QR_CAA(S,m,l,delta,A,sigma,rank));
      break;
  }
  PetscCall(PetscLogEventEnd(BV_SVDAndRank,S,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCISSResizeBases - Resize the bases involved in CISS solvers when the $L$ parameter grows.

   Logically Collective

   Input Parameters:
+  S       - basis of $L\cdot M$ columns
.  V       - basis of $L$ columns (may be associated to subcommunicators)
.  Y       - basis of `npoints`$\cdot L$ columns
.  Lold    - old value of $L$
.  Lnew    - new value of $L$
.  M       - the moment size
-  npoints - number of integration points

   Level: developer

.seealso: [](sec:bv), `BVResize()`
@*/
PetscErrorCode BVCISSResizeBases(BV S,BV V,BV Y,PetscInt Lold,PetscInt Lnew,PetscInt M,PetscInt npoints)
{
  PetscInt       i,j;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,BV_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidLogicalCollectiveInt(S,Lold,4);
  PetscValidLogicalCollectiveInt(S,Lnew,5);
  PetscValidLogicalCollectiveInt(S,M,6);
  PetscValidLogicalCollectiveInt(S,npoints,7);

  PetscCall(BVResize(S,Lnew*M,PETSC_FALSE));
  PetscCall(BVResize(V,Lnew,PETSC_TRUE));
  PetscCall(BVResize(Y,Lnew*npoints,PETSC_TRUE));
  /* columns of Y are interleaved */
  for (i=npoints-1;i>=0;i--) {
    for (j=Lold-1;j>=0;j--) PetscCall(BVCopyColumn(Y,i*Lold+j,i*Lnew+j));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
