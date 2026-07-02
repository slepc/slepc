/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepccontour.h>
#include <slepcblaslapack.h>

/*
   SlepcContourDataCreate - Create a contour data structure.

   Input Parameters:
   n - the number of integration points
   npart - number of partitions for the subcommunicator
   parent - parent object
*/
PetscErrorCode SlepcContourDataCreate(PetscInt n,PetscInt npart,PetscObject parent,SlepcContourData *contour)
{
  PetscFunctionBegin;
  PetscCall(PetscNew(contour));
  (*contour)->parent = parent;
  PetscCall(PetscSubcommCreate(PetscObjectComm(parent),&(*contour)->subcomm));
  PetscCall(PetscSubcommSetNumber((*contour)->subcomm,npart));
  PetscCall(PetscSubcommSetType((*contour)->subcomm,PETSC_SUBCOMM_INTERLACED));
  (*contour)->npoints = n / npart;
  if (n%npart > (*contour)->subcomm->color) (*contour)->npoints++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcContourDataReset - Resets the KSP objects in a contour data structure,
   and destroys any objects whose size depends on the problem size.
*/
PetscErrorCode SlepcContourDataReset(SlepcContourData contour)
{
  PetscInt       i;

  PetscFunctionBegin;
  if (contour->ksp) {
    for (i=0;i<contour->npoints;i++) PetscCall(KSPReset(contour->ksp[i]));
  }
  if (contour->pA) {
    PetscCall(MatDestroyMatrices(contour->nmat,&contour->pA));
    PetscCall(MatDestroyMatrices(contour->nmat,&contour->pP));
    contour->nmat = 0;
  }
  PetscCall(VecScatterDestroy(&contour->scatterin));
  PetscCall(VecDestroy(&contour->xsub));
  PetscCall(VecDestroy(&contour->xdup));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcContourDataDestroy - Destroys the contour data structure.
*/
PetscErrorCode SlepcContourDataDestroy(SlepcContourData *contour)
{
  PetscInt       i;

  PetscFunctionBegin;
  if (!(*contour)) PetscFunctionReturn(PETSC_SUCCESS);
  if ((*contour)->ksp) {
    for (i=0;i<(*contour)->npoints;i++) PetscCall(KSPDestroy(&(*contour)->ksp[i]));
    PetscCall(PetscFree((*contour)->ksp));
  }
  PetscCall(PetscSubcommDestroy(&(*contour)->subcomm));
  PetscCall(PetscFree((*contour)));
  *contour = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcContourRedundantMat - Creates redundant copies of the passed matrices in the subcomm.

   Input Parameters:
   nmat - the number of matrices
   A    - array of matrices
   P    - array of matrices (preconditioner)
*/
PetscErrorCode SlepcContourRedundantMat(SlepcContourData contour,PetscInt nmat,Mat *A,Mat *P)
{
  PetscInt       i;
  MPI_Comm       child;

  PetscFunctionBegin;
  if (contour->pA) {
    PetscCall(MatDestroyMatrices(contour->nmat,&contour->pA));
    PetscCall(MatDestroyMatrices(contour->nmat,&contour->pP));
    contour->nmat = 0;
  }
  if (contour->subcomm && contour->subcomm->n != 1) {
    PetscCall(PetscSubcommGetChild(contour->subcomm,&child));
    PetscCall(PetscCalloc1(nmat,&contour->pA));
    for (i=0;i<nmat;i++) PetscCall(MatCreateRedundantMatrix(A[i],contour->subcomm->n,child,MAT_INITIAL_MATRIX,&contour->pA[i]));
    if (P) {
      PetscCall(PetscCalloc1(nmat,&contour->pP));
      for (i=0;i<nmat;i++) PetscCall(MatCreateRedundantMatrix(P[i],contour->subcomm->n,child,MAT_INITIAL_MATRIX,&contour->pP[i]));
    }
    contour->nmat = nmat;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcContourScatterCreate - Creates a scatter context to communicate between a
   regular vector and a vector xdup that can hold one duplicate per each subcommunicator
   on the contiguous parent communicator. Also creates auxiliary vectors xdup and xsub
   (the latter with the same layout as the redundant matrices in the subcommunicator).

   Input Parameters:
   v - the regular vector from which dimensions are taken
*/
PetscErrorCode SlepcContourScatterCreate(SlepcContourData contour,Vec v)
{
  IS             is1,is2;
  PetscInt       i,j,k,m,mstart,mend,mlocal;
  PetscInt       *idx1,*idx2,mloc_sub;
  MPI_Comm       contpar,parent;

  PetscFunctionBegin;
  PetscCall(VecDestroy(&contour->xsub));
  PetscCall(MatCreateVecsEmpty(contour->pA[0],&contour->xsub,NULL));

  PetscCall(VecDestroy(&contour->xdup));
  PetscCall(MatGetLocalSize(contour->pA[0],&mloc_sub,NULL));
  PetscCall(PetscSubcommGetContiguousParent(contour->subcomm,&contpar));
  PetscCall(VecCreate(contpar,&contour->xdup));
  PetscCall(VecSetSizes(contour->xdup,mloc_sub,PETSC_DECIDE));
  PetscCall(VecSetType(contour->xdup,((PetscObject)v)->type_name));

  PetscCall(VecScatterDestroy(&contour->scatterin));
  PetscCall(VecGetSize(v,&m));
  PetscCall(VecGetOwnershipRange(v,&mstart,&mend));
  mlocal = mend - mstart;
  PetscCall(PetscMalloc2(contour->subcomm->n*mlocal,&idx1,contour->subcomm->n*mlocal,&idx2));
  j = 0;
  for (k=0;k<contour->subcomm->n;k++) {
    for (i=mstart;i<mend;i++) {
      idx1[j]   = i;
      idx2[j++] = i + m*k;
    }
  }
  PetscCall(PetscSubcommGetParent(contour->subcomm,&parent));
  PetscCall(ISCreateGeneral(parent,contour->subcomm->n*mlocal,idx1,PETSC_COPY_VALUES,&is1));
  PetscCall(ISCreateGeneral(parent,contour->subcomm->n*mlocal,idx2,PETSC_COPY_VALUES,&is2));
  PetscCall(VecScatterCreate(v,is1,contour->xdup,is2,&contour->scatterin));
  PetscCall(ISDestroy(&is1));
  PetscCall(ISDestroy(&is2));
  PetscCall(PetscFree2(idx1,idx2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcCISS_isGhost - Determine if any of the computed eigenpairs are spurious.

   Input Parameters:
   X - the matrix of eigenvectors (MATSEQDENSE)
   n - the number of columns to consider
   sigma - the singular values
   thresh - threshold to decide whether a value is spurious

   Output Parameter:
   fl - array of n booleans
*/
PetscErrorCode SlepcCISS_isGhost(Mat X,PetscInt n,PetscReal *sigma,PetscReal thresh,PetscBool *fl)
{
  const PetscScalar *pX;
  PetscInt          i,j,m,ld;
  PetscReal         *tau,s1,s2,tau_max=0.0;

  PetscFunctionBegin;
  PetscCall(MatGetSize(X,&m,NULL));
  PetscCall(MatDenseGetLDA(X,&ld));
  PetscCall(PetscMalloc1(n,&tau));
  PetscCall(MatDenseGetArrayRead(X,&pX));
  for (j=0;j<n;j++) {
    s1 = 0.0;
    s2 = 0.0;
    for (i=0;i<m;i++) {
      s1 += PetscAbsScalar(PetscPowScalarInt(pX[i+j*ld],2));
      s2 += PetscPowRealInt(PetscAbsScalar(pX[i+j*ld]),2)/sigma[i];
    }
    tau[j] = s1/s2;
    tau_max = PetscMax(tau_max,tau[j]);
  }
  PetscCall(MatDenseRestoreArrayRead(X,&pX));
  for (j=0;j<n;j++) fl[j] = (tau[j]>=thresh*tau_max)? PETSC_TRUE: PETSC_FALSE;
  PetscCall(PetscFree(tau));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcCISS_BH_SVD - Compute SVD of block Hankel matrix and its rank.

   Input Parameters:
   H  - block Hankel matrix obtained via CISS_BlockHankel()
   ml - dimension of rows and columns, equal to M*L
   delta - the tolerance used to determine the rank

   Output Parameters:
   sigma - computed singular values
   rank  - the rank of H
*/
PetscErrorCode SlepcCISS_BH_SVD(PetscScalar *H,PetscInt ml,PetscReal delta,PetscReal *sigma,PetscInt *rank)
{
  PetscInt       i;
  PetscBLASInt   m,n,lda,ldu,ldvt,lwork;
  PetscScalar    *work;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(5*ml,&work));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc1(5*ml,&rwork));
#endif
  PetscCall(PetscBLASIntCast(ml,&m));
  n = m; lda = m; ldu = m; ldvt = m; lwork = 5*m;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if defined(PETSC_USE_COMPLEX)
  PetscCallLAPACKInfo("LAPACKgesvd",LAPACKgesvd_("N","N",&m,&n,H,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,rwork,&info));
#else
  PetscCallLAPACKInfo("LAPACKgesvd",LAPACKgesvd_("N","N",&m,&n,H,&lda,sigma,NULL,&ldu,NULL,&ldvt,work,&lwork,&info));
#endif
  PetscCall(PetscFPTrapPop());
  (*rank) = 0;
  for (i=0;i<ml;i++) {
    if (sigma[i]/PetscMax(sigma[0],1.0)>delta) (*rank)++;
  }
  PetscCall(PetscFree(work));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}
