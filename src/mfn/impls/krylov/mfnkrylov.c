/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix function solver: "krylov"

   Method: Arnoldi with Eiermann-Ernst restart

   Algorithm:

       Build Arnoldi approximations using f(H) for the Hessenberg matrix H,
       restart by discarding the Krylov basis but keeping H.

   References:

       [1] M. Eiermann and O. Ernst, "A restarted Krylov subspace method
           for the evaluation of matrix functions", SIAM J. Numer. Anal.
           44(6):2481-2504, 2006.
*/

#include <slepc/private/mfnimpl.h>
#include <slepcblaslapack.h>

static PetscErrorCode MFNSetUp_Krylov(MFN mfn)
{
  PetscInt       N;

  PetscFunctionBegin;
  PetscCall(MatGetSize(mfn->A,&N,NULL));
  if (mfn->ncv==PETSC_DETERMINE) mfn->ncv = PetscMin(30,N);
  if (mfn->max_it==PETSC_DETERMINE) mfn->max_it = 100;
  PetscCall(MFNAllocateSolution(mfn,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNSolve_Krylov(MFN mfn,Vec b,Vec x)
{
  PetscInt          n=0,m,ld,ldh,j;
  PetscBLASInt      m_,inc=1;
  Mat               M,G=NULL,H=NULL;
  Vec               F=NULL;
  PetscScalar       *marray,*farray,*harray;
  const PetscScalar *garray;
  PetscReal         beta,betaold=0.0,nrm=1.0;
  PetscBool         breakdown;

  PetscFunctionBegin;
  m  = mfn->ncv;
  ld = m+1;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ld,m,NULL,&M));
  PetscCall(MatDenseGetArray(M,&marray));

  /* set initial vector to b/||b|| */
  PetscCall(BVInsertVec(mfn->V,0,b));
  PetscCall(BVScaleColumn(mfn->V,0,1.0/mfn->bnorm));
  PetscCall(VecSet(x,0.0));

  /* Restart loop */
  while (mfn->reason == MFN_CONVERGED_ITERATING) {
    mfn->its++;

    /* compute Arnoldi factorization */
    PetscCall(BVMatArnoldi(mfn->V,mfn->transpose_solve?mfn->AT:mfn->A,M,0,&m,&beta,&breakdown));

    /* save previous Hessenberg matrix in G; allocate new storage for H and f(H) */
    if (mfn->its>1) { G = H; H = NULL; }
    ldh = n+m;
    PetscCall(MFN_CreateVec(ldh,&F));
    PetscCall(MFN_CreateDenseMat(ldh,&H));

    /* glue together the previous H and the new H obtained with Arnoldi */
    PetscCall(MatDenseGetArray(H,&harray));
    for (j=0;j<m;j++) PetscCall(PetscArraycpy(harray+n+(j+n)*ldh,marray+j*ld,m));
    if (mfn->its>1) {
      PetscCall(MatDenseGetArrayRead(G,&garray));
      for (j=0;j<n;j++) PetscCall(PetscArraycpy(harray+j*ldh,garray+j*n,n));
      PetscCall(MatDenseRestoreArrayRead(G,&garray));
      PetscCall(MatDestroy(&G));
      harray[n+(n-1)*ldh] = betaold;
    }
    PetscCall(MatDenseRestoreArray(H,&harray));

    if (mfn->its==1) {
      /* set symmetry flag of H from A */
      PetscCall(MatPropagateSymmetryOptions(mfn->A,H));
    }

    /* evaluate f(H) */
    PetscCall(FNEvaluateFunctionMatVec(mfn->fn,H,F));

    /* x += ||b||*V*f(H)*e_1 */
    PetscCall(VecGetArray(F,&farray));
    PetscCall(PetscBLASIntCast(m,&m_));
    nrm = BLASnrm2_(&m_,farray+n,&inc);   /* relative norm of the update ||u||/||b|| */
    PetscCall(MFNMonitor(mfn,mfn->its,nrm));
    for (j=0;j<m;j++) farray[j+n] *= mfn->bnorm;
    PetscCall(BVSetActiveColumns(mfn->V,0,m));
    PetscCall(BVMultVec(mfn->V,1.0,1.0,x,farray+n));
    PetscCall(VecRestoreArray(F,&farray));

    /* check convergence */
    if (mfn->its >= mfn->max_it) mfn->reason = MFN_DIVERGED_ITS;
    if (mfn->its>1) {
      if (m<mfn->ncv || breakdown || beta==0.0 || nrm<mfn->tol) mfn->reason = MFN_CONVERGED_TOL;
    }

    /* restart with vector v_{m+1} */
    if (mfn->reason == MFN_CONVERGED_ITERATING) {
      PetscCall(BVCopyColumn(mfn->V,m,0));
      n += m;
      betaold = beta;
    }
  }

  PetscCall(MatDestroy(&H));
  PetscCall(MatDestroy(&G));
  PetscCall(VecDestroy(&F));
  PetscCall(MatDenseRestoreArray(M,&marray));
  PetscCall(MatDestroy(&M));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   MFNKRYLOV - MFNKRYLOV = "krylov" - A restarted Krylov method for the application
   of a matrix function to a vector.

   Note:
   This solver builds Arnoldi approximations using $f(H)$ for the Hessenberg matrix
   $H$, and restarts by discarding the Krylov basis but keeping $H$, as proposed by
   {cite:t}`Eie06`.

   Level: beginner

.seealso: [](ch:mfn), `MFN`, `MFNType`, `MFNSetType()`
M*/
SLEPC_EXTERN PetscErrorCode MFNCreate_Krylov(MFN mfn)
{
  PetscFunctionBegin;
  mfn->ops->solve          = MFNSolve_Krylov;
  mfn->ops->setup          = MFNSetUp_Krylov;
  PetscFunctionReturn(PETSC_SUCCESS);
}
