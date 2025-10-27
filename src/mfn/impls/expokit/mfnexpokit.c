/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix function solver: "expokit"

   Method: Arnoldi method tailored for the matrix exponential

   Algorithm:

       Uses Arnoldi relations to compute exp(t_step*A)*v_last for
       several time steps.

   References:

       [1] R. Sidje, "Expokit: a software package for computing matrix
           exponentials", ACM Trans. Math. Softw. 24(1):130-156, 1998.
*/

#include <slepc/private/mfnimpl.h>

static PetscErrorCode MFNSetUp_Expokit(MFN mfn)
{
  PetscInt       N;
  PetscBool      isexp;

  PetscFunctionBegin;
  PetscCall(MatGetSize(mfn->A,&N,NULL));
  if (mfn->ncv==PETSC_DETERMINE) mfn->ncv = PetscMin(30,N);
  if (mfn->max_it==PETSC_DETERMINE) mfn->max_it = 100;
  PetscCall(MFNAllocateSolution(mfn,2));

  PetscCall(PetscObjectTypeCompare((PetscObject)mfn->fn,FNEXP,&isexp));
  PetscCheck(isexp,PETSC_COMM_SELF,PETSC_ERR_SUP,"This solver only supports the exponential function");
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFNSolve_Expokit(MFN mfn,Vec b,Vec x)
{
  PetscInt          mxstep,mxrej,m,mb,ld,i,j,ireject,mx,k1;
  Vec               v,r;
  Mat               H,M=NULL,K=NULL;
  FN                fn;
  PetscScalar       *Harray,*B,*F,*betaF,t,sgn,sfactor;
  const PetscScalar *pK;
  PetscReal         anorm,avnorm,tol,err_loc,rndoff,t_out,t_new,t_now,t_step;
  PetscReal         xm,fact,s,p1,p2,beta,beta2,gamma,delta;
  PetscBool         breakdown;

  PetscFunctionBegin;
  m   = mfn->ncv;
  tol = mfn->tol;
  mxstep = mfn->max_it;
  mxrej = 10;
  gamma = 0.9;
  delta = 1.2;
  mb    = m;
  PetscCall(FNGetScale(mfn->fn,&t,&sfactor));
  PetscCall(FNDuplicate(mfn->fn,PetscObjectComm((PetscObject)mfn->fn),&fn));
  PetscCall(FNSetScale(fn,1.0,1.0));
  t_out = PetscAbsScalar(t);
  t_now = 0.0;
  PetscCall(MatNorm(mfn->A,NORM_INFINITY,&anorm));
  rndoff = anorm*PETSC_MACHINE_EPSILON;

  k1 = 2;
  xm = 1.0/(PetscReal)m;
  beta = mfn->bnorm;
  fact = PetscPowRealInt((m+1)/2.72,m+1)*PetscSqrtReal(2.0*PETSC_PI*(m+1));
  t_new = (1.0/anorm)*PetscPowReal((fact*tol)/(4.0*beta*anorm),xm);
  s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_new))-1);
  t_new = PetscCeilReal(t_new/s)*s;
  sgn = t/PetscAbsScalar(t);

  PetscCall(VecCopy(b,x));
  ld = m+2;
  PetscCall(PetscCalloc2(m+1,&betaF,ld*ld,&B));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ld,ld,NULL,&H));
  PetscCall(MatDenseGetArray(H,&Harray));

  while (mfn->reason == MFN_CONVERGED_ITERATING) {
    mfn->its++;
    if (PetscIsInfOrNanReal(t_new)) t_new = PETSC_MAX_REAL;
    t_step = PetscMin(t_out-t_now,t_new);
    PetscCall(BVInsertVec(mfn->V,0,x));
    PetscCall(BVScaleColumn(mfn->V,0,1.0/beta));
    PetscCall(BVMatArnoldi(mfn->V,mfn->transpose_solve?mfn->AT:mfn->A,H,0,&mb,&beta2,&breakdown));
    if (breakdown) {
      k1 = 0;
      t_step = t_out-t_now;
    }
    if (k1!=0) {
      Harray[m+1+ld*m] = 1.0;
      PetscCall(BVGetColumn(mfn->V,m,&v));
      PetscCall(BVGetColumn(mfn->V,m+1,&r));
      PetscCall(MatMult(mfn->transpose_solve?mfn->AT:mfn->A,v,r));
      PetscCall(BVRestoreColumn(mfn->V,m,&v));
      PetscCall(BVRestoreColumn(mfn->V,m+1,&r));
      PetscCall(BVNormColumn(mfn->V,m+1,NORM_2,&avnorm));
    }
    PetscCall(PetscArraycpy(B,Harray,ld*ld));

    ireject = 0;
    while (ireject <= mxrej) {
      mx = mb + k1;
      for (i=0;i<mx;i++) {
        for (j=0;j<mx;j++) {
          Harray[i+j*ld] = sgn*B[i+j*ld]*t_step;
        }
      }
      PetscCall(MFN_CreateDenseMat(mx,&M));
      PetscCall(MFN_CreateDenseMat(mx,&K));
      PetscCall(MatDenseGetArray(M,&F));
      for (j=0;j<mx;j++) PetscCall(PetscArraycpy(F+j*mx,Harray+j*ld,mx));
      PetscCall(MatDenseRestoreArray(M,&F));
      PetscCall(FNEvaluateFunctionMat(fn,M,K));

      if (k1==0) {
        err_loc = tol;
        break;
      } else {
        PetscCall(MatDenseGetArrayRead(K,&pK));
        p1 = PetscAbsScalar(beta*pK[m]);
        p2 = PetscAbsScalar(beta*pK[m+1]*avnorm);
        PetscCall(MatDenseRestoreArrayRead(K,&pK));
        if (p1 > 10*p2) {
          err_loc = p2;
          xm = 1.0/(PetscReal)m;
        } else if (p1 > p2) {
          err_loc = (p1*p2)/(p1-p2);
          xm = 1.0/(PetscReal)m;
        } else {
          err_loc = p1;
          xm = 1.0/(PetscReal)(m-1);
        }
      }
      if (err_loc <= delta*t_step*tol) break;
      else {
        t_step = gamma*t_step*PetscPowReal(t_step*tol/err_loc,xm);
        s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_step))-1);
        t_step = PetscCeilReal(t_step/s)*s;
        ireject = ireject+1;
      }
    }

    mx = mb + PetscMax(0,k1-1);
    PetscCall(MatDenseGetArrayRead(K,&pK));
    for (j=0;j<mx;j++) betaF[j] = beta*pK[j];
    PetscCall(MatDenseRestoreArrayRead(K,&pK));
    PetscCall(BVSetActiveColumns(mfn->V,0,mx));
    PetscCall(BVMultVec(mfn->V,1.0,0.0,x,betaF));
    PetscCall(VecNorm(x,NORM_2,&beta));

    t_now = t_now+t_step;
    if (t_now>=t_out) mfn->reason = MFN_CONVERGED_TOL;
    else {
      t_new = gamma*t_step*PetscPowReal((t_step*tol)/err_loc,xm);
      s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_new))-1);
      t_new = PetscCeilReal(t_new/s)*s;
    }
    err_loc = PetscMax(err_loc,rndoff);
    if (mfn->its==mxstep) mfn->reason = MFN_DIVERGED_ITS;
    PetscCall(MFNMonitor(mfn,mfn->its,err_loc));
  }
  PetscCall(VecScale(x,sfactor));

  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&K));
  PetscCall(FNDestroy(&fn));
  PetscCall(MatDenseRestoreArray(H,&Harray));
  PetscCall(MatDestroy(&H));
  PetscCall(PetscFree2(betaF,B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   MFNEXPOKIT - MFNEXPOKIT = "expokit" - An Arnoldi method tailored for the
   matrix exponential, as implemented in the Expokit package.

   Notes:
   This solver implements the method of Expokit natively in SLEPc, rather than
   interfacing to the original Expokit software {cite:p}`Sid98`.

   The method uses Arnoldi relations to compute $\exp(\alpha A)v_\mathrm{last}$
   for several time steps, where $\alpha$ is the time step size.

   Level: beginner

.seealso: [](ch:mfn), `MFN`, `MFNType`, `MFNSetType()`
M*/
SLEPC_EXTERN PetscErrorCode MFNCreate_Expokit(MFN mfn)
{
  PetscFunctionBegin;
  mfn->ops->solve          = MFNSolve_Expokit;
  mfn->ops->setup          = MFNSetUp_Expokit;
  PetscFunctionReturn(PETSC_SUCCESS);
}
