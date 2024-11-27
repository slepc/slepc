/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "power"

   Method: Power Iteration

   Algorithm:

       This solver implements the power iteration for finding dominant
       eigenpairs. It also includes the following well-known methods:
       - Inverse Iteration: when used in combination with shift-and-invert
         spectral transformation.
       - Rayleigh Quotient Iteration (RQI): also with shift-and-invert plus
         a variable shift.

       It can also be used for nonlinear inverse iteration on the problem
       A(x)*x=lambda*B(x)*x, where A and B are not constant but depend on x.

   References:

       [1] "Single Vector Iteration Methods in SLEPc", SLEPc Technical Report
           STR-2, available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <slepcblaslapack.h>
/* petsc headers */
#include <petscdm.h>
#include <petscsnes.h>

static PetscErrorCode EPSPowerFormFunction_Update(SNES,Vec,Vec,void*);
static PetscErrorCode EPSSolve_Power(EPS);
static PetscErrorCode EPSSolve_TS_Power(EPS);

typedef struct {
  EPSPowerShiftType shift_type;
  PetscBool         nonlinear;
  PetscBool         update;
  PetscBool         sign_normalization;
  SNES              snes;
  PetscErrorCode    (*formFunctionB)(SNES,Vec,Vec,void*);
  void              *formFunctionBctx;
  PetscErrorCode    (*formFunctionA)(SNES,Vec,Vec,void*);
  void              *formFunctionActx;
  PetscErrorCode    (*formFunctionAB)(SNES,Vec,Vec,Vec,void*);
  PetscErrorCode    (*formNorm)(SNES,Vec,PetscReal*,void*);
  void              *formNormCtx;
  PetscInt          idx;  /* index of the first nonzero entry in the iteration vector */
  PetscMPIInt       p;    /* process id of the owner of idx */
  PetscReal         norm0; /* norm of initial vector */
} EPS_POWER;

static PetscErrorCode SNESMonitor_PowerUpdate(SNES snes,PetscInt its,PetscReal fnorm,void *ctx)
{
  EPS            eps;

  PetscFunctionBegin;
  PetscCall(PetscObjectQuery((PetscObject)snes,"eps",(PetscObject *)&eps));
  PetscCheck(eps,PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_NULL,"No composed EPS");
  /* Call EPS monitor on each SNES iteration */
  PetscCall(EPSMonitor(eps,its,eps->nconv,eps->eigr,eps->eigi,eps->errest,PetscMin(eps->nconv+1,eps->nev)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  STMatMode      mode;
  Mat            A,B,P;
  Vec            res;
  PetscContainer container;
  PetscErrorCode (*formFunctionA)(SNES,Vec,Vec,void*);
  PetscErrorCode (*formJacobianA)(SNES,Vec,Mat,Mat,void*);
  PetscErrorCode (*formNorm)(SNES,Vec,PetscReal*,void*);
  void           *ctx;

  PetscFunctionBegin;
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev");
  } else eps->ncv = eps->nev;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) {
    /* SNES will directly return the solution for us, and we need to do only one iteration */
    if (power->nonlinear && power->update) eps->max_it = 1;
    else eps->max_it = PetscMax(1000*eps->nev,100*eps->n);
  }
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which==EPS_LARGEST_MAGNITUDE || eps->which==EPS_TARGET_MAGNITUDE,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only largest magnitude or target magnitude eigenvalues");
  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) {
    PetscCheck(!power->nonlinear,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Variable shifts not allowed in nonlinear problems");
    EPSCheckSinvertCayleyCondition(eps,PETSC_TRUE," (with variable shifts)");
    PetscCall(STGetMatMode(eps->st,&mode));
    PetscCheck(mode!=ST_MATMODE_INPLACE,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"ST matrix mode inplace does not work with variable shifts");
  }
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_THRESHOLD);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION);
  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(EPS_SetInnerProduct(eps));

  if (power->nonlinear) {
    PetscCheck(eps->nev==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Nonlinear inverse iteration cannot compute more than one eigenvalue");
    PetscCall(EPSSetWorkVecs(eps,3));
    PetscCheck(!power->update || eps->max_it==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"More than one iteration is not allowed for Newton eigensolver (SNES)");

    /* Set up SNES solver */
    /* If SNES was setup, we need to reset it so that it will be consistent with the current EPS */
    if (power->snes) PetscCall(SNESReset(power->snes));
    else PetscCall(EPSPowerGetSNES(eps,&power->snes));

    /* For nonlinear solver (Newton), we should scale the initial vector back.
       The initial vector will be scaled in EPSSetUp. */
    if (eps->IS) PetscCall(VecNorm(eps->IS[0],NORM_2,&power->norm0));

    PetscCall(EPSGetOperators(eps,&A,&B));

    /* form function */
    PetscCall(PetscObjectQueryFunction((PetscObject)A,"formFunction",&formFunctionA));
    PetscCheck(formFunctionA,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER,"For nonlinear inverse iteration you must compose a callback function 'formFunction' in matrix A");
    PetscCall(PetscObjectQuery((PetscObject)A,"formFunctionCtx",(PetscObject*)&container));
    if (container) PetscCall(PetscContainerGetPointer(container,&ctx));
    else ctx = NULL;
    PetscCall(MatCreateVecs(A,&res,NULL));
    power->formFunctionA = formFunctionA;
    power->formFunctionActx = ctx;
    if (power->update) {
      PetscCall(SNESSetFunction(power->snes,res,EPSPowerFormFunction_Update,ctx));
      PetscCall(PetscObjectQueryFunction((PetscObject)A,"formFunctionAB",&power->formFunctionAB));
      PetscCall(SNESMonitorSet(power->snes,SNESMonitor_PowerUpdate,NULL,NULL));
    }
    else PetscCall(SNESSetFunction(power->snes,res,formFunctionA,ctx));
    PetscCall(VecDestroy(&res));

    /* form Jacobian */
    PetscCall(PetscObjectQueryFunction((PetscObject)A,"formJacobian",&formJacobianA));
    PetscCheck(formJacobianA,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER,"For nonlinear inverse iteration you must compose a callback function 'formJacobian' in matrix A");
    PetscCall(PetscObjectQuery((PetscObject)A,"formJacobianCtx",(PetscObject*)&container));
    if (container) PetscCall(PetscContainerGetPointer(container,&ctx));
    else ctx = NULL;
    /* This allows users to compute a different preconditioning matrix than A.
     * It is useful when A and B are shell matrices.
     */
    PetscCall(STGetPreconditionerMat(eps->st,&P));
    /* If users do not provide a matrix, we simply use A */
    PetscCall(SNESSetJacobian(power->snes,A,P? P:A,formJacobianA,ctx));
    PetscCall(SNESSetFromOptions(power->snes));
    PetscCall(SNESSetUp(power->snes));

    ctx = NULL;
    formNorm = NULL;
    if (B) {
      PetscCall(PetscObjectQueryFunction((PetscObject)B,"formFunction",&power->formFunctionB));
      PetscCall(PetscObjectQuery((PetscObject)B,"formFunctionCtx",(PetscObject*)&container));
      if (power->formFunctionB && container) PetscCall(PetscContainerGetPointer(container,&power->formFunctionBctx));
      else power->formFunctionBctx = NULL;

      /* form norm */
      PetscCall(PetscObjectQueryFunction((PetscObject)B,"formNorm",&formNorm));
      if (formNorm) {
        PetscCall(PetscObjectQuery((PetscObject)B,"formNormCtx",(PetscObject*)&container));
        if (container) PetscCall(PetscContainerGetPointer(container,&ctx));
      }
    }
    power->formNorm = formNorm;
    power->formNormCtx = ctx;
  } else {
    if (eps->twosided) PetscCall(EPSSetWorkVecs(eps,3));
    else PetscCall(EPSSetWorkVecs(eps,2));
    PetscCall(DSSetType(eps->ds,DSNHEP));
    PetscCall(DSAllocate(eps->ds,eps->nev));
  }
  /* dispatch solve method */
  if (eps->twosided) {
    PetscCheck(!power->nonlinear,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Nonlinear inverse iteration does not have two-sided variant");
    PetscCheck(power->shift_type!=EPS_POWER_SHIFT_WILKINSON,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Two-sided variant does not support Wilkinson shifts");
    eps->ops->solve = EPSSolve_TS_Power;
  } else eps->ops->solve = EPSSolve_Power;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Find the index of the first nonzero entry of x, and the process that owns it.
*/
static PetscErrorCode FirstNonzeroIdx(Vec x,PetscInt *idx,PetscMPIInt *p)
{
  PetscInt          i,first,last,N;
  PetscLayout       map;
  const PetscScalar *xx;

  PetscFunctionBegin;
  PetscCall(VecGetSize(x,&N));
  PetscCall(VecGetOwnershipRange(x,&first,&last));
  PetscCall(VecGetArrayRead(x,&xx));
  for (i=first;i<last;i++) {
    if (PetscAbsScalar(xx[i-first])>10*PETSC_MACHINE_EPSILON) break;
  }
  PetscCall(VecRestoreArrayRead(x,&xx));
  if (i==last) i=N;
  PetscCallMPI(MPIU_Allreduce(&i,idx,1,MPIU_INT,MPI_MIN,PetscObjectComm((PetscObject)x)));
  PetscCheck(*idx!=N,PetscObjectComm((PetscObject)x),PETSC_ERR_PLIB,"Zero vector found");
  PetscCall(VecGetLayout(x,&map));
  PetscCall(PetscLayoutFindOwner(map,*idx,p));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Normalize a vector x with respect to a given norm as well as, optionally, the
   sign of the first nonzero entry (assumed to be idx in process p).
*/
static PetscErrorCode Normalize(Vec x,PetscReal norm,PetscInt idx,PetscMPIInt p,PetscBool sign_normalization,PetscScalar *sign)
{
  PetscScalar       alpha=1.0;
  PetscInt          first,last;
  PetscReal         absv;
  const PetscScalar *xx;

  PetscFunctionBegin;
  if (sign_normalization) {
    PetscCall(VecGetOwnershipRange(x,&first,&last));
    if (idx>=first && idx<last) {
      PetscCall(VecGetArrayRead(x,&xx));
      absv = PetscAbsScalar(xx[idx-first]);
      if (absv>10*PETSC_MACHINE_EPSILON) alpha = xx[idx-first]/absv;
      PetscCall(VecRestoreArrayRead(x,&xx));
    }
    PetscCallMPI(MPI_Bcast(&alpha,1,MPIU_SCALAR,p,PetscObjectComm((PetscObject)x)));
  }
  if (sign) *sign = alpha;
  alpha *= norm;
  PetscCall(VecScale(x,1.0/alpha));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerUpdateFunctionB(EPS eps,Vec x,Vec Bx)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Mat            B;

  PetscFunctionBegin;
  PetscCall(STResetMatrixState(eps->st));
  PetscCall(EPSGetOperators(eps,NULL,&B));
  if (B) {
    if (power->formFunctionB) PetscCall((*power->formFunctionB)(power->snes,x,Bx,power->formFunctionBctx));
    else PetscCall(MatMult(B,x,Bx));
  } else PetscCall(VecCopy(x,Bx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerUpdateFunctionA(EPS eps,Vec x,Vec Ax)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Mat            A;

  PetscFunctionBegin;
  PetscCall(STResetMatrixState(eps->st));
  PetscCall(EPSGetOperators(eps,&A,NULL));
  PetscCheck(A,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_NULL,"Matrix A is required for an eigenvalue problem");
  if (power->formFunctionA) PetscCall((*power->formFunctionA)(power->snes,x,Ax,power->formFunctionActx));
  else PetscCall(MatMult(A,x,Ax));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerFormFunction_Update(SNES snes,Vec x,Vec y,void *ctx)
{
  EPS            eps;
  PetscReal      bx;
  Vec            Bx;
  PetscScalar    sign;
  EPS_POWER      *power;

  PetscFunctionBegin;
  PetscCall(PetscObjectQuery((PetscObject)snes,"eps",(PetscObject *)&eps));
  PetscCheck(eps,PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_NULL,"No composed EPS");
  PetscCall(STResetMatrixState(eps->st));
  power = (EPS_POWER*)eps->data;
  Bx = eps->work[2];
  if (power->formFunctionAB) PetscCall((*power->formFunctionAB)(snes,x,y,Bx,ctx));
  else {
    PetscCall(EPSPowerUpdateFunctionA(eps,x,y));
    PetscCall(EPSPowerUpdateFunctionB(eps,x,Bx));
  }
  if (power->formNorm) PetscCall((*power->formNorm)(snes,Bx,&bx,power->formNormCtx));
  else PetscCall(VecNorm(Bx,NORM_2,&bx));
  PetscCall(Normalize(Bx,bx,power->idx,power->p,power->sign_normalization,&sign));
  PetscCall(VecAXPY(y,-1.0,Bx));
  /* Keep tracking eigenvalue update. It would be useful when we want to monitor solver progress via snes monitor. */
  eps->eigr[(eps->nconv < eps->nev)? eps->nconv:(eps->nconv-1)] = 1.0/(bx*sign);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Use SNES to compute y = A^{-1}*B*x for the nonlinear problem
*/
static PetscErrorCode EPSPowerApply_SNES(EPS eps,Vec x,Vec y)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Vec            Bx;

  PetscFunctionBegin;
  PetscCall(VecCopy(x,y));
  if (power->update) PetscCall(SNESSolve(power->snes,NULL,y));
  else {
    Bx = eps->work[2];
    PetscCall(SNESSolve(power->snes,Bx,y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Use nonlinear inverse power to compute an initial guess
*/
static PetscErrorCode EPSPowerComputeInitialGuess_Update(EPS eps)
{
  EPS            powereps;
  Mat            A,B,P;
  Vec            v1,v2;
  SNES           snes;
  DM             dm,newdm;
  PetscBool      sign_normalization;

  PetscFunctionBegin;
  PetscCall(EPSCreate(PetscObjectComm((PetscObject)eps),&powereps));
  PetscCall(EPSGetOperators(eps,&A,&B));
  PetscCall(EPSSetType(powereps,EPSPOWER));
  PetscCall(EPSSetOperators(powereps,A,B));
  PetscCall(EPSSetTolerances(powereps,1e-6,4));
  PetscCall(EPSSetOptionsPrefix(powereps,((PetscObject)eps)->prefix));
  PetscCall(EPSAppendOptionsPrefix(powereps,"init_"));
  PetscCall(EPSSetProblemType(powereps,EPS_GNHEP));
  PetscCall(EPSSetWhichEigenpairs(powereps,EPS_TARGET_MAGNITUDE));
  PetscCall(EPSPowerSetNonlinear(powereps,PETSC_TRUE));
  PetscCall(STGetPreconditionerMat(eps->st,&P));
  /* attach dm to initial solve */
  PetscCall(EPSPowerGetSNES(eps,&snes));
  PetscCall(SNESGetDM(snes,&dm));
  /* use  dmshell to temporarily store snes context */
  PetscCall(DMCreate(PetscObjectComm((PetscObject)eps),&newdm));
  PetscCall(DMSetType(newdm,DMSHELL));
  PetscCall(DMSetUp(newdm));
  PetscCall(DMCopyDMSNES(dm,newdm));
  PetscCall(EPSPowerGetSNES(powereps,&snes));
  PetscCall(SNESSetDM(snes,dm));
  PetscCall(EPSPowerGetSignNormalization(eps,&sign_normalization));
  PetscCall(EPSPowerSetSignNormalization(powereps,sign_normalization));
  PetscCall(EPSSetFromOptions(powereps));
  if (P) PetscCall(STSetPreconditionerMat(powereps->st,P));
  PetscCall(EPSSolve(powereps));
  PetscCall(BVGetColumn(eps->V,0,&v2));
  PetscCall(BVGetColumn(powereps->V,0,&v1));
  PetscCall(VecCopy(v1,v2));
  PetscCall(BVRestoreColumn(powereps->V,0,&v1));
  PetscCall(BVRestoreColumn(eps->V,0,&v2));
  PetscCall(EPSDestroy(&powereps));
  /* restore context back to the old nonlinear solver */
  PetscCall(DMCopyDMSNES(newdm,dm));
  PetscCall(DMDestroy(&newdm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_Power(EPS eps)
{
  EPS_POWER           *power = (EPS_POWER*)eps->data;
  PetscInt            k,ld;
  Vec                 v,y,e,Bx;
  Mat                 A;
  KSP                 ksp;
  PetscReal           relerr,norm,norm1,rt1,rt2,cs1;
  PetscScalar         theta,rho,delta,sigma,alpha2,beta1,sn1,*T,sign;
  PetscBool           breakdown;
  KSPConvergedReason  reason;
  SNESConvergedReason snesreason;

  PetscFunctionBegin;
  e = eps->work[0];
  y = eps->work[1];
  if (power->nonlinear) Bx = eps->work[2];
  else Bx = NULL;

  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) PetscCall(STGetKSP(eps->st,&ksp));
  if (power->nonlinear) {
    PetscCall(PetscObjectCompose((PetscObject)power->snes,"eps",(PetscObject)eps));
    /* Compute an initial guess only when users do not provide one */
    if (power->update && !eps->nini) PetscCall(EPSPowerComputeInitialGuess_Update(eps));
  } else PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  if (!power->update) PetscCall(EPSGetStartVector(eps,0,NULL));
  if (power->nonlinear) {
    PetscCall(BVGetColumn(eps->V,0,&v));
    if (eps->nini) {
      /* We scale the initial vector back if the initial vector was provided by users */
      PetscCall(VecScale(v,power->norm0));
    }
    /* Do a couple things:
     * 1) Determine the first non-zero index for Bx and the proc that owns it. This will be
     *    used if performing normalization by the sign of the first nonzero element of Bx.
     * 2) Scale Bx by its norm. For non-update power iterations, Bx (in code terms) is used
     *    as the RHS argument to SNESSolve. And recall that the formula for generalized
     *    inverse power iterations in this case is: (Ax)_n = (Bx)_{n-1}/|(Bx)_{n-1}| (in
     *    math terms)
     */
    PetscCall(EPSPowerUpdateFunctionB(eps,v,Bx));
    PetscCall(BVRestoreColumn(eps->V,0,&v));
    if (power->formNorm) PetscCall((*power->formNorm)(power->snes,Bx,&norm,power->formNormCtx));
    else PetscCall(VecNorm(Bx,NORM_2,&norm));
    PetscCall(FirstNonzeroIdx(Bx,&power->idx,&power->p));
    PetscCall(Normalize(Bx,norm,power->idx,power->p,power->sign_normalization,NULL));
  }

  PetscCall(STGetShift(eps->st,&sigma));    /* original shift */
  rho = sigma;

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    k = eps->nconv;

    /* y = OP v */
    PetscCall(BVGetColumn(eps->V,k,&v));
    if (power->nonlinear) PetscCall(EPSPowerApply_SNES(eps,v,y));
    else PetscCall(STApply(eps->st,v,y));
    PetscCall(BVRestoreColumn(eps->V,k,&v));

    /* purge previously converged eigenvectors */
    if (PetscUnlikely(power->nonlinear)) {
      /* We do not need to call this for Newton eigensolver since eigenvalue is
       * updated in function evaluations.
       */
      if (!power->update) {
        PetscCall(EPSPowerUpdateFunctionB(eps,y,Bx));
        if (power->formNorm) PetscCall((*power->formNorm)(power->snes,Bx,&norm,power->formNormCtx));
        else PetscCall(VecNorm(Bx,NORM_2,&norm));
        PetscCall(Normalize(Bx,norm,power->idx,power->p,power->sign_normalization,&sign));
      }
    } else {
      PetscCall(DSGetArray(eps->ds,DS_MAT_A,&T));
      PetscCall(BVSetActiveColumns(eps->V,0,k));
      PetscCall(BVOrthogonalizeVec(eps->V,y,T+k*ld,&norm,NULL));
    }

    /* theta = (v,y)_B */
    PetscCall(BVSetActiveColumns(eps->V,k,k+1));
    PetscCall(BVDotVec(eps->V,y,&theta));
    if (!power->nonlinear) {
      T[k+k*ld] = theta;
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&T));
    }

    /* Eigenvalue is already stored in function evaluations.
     * Assign eigenvalue to theta to make the rest of the code consistent
     */
    if (power->update) theta = eps->eigr[eps->nconv];
    else if (power->nonlinear) theta = 1.0/(norm*sign); /* Eigenvalue: 1/|Bx| */

    if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) { /* direct & inverse iteration */

      /* approximate eigenvalue is the Rayleigh quotient */
      eps->eigr[eps->nconv] = theta;

      /**
       * If the Newton method (update, SNES) is used, we do not compute "relerr"
       * since SNES determines the convergence.
       */
      if (PetscUnlikely(power->update)) relerr = 0.;
      else {
        /* compute relative error as ||y-theta v||_2/|theta| */
        PetscCall(VecCopy(y,e));
        PetscCall(BVGetColumn(eps->V,k,&v));
        PetscCall(VecAXPY(e,power->nonlinear?-1.0:-theta,v));
        PetscCall(BVRestoreColumn(eps->V,k,&v));
        PetscCall(VecNorm(e,NORM_2,&relerr));
        if (PetscUnlikely(power->nonlinear)) relerr *= PetscAbsScalar(theta);
        else relerr /= PetscAbsScalar(theta);
      }

    } else {  /* RQI */

      /* delta = ||y||_B */
      delta = norm;

      /* compute relative error */
      if (rho == 0.0) relerr = PETSC_MAX_REAL;
      else relerr = 1.0 / (norm*PetscAbsScalar(rho));

      /* approximate eigenvalue is the shift */
      eps->eigr[eps->nconv] = rho;

      /* compute new shift */
      if (relerr<eps->tol) {
        rho = sigma;  /* if converged, restore original shift */
        PetscCall(STSetShift(eps->st,rho));
      } else {
        rho = rho + PetscConj(theta)/(delta*delta);  /* Rayleigh quotient R(v) */
        if (power->shift_type == EPS_POWER_SHIFT_WILKINSON) {
          /* beta1 is the norm of the residual associated with R(v) */
          PetscCall(BVGetColumn(eps->V,k,&v));
          PetscCall(VecAXPY(v,-PetscConj(theta)/(delta*delta),y));
          PetscCall(BVRestoreColumn(eps->V,k,&v));
          PetscCall(BVScaleColumn(eps->V,k,1.0/delta));
          PetscCall(BVNormColumn(eps->V,k,NORM_2,&norm1));
          beta1 = norm1;

          /* alpha2 = (e'*A*e)/(beta1*beta1), where e is the residual */
          PetscCall(STGetMatrix(eps->st,0,&A));
          PetscCall(BVGetColumn(eps->V,k,&v));
          PetscCall(MatMult(A,v,e));
          PetscCall(VecDot(v,e,&alpha2));
          PetscCall(BVRestoreColumn(eps->V,k,&v));
          alpha2 = alpha2 / (beta1 * beta1);

          /* choose the eigenvalue of [rho beta1; beta1 alpha2] closest to rho */
          PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
          PetscCallBLAS("LAPACKlaev2",LAPACKlaev2_(&rho,&beta1,&alpha2,&rt1,&rt2,&cs1,&sn1));
          PetscCall(PetscFPTrapPop());
          if (PetscAbsScalar(rt1-rho) < PetscAbsScalar(rt2-rho)) rho = rt1;
          else rho = rt2;
        }
        /* update operator according to new shift */
        PetscCall(KSPSetErrorIfNotConverged(ksp,PETSC_FALSE));
        PetscCall(STSetShift(eps->st,rho));
        PetscCall(KSPGetConvergedReason(ksp,&reason));
        if (reason) {
          PetscCall(PetscInfo(eps,"Factorization failed, repeat with a perturbed shift\n"));
          rho *= 1+10*PETSC_MACHINE_EPSILON;
          PetscCall(STSetShift(eps->st,rho));
          PetscCall(KSPGetConvergedReason(ksp,&reason));
          PetscCheck(!reason,PetscObjectComm((PetscObject)ksp),PETSC_ERR_CONV_FAILED,"Second factorization failed");
        }
        PetscCall(KSPSetErrorIfNotConverged(ksp,PETSC_TRUE));
      }
    }
    eps->errest[eps->nconv] = relerr;

    /* normalize */
    if (!power->nonlinear) PetscCall(Normalize(y,norm,power->idx,power->p,power->sign_normalization,NULL));
    PetscCall(BVInsertVec(eps->V,k,y));

    if (PetscUnlikely(power->update)) {
      PetscCall(SNESGetConvergedReason(power->snes,&snesreason));
      /* For Newton eigensolver, we are ready to return once SNES converged. */
      if (snesreason>0) eps->nconv = 1;
    } else if (PetscUnlikely(relerr<eps->tol)) {   /* accept eigenpair */
      eps->nconv = eps->nconv + 1;
      if (eps->nconv<eps->nev) {
        PetscCall(EPSGetStartVector(eps,eps->nconv,&breakdown));
        if (breakdown) {
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
          break;
        }
      }
    }
    /* For Newton eigensolver, monitor will be called from SNES monitor */
    if (!power->update) PetscCall(EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,PetscMin(eps->nconv+1,eps->nev)));
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx));

    /**
     * When a customized stopping test is used, and reason can be set to be converged (EPS_CONVERGED_USER).
     * In this case, we need to increase eps->nconv to "1" so users can retrieve the solution.
     */
    if (PetscUnlikely(power->nonlinear && eps->reason>0)) eps->nconv = 1;
  }

  if (power->nonlinear) PetscCall(PetscObjectCompose((PetscObject)power->snes,"eps",NULL));
  else {
    PetscCall(DSSetDimensions(eps->ds,eps->nconv,0,0));
    PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_TS_Power(EPS eps)
{
  EPS_POWER          *power = (EPS_POWER*)eps->data;
  PetscInt           k,ld;
  Vec                v,w,y,e,z;
  KSP                ksp;
  PetscReal          relerr=1.0,relerrl,delta;
  PetscScalar        theta,rho,alpha,sigma;
  PetscBool          breakdown,breakdownl;
  KSPConvergedReason reason;

  PetscFunctionBegin;
  e = eps->work[0];
  v = eps->work[1];
  w = eps->work[2];

  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) PetscCall(STGetKSP(eps->st,&ksp));
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(EPSGetStartVector(eps,0,NULL));
  PetscCall(EPSGetLeftStartVector(eps,0,NULL));
  PetscCall(BVBiorthonormalizeColumn(eps->V,eps->W,0,NULL));
  PetscCall(BVCopyVec(eps->V,0,v));
  PetscCall(BVCopyVec(eps->W,0,w));
  PetscCall(STGetShift(eps->st,&sigma));    /* original shift */
  rho = sigma;

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    k = eps->nconv;

    /* y = OP v, z = OP' w */
    PetscCall(BVGetColumn(eps->V,k,&y));
    PetscCall(STApply(eps->st,v,y));
    PetscCall(BVRestoreColumn(eps->V,k,&y));
    PetscCall(BVGetColumn(eps->W,k,&z));
    PetscCall(STApplyHermitianTranspose(eps->st,w,z));
    PetscCall(BVRestoreColumn(eps->W,k,&z));

    /* purge previously converged eigenvectors */
    PetscCall(BVBiorthogonalizeColumn(eps->V,eps->W,k));

    /* theta = (w,y)_B */
    PetscCall(BVSetActiveColumns(eps->V,k,k+1));
    PetscCall(BVDotVec(eps->V,w,&theta));
    theta = PetscConj(theta);

    if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) { /* direct & inverse iteration */

      /* approximate eigenvalue is the Rayleigh quotient */
      eps->eigr[eps->nconv] = theta;

      /* compute relative errors as ||y-theta v||_2/|theta| and ||z-conj(theta) w||_2/|theta|*/
      PetscCall(BVCopyVec(eps->V,k,e));
      PetscCall(VecAXPY(e,-theta,v));
      PetscCall(VecNorm(e,NORM_2,&relerr));
      PetscCall(BVCopyVec(eps->W,k,e));
      PetscCall(VecAXPY(e,-PetscConj(theta),w));
      PetscCall(VecNorm(e,NORM_2,&relerrl));
      relerr = PetscMax(relerr,relerrl)/PetscAbsScalar(theta);
    }

    /* normalize */
    PetscCall(BVSetActiveColumns(eps->V,k,k+1));
    PetscCall(BVGetColumn(eps->W,k,&z));
    PetscCall(BVDotVec(eps->V,z,&alpha));
    PetscCall(BVRestoreColumn(eps->W,k,&z));
    delta = PetscSqrtReal(PetscAbsScalar(alpha));
    PetscCheck(delta!=0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_CONV_FAILED,"Breakdown in two-sided Power/RQI");
    PetscCall(BVScaleColumn(eps->V,k,1.0/PetscConj(alpha/delta)));
    PetscCall(BVScaleColumn(eps->W,k,1.0/delta));
    PetscCall(BVCopyVec(eps->V,k,v));
    PetscCall(BVCopyVec(eps->W,k,w));

    if (power->shift_type == EPS_POWER_SHIFT_RAYLEIGH) { /* RQI */

      /* compute relative error */
      if (rho == 0.0) relerr = PETSC_MAX_REAL;
      else relerr = 1.0 / PetscAbsScalar(delta*rho);

      /* approximate eigenvalue is the shift */
      eps->eigr[eps->nconv] = rho;

      /* compute new shift */
      if (relerr<eps->tol) {
        rho = sigma;  /* if converged, restore original shift */
        PetscCall(STSetShift(eps->st,rho));
      } else {
        rho = rho + PetscConj(theta)/(delta*delta);  /* Rayleigh quotient R(v) */
        /* update operator according to new shift */
        PetscCall(KSPSetErrorIfNotConverged(ksp,PETSC_FALSE));
        PetscCall(STSetShift(eps->st,rho));
        PetscCall(KSPGetConvergedReason(ksp,&reason));
        if (reason) {
          PetscCall(PetscInfo(eps,"Factorization failed, repeat with a perturbed shift\n"));
          rho *= 1+10*PETSC_MACHINE_EPSILON;
          PetscCall(STSetShift(eps->st,rho));
          PetscCall(KSPGetConvergedReason(ksp,&reason));
          PetscCheck(!reason,PetscObjectComm((PetscObject)ksp),PETSC_ERR_CONV_FAILED,"Second factorization failed");
        }
        PetscCall(KSPSetErrorIfNotConverged(ksp,PETSC_TRUE));
      }
    }
    eps->errest[eps->nconv] = relerr;

    /* if relerr<tol, accept eigenpair */
    if (relerr<eps->tol) {
      eps->nconv = eps->nconv + 1;
      if (eps->nconv<eps->nev) {
        PetscCall(EPSGetStartVector(eps,eps->nconv,&breakdown));
        PetscCall(EPSGetLeftStartVector(eps,eps->nconv,&breakdownl));
        if (breakdown || breakdownl) {
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
          break;
        }
        PetscCall(BVBiorthonormalizeColumn(eps->V,eps->W,eps->nconv,NULL));
      }
    }
    PetscCall(EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,PetscMin(eps->nconv+1,eps->nev)));
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx));
  }

  PetscCall(DSSetDimensions(eps->ds,eps->nconv,0,0));
  PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSStopping_Power(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  SNESConvergedReason snesreason;

  PetscFunctionBegin;
  if (PetscUnlikely(power->update)) {
    PetscCall(SNESGetConvergedReason(power->snes,&snesreason));
    if (snesreason < 0) {
      *reason = EPS_DIVERGED_BREAKDOWN;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  PetscCall(EPSStoppingBasic(eps,its,max_it,nconv,nev,reason,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBackTransform_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear) eps->eigr[0] = 1.0/eps->eigr[0];
  else if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) PetscCall(EPSBackTransform_Default(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_Power(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_POWER         *power = (EPS_POWER*)eps->data;
  PetscBool         flg,val;
  EPSPowerShiftType shift;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Power Options");

    PetscCall(PetscOptionsEnum("-eps_power_shift_type","Shift type","EPSPowerSetShiftType",EPSPowerShiftTypes,(PetscEnum)power->shift_type,(PetscEnum*)&shift,&flg));
    if (flg) PetscCall(EPSPowerSetShiftType(eps,shift));

    PetscCall(PetscOptionsBool("-eps_power_nonlinear","Use nonlinear inverse iteration","EPSPowerSetNonlinear",power->nonlinear,&val,&flg));
    if (flg) PetscCall(EPSPowerSetNonlinear(eps,val));

    PetscCall(PetscOptionsBool("-eps_power_update","Update residual monolithically","EPSPowerSetUpdate",power->update,&val,&flg));
    if (flg) PetscCall(EPSPowerSetUpdate(eps,val));

    PetscCall(PetscOptionsBool("-eps_power_sign_normalization","Normalize Bx with sign of first nonzero entry","EPSPowerSetSignNormalization",power->sign_normalization,&power->sign_normalization,&flg));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerSetShiftType_Power(EPS eps,EPSPowerShiftType shift)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  switch (shift) {
    case EPS_POWER_SHIFT_CONSTANT:
    case EPS_POWER_SHIFT_RAYLEIGH:
    case EPS_POWER_SHIFT_WILKINSON:
      if (power->shift_type != shift) {
        power->shift_type = shift;
        eps->state = EPS_STATE_INITIAL;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid shift type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerSetShiftType - Sets the type of shifts used during the power
   iteration. This can be used to emulate the Rayleigh Quotient Iteration
   (RQI) method.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  shift - the type of shift

   Options Database Key:
.  -eps_power_shift_type - Sets the shift type (either 'constant' or
                           'rayleigh' or 'wilkinson')

   Notes:
   By default, shifts are constant (EPS_POWER_SHIFT_CONSTANT) and the iteration
   is the simple power method (or inverse iteration if a shift-and-invert
   transformation is being used).

   A variable shift can be specified (EPS_POWER_SHIFT_RAYLEIGH or
   EPS_POWER_SHIFT_WILKINSON). In this case, the iteration behaves rather like
   a cubic converging method such as RQI.

   Level: advanced

.seealso: EPSPowerGetShiftType(), STSetShift(), EPSPowerShiftType
@*/
PetscErrorCode EPSPowerSetShiftType(EPS eps,EPSPowerShiftType shift)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,shift,2);
  PetscTryMethod(eps,"EPSPowerSetShiftType_C",(EPS,EPSPowerShiftType),(eps,shift));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerGetShiftType_Power(EPS eps,EPSPowerShiftType *shift)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *shift = power->shift_type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerGetShiftType - Gets the type of shifts used during the power
   iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  shift - the type of shift

   Level: advanced

.seealso: EPSPowerSetShiftType(), EPSPowerShiftType
@*/
PetscErrorCode EPSPowerGetShiftType(EPS eps,EPSPowerShiftType *shift)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(shift,2);
  PetscUseMethod(eps,"EPSPowerGetShiftType_C",(EPS,EPSPowerShiftType*),(eps,shift));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerSetNonlinear_Power(EPS eps,PetscBool nonlinear)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear != nonlinear) {
    power->nonlinear = nonlinear;
    eps->useds = PetscNot(nonlinear);
    eps->ops->setupsort = nonlinear? NULL: EPSSetUpSort_Default;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerSetNonlinear - Sets a flag to indicate that the problem is nonlinear.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  nonlinear - whether the problem is nonlinear or not

   Options Database Key:
.  -eps_power_nonlinear - Sets the nonlinear flag

   Notes:
   If this flag is set, the solver assumes that the problem is nonlinear,
   that is, the operators that define the eigenproblem are not constant
   matrices, but depend on the eigenvector, A(x)*x=lambda*B(x)*x. This is
   different from the case of nonlinearity with respect to the eigenvalue
   (use the NEP solver class for this kind of problems).

   The way in which nonlinear operators are specified is very similar to
   the case of PETSc's SNES solver. The difference is that the callback
   functions are provided via composed functions "formFunction" and
   "formJacobian" in each of the matrix objects passed as arguments of
   EPSSetOperators(). The application context required for these functions
   can be attached via a composed PetscContainer.

   Level: advanced

.seealso: EPSPowerGetNonlinear(), EPSSetOperators()
@*/
PetscErrorCode EPSPowerSetNonlinear(EPS eps,PetscBool nonlinear)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,nonlinear,2);
  PetscTryMethod(eps,"EPSPowerSetNonlinear_C",(EPS,PetscBool),(eps,nonlinear));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerGetNonlinear_Power(EPS eps,PetscBool *nonlinear)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *nonlinear = power->nonlinear;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerGetNonlinear - Returns a flag indicating if the problem is nonlinear.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  nonlinear - the nonlinear flag

   Level: advanced

.seealso: EPSPowerSetUpdate(), EPSPowerSetNonlinear()
@*/
PetscErrorCode EPSPowerGetNonlinear(EPS eps,PetscBool *nonlinear)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(nonlinear,2);
  PetscUseMethod(eps,"EPSPowerGetNonlinear_C",(EPS,PetscBool*),(eps,nonlinear));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerSetUpdate_Power(EPS eps,PetscBool update)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  PetscCheck(power->nonlinear,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"This option does not make sense for linear problems");
  power->update = update;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerSetUpdate - Sets a flag to indicate that the residual is updated monolithically
   for nonlinear problems. This potentially has a better convergence.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  update - whether the residual is updated monolithically or not

   Options Database Key:
.  -eps_power_update - Sets the update flag

   Level: advanced

.seealso: EPSPowerGetUpdate(), EPSPowerGetNonlinear(), EPSSetOperators()
@*/
PetscErrorCode EPSPowerSetUpdate(EPS eps,PetscBool update)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,update,2);
  PetscTryMethod(eps,"EPSPowerSetUpdate_C",(EPS,PetscBool),(eps,update));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerGetUpdate_Power(EPS eps,PetscBool *update)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *update = power->update;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerGetUpdate - Returns a flag indicating if the residual is updated monolithically
   for nonlinear problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  update - the update flag

   Level: advanced

.seealso: EPSPowerSetUpdate(), EPSPowerSetNonlinear()
@*/
PetscErrorCode EPSPowerGetUpdate(EPS eps,PetscBool *update)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(update,2);
  PetscUseMethod(eps,"EPSPowerGetUpdate_C",(EPS,PetscBool*),(eps,update));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerSetSignNormalization_Power(EPS eps,PetscBool sign_normalization)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  power->sign_normalization = sign_normalization;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerSetSignNormalization - Sets a flag to indicate whether the Bx vector
   should be normalized by the sign of the first non-zero element in the vector.
   E.g., if this is true, the post-normalization value of the first non-zero element
   in the vector is guaranteed to be positive.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  sign_normalization - whether Bx should be multiplied by the sign of the first non-zero
                        element when performing normalization steps

   Options Database Key:
.  -eps_power_sign_normalization - Sets the sign normalization flag

   Level: advanced

.seealso: EPSPowerGetSignNormalization()
@*/
PetscErrorCode EPSPowerSetSignNormalization(EPS eps,PetscBool sign_normalization)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,sign_normalization,2);
  PetscTryMethod(eps,"EPSPowerSetSignNormalization_C",(EPS,PetscBool),(eps,sign_normalization));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerGetSignNormalization_Power(EPS eps,PetscBool *sign_normalization)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *sign_normalization = power->sign_normalization;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerGetSignNormalization - Returns a flag indicating whether the Bx vector
   is normalized by the sign of the first non-zero element in the vector. E.g.,
   if this is true, the post-normalization value of the first non-zero element in
   the vector is guaranteed to be positive.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  sign_normalization - the sign normalization flag

   Level: advanced

.seealso: EPSPowerSetSignNormalization()
@*/
PetscErrorCode EPSPowerGetSignNormalization(EPS eps,PetscBool *sign_normalization)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(sign_normalization,2);
  PetscUseMethod(eps,"EPSPowerGetSignNormalization_C",(EPS,PetscBool*),(eps,sign_normalization));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerSetSNES_Power(EPS eps,SNES snes)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)snes));
  PetscCall(SNESDestroy(&power->snes));
  power->snes = snes;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerSetSNES - Associate a nonlinear solver object (SNES) to the
   eigenvalue solver (to be used in nonlinear inverse iteration).

   Collective

   Input Parameters:
+  eps  - the eigenvalue solver
-  snes - the nonlinear solver object

   Level: advanced

.seealso: EPSPowerGetSNES()
@*/
PetscErrorCode EPSPowerSetSNES(EPS eps,SNES snes)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(snes,SNES_CLASSID,2);
  PetscCheckSameComm(eps,1,snes,2);
  PetscTryMethod(eps,"EPSPowerSetSNES_C",(EPS,SNES),(eps,snes));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSPowerGetSNES_Power(EPS eps,SNES *snes)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (!power->snes) {
    PetscCall(SNESCreate(PetscObjectComm((PetscObject)eps),&power->snes));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)power->snes,(PetscObject)eps,1));
    PetscCall(SNESSetOptionsPrefix(power->snes,((PetscObject)eps)->prefix));
    PetscCall(SNESAppendOptionsPrefix(power->snes,"eps_power_"));
    PetscCall(PetscObjectSetOptions((PetscObject)power->snes,((PetscObject)eps)->options));
  }
  *snes = power->snes;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSPowerGetSNES - Retrieve the nonlinear solver object (SNES) associated
   with the eigenvalue solver.

   Not Collective

   Input Parameter:
.  eps - the eigenvalue solver

   Output Parameter:
.  snes - the nonlinear solver object

   Level: advanced

.seealso: EPSPowerSetSNES()
@*/
PetscErrorCode EPSPowerGetSNES(EPS eps,SNES *snes)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(snes,2);
  PetscUseMethod(eps,"EPSPowerGetSNES_C",(EPS,SNES*),(eps,snes));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->snes) PetscCall(SNESReset(power->snes));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear) PetscCall(SNESDestroy(&power->snes));
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetShiftType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetShiftType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetNonlinear_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetNonlinear_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetUpdate_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetUpdate_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSignNormalization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSignNormalization_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSNES_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSNES_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_Power(EPS eps,PetscViewer viewer)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (power->nonlinear) {
      if (power->sign_normalization) PetscCall(PetscViewerASCIIPrintf(viewer,"  normalizing Bx by the sign of the first nonzero element\n"));
      else PetscCall(PetscViewerASCIIPrintf(viewer,"  not normalizing Bx by the sign of the first nonzero element\n"));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  using nonlinear inverse iteration\n"));
      if (power->update) PetscCall(PetscViewerASCIIPrintf(viewer,"  updating the residual monolithically\n"));
      if (!power->snes) PetscCall(EPSPowerGetSNES(eps,&power->snes));
      PetscCall(PetscViewerASCIIPushTab(viewer));
      PetscCall(SNESView(power->snes,viewer));
      PetscCall(PetscViewerASCIIPopTab(viewer));
    } else PetscCall(PetscViewerASCIIPrintf(viewer,"  %s shifts\n",EPSPowerShiftTypes[power->shift_type]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSComputeVectors_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (eps->twosided) {
    PetscCall(EPSComputeVectors_Twosided(eps));
    PetscCall(BVNormalize(eps->V,NULL));
    PetscCall(BVNormalize(eps->W,NULL));
  } else if (!power->nonlinear) PetscCall(EPSComputeVectors_Schur(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_Power(EPS eps)
{
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  if (power->nonlinear) {
    eps->categ=EPS_CATEGORY_PRECOND;
    PetscCall(STGetKSP(eps->st,&ksp));
    /* Set ST as STPRECOND so it can carry one preconditioning matrix
     * It is useful when A and B are shell matrices
     */
    PetscCall(STSetType(eps->st,STPRECOND));
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PCSetType(pc,PCNONE));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Power(EPS eps)
{
  EPS_POWER      *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->setup          = EPSSetUp_Power;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_Power;
  eps->ops->reset          = EPSReset_Power;
  eps->ops->destroy        = EPSDestroy_Power;
  eps->ops->view           = EPSView_Power;
  eps->ops->backtransform  = EPSBackTransform_Power;
  eps->ops->computevectors = EPSComputeVectors_Power;
  eps->ops->setdefaultst   = EPSSetDefaultST_Power;
  eps->stopping            = EPSStopping_Power;
  ctx->sign_normalization  = PETSC_TRUE;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetShiftType_C",EPSPowerSetShiftType_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetShiftType_C",EPSPowerGetShiftType_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetNonlinear_C",EPSPowerSetNonlinear_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetNonlinear_C",EPSPowerGetNonlinear_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetUpdate_C",EPSPowerSetUpdate_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetUpdate_C",EPSPowerGetUpdate_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSignNormalization_C",EPSPowerSetSignNormalization_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSignNormalization_C",EPSPowerGetSignNormalization_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSNES_C",EPSPowerSetSNES_Power));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSNES_C",EPSPowerGetSNES_Power));
  PetscFunctionReturn(PETSC_SUCCESS);
}
