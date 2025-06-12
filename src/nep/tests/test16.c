/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates use of NEPSetEigenvalueComparison().\n\n"
  "This is a simplified version of ex20.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions.\n";

/*
   Solve 1-D PDE
            -u'' = lambda*u
   on [0,1] subject to
            u(0)=0, u'(1)=u(1)*lambda*kappa/(kappa-lambda)
*/

#include <slepcnep.h>

/*
   User-defined routines
*/
PetscErrorCode FormFunction(NEP,PetscScalar,Mat,Mat,void*);
PetscErrorCode MyEigenSort(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);

/*
   User-defined application context
*/
typedef struct {
  PetscScalar kappa;   /* ratio between stiffness of spring and attached mass */
  PetscReal   h;       /* mesh spacing */
} ApplicationCtx;

int main(int argc,char **argv)
{
  NEP            nep;             /* nonlinear eigensolver context */
  Mat            F;               /* Function matrix */
  ApplicationCtx ctx;             /* user-defined context */
  PetscScalar    target;
  RG             rg;
  PetscInt       n=128;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Nonlinear Eigenproblem, n=%" PetscInt_FMT "\n\n",n));
  ctx.h = 1.0/(PetscReal)n;
  ctx.kappa = 1.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               Prepare nonlinear eigensolver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&F));
  PetscCall(MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(F));
  PetscCall(MatSeqAIJSetPreallocation(F,3,NULL));
  PetscCall(MatMPIAIJSetPreallocation(F,3,NULL,1,NULL));
  PetscCall(NEPSetFunction(nep,F,F,FormFunction,&ctx));

  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,2.0,400.0,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,2.0,400.0,0,0));
#endif
  PetscCall(NEPSetTarget(nep,25.0));
  PetscCall(NEPSetEigenvalueComparison(nep,MyEigenSort,&target));
  PetscCall(NEPSetTolerances(nep,PETSC_SMALL,PETSC_CURRENT));
  PetscCall(NEPSetFromOptions(nep));
  PetscCall(NEPGetTarget(nep,&target));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Solve the eigensystem and display the solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSolve(nep));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&F));
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormFunction - Computes Function matrix  T(lambda)

   Input Parameters:
.  nep    - the NEP context
.  lambda - the scalar argument
.  ctx    - optional user-defined context, as set by NEPSetFunction()

   Output Parameters:
.  fun - Function matrix
.  B   - optionally different preconditioning matrix
*/
PetscErrorCode FormFunction(NEP nep,PetscScalar lambda,Mat fun,Mat B,void *ctx)
{
  ApplicationCtx *user = (ApplicationCtx*)ctx;
  PetscScalar    A[3],c,d;
  PetscReal      h;
  PetscInt       i,n,j[3],Istart,Iend;
  PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;

  PetscFunctionBeginUser;
  /*
     Compute Function entries and insert into matrix
  */
  PetscCall(MatGetSize(fun,&n,NULL));
  PetscCall(MatGetOwnershipRange(fun,&Istart,&Iend));
  if (Istart==0) FirstBlock=PETSC_TRUE;
  if (Iend==n) LastBlock=PETSC_TRUE;
  h = user->h;
  c = user->kappa/(lambda-user->kappa);
  d = n;

  /*
     Interior grid points
  */
  for (i=(FirstBlock? Istart+1: Istart);i<(LastBlock? Iend-1: Iend);i++) {
    j[0] = i-1; j[1] = i; j[2] = i+1;
    A[0] = A[2] = -d-lambda*h/6.0; A[1] = 2.0*(d-lambda*h/3.0);
    PetscCall(MatSetValues(fun,1,&i,3,j,A,INSERT_VALUES));
  }

  /*
     Boundary points
  */
  if (FirstBlock) {
    i = 0;
    j[0] = 0; j[1] = 1;
    A[0] = 2.0*(d-lambda*h/3.0); A[1] = -d-lambda*h/6.0;
    PetscCall(MatSetValues(fun,1,&i,2,j,A,INSERT_VALUES));
  }

  if (LastBlock) {
    i = n-1;
    j[0] = n-2; j[1] = n-1;
    A[0] = -d-lambda*h/6.0; A[1] = d-lambda*h/3.0+lambda*c;
    PetscCall(MatSetValues(fun,1,&i,2,j,A,INSERT_VALUES));
  }

  /*
     Assemble matrix
  */
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  if (fun != B) {
    PetscCall(MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, eigenvalues are sorted with respect to the target,
    but those on the right of the target are preferred.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscReal   a,b;
  PetscScalar target = *(PetscScalar*)ctx;

  PetscFunctionBeginUser;
  if (PetscRealPart(ar-target)<0.0 && PetscRealPart(br-target)>0.0) *r = 1;
  else {
    a = SlepcAbsEigenvalue(ar-target,ai);
    b = SlepcAbsEigenvalue(br-target,bi);
    if (a>b) *r = 1;
    else if (a<b) *r = -1;
    else *r = 0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -nep_nev 4 -nep_ncv 8 -terse
      requires: double !complex

TEST*/
