/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the use of a user-defined stopping test.\n\n"
  "This is based on ex22.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions.\n"
  "  -tau <tau>, where <tau> is the delay parameter.\n\n";

/*
   Solve parabolic partial differential equation with time delay tau

            u_t = u_xx + a*u(t) + b*u(t-tau)
            u(0,t) = u(pi,t) = 0

   with a = 20 and b(x) = -4.1+x*(1-exp(x-pi)).

   Discretization leads to a DDE of dimension n

            -u' = A*u(t) + B*u(t-tau)

   which results in the nonlinear eigenproblem

            (-lambda*I + A + exp(-tau*lambda)*B)*u = 0
*/

#include <slepcnep.h>

/*
   User-defined routines
*/
PetscErrorCode MyStoppingTest(NEP,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*);

typedef struct {
  PetscInt    lastnconv;      /* last value of nconv; used in stopping test */
  PetscInt    nreps;          /* number of repetitions of nconv; used in stopping test */
} CTX_DELAY;

int main(int argc,char **argv)
{
  NEP            nep;
  Mat            Id,A,B;
  FN             f1,f2,f3;
  RG             rg;
  CTX_DELAY      *ctx;
  Mat            mats[3];
  FN             funs[3];
  PetscScalar    coeffs[2],b;
  PetscInt       n=128,Istart,Iend,i,mpd;
  PetscReal      tau=0.001,h,a=20,xi;
  PetscBool      terse;
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau",&tau,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Delay Eigenproblem, n=%" PetscInt_FMT ", tau=%g\n\n",n,(double)tau));
  h = PETSC_PI/(PetscReal)(n+1);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear eigensolver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));

  /* Identity matrix */
  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,1.0,&Id));
  PetscCall(MatSetOption(Id,MAT_HERMITIAN,PETSC_TRUE));

  /* A = 1/h^2*tridiag(1,-2,1) + a*I */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,1.0/(h*h),INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,1.0/(h*h),INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,-2.0/(h*h)+a,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));

  /* B = diag(b(xi)) */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));
  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    xi = (i+1)*h;
    b = -4.1+xi*(1.0-PetscExpReal(xi-PETSC_PI));
    PetscCall(MatSetValues(B,1,&i,1,&i,&b,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(B,MAT_HERMITIAN,PETSC_TRUE));

  /* Functions: f1=-lambda, f2=1.0, f3=exp(-tau*lambda) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f1));
  PetscCall(FNSetType(f1,FNRATIONAL));
  coeffs[0] = -1.0; coeffs[1] = 0.0;
  PetscCall(FNRationalSetNumerator(f1,2,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f2));
  PetscCall(FNSetType(f2,FNRATIONAL));
  coeffs[0] = 1.0;
  PetscCall(FNRationalSetNumerator(f2,1,coeffs));

  PetscCall(FNCreate(PETSC_COMM_WORLD,&f3));
  PetscCall(FNSetType(f3,FNEXP));
  PetscCall(FNSetScale(f3,-tau,1.0));

  /* Set the split operator */
  mats[0] = A;  funs[0] = f2;
  mats[1] = Id; funs[1] = f1;
  mats[2] = B;  funs[2] = f3;
  PetscCall(NEPSetSplitOperator(nep,3,mats,funs,SUBSET_NONZERO_PATTERN));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if PetscDefined(USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,5,20,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,5,20,-0.0,0.0));
#endif
  PetscCall(NEPSetTarget(nep,15.0));
  PetscCall(NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE));

  /*
     Set solver options. In particular, we must allocate sufficient
     storage for all eigenpairs that may converge (ncv). This is
     application-dependent.
  */
  mpd = 40;
  PetscCall(NEPSetDimensions(nep,2*mpd,3*mpd,mpd));
  PetscCall(NEPSetTolerances(nep,PETSC_CURRENT,2000));
  PetscCall(PetscNew(&ctx));
  ctx->lastnconv = 0;
  ctx->nreps     = 0;
  PetscCall(NEPSetStoppingTestFunction(nep,MyStoppingTest,(void*)ctx,NULL));

  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSolve(nep));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(NEPConvergedReasonView(nep,viewer));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (!terse) PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,viewer));
  PetscCall(PetscViewerPopFormat(viewer));

  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&Id));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(FNDestroy(&f1));
  PetscCall(FNDestroy(&f2));
  PetscCall(FNDestroy(&f3));
  PetscCall(PetscFree(ctx));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
    Function for user-defined stopping test.

    Ignores the value of nev. It only takes into account the number of
    eigenpairs that have converged in recent outer iterations (restarts);
    if no new eigenvalues have converged in the last few restarts,
    we stop the iteration, assuming that no more eigenvalues are present
    inside the region.
*/
PetscErrorCode MyStoppingTest(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ptr)
{
  CTX_DELAY      *ctx = (CTX_DELAY*)ptr;

  PetscFunctionBeginUser;
  /* check usual termination conditions, but ignoring the case nconv>=nev */
  PetscCall(NEPStoppingBasic(nep,its,max_it,nconv,PETSC_INT_MAX,reason,NULL));
  if (*reason==NEP_CONVERGED_ITERATING) {
    /* check if nconv is the same as before */
    if (nconv==ctx->lastnconv) ctx->nreps++;
    else {
      ctx->lastnconv = nconv;
      ctx->nreps     = 0;
    }
    /* check if no eigenvalues converged in last 10 restarts */
    if (nconv && ctx->nreps>10) *reason = NEP_CONVERGED_USER;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -terse

TEST*/
