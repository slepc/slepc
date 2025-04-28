/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Simple nonlinear eigenproblem using the NLEIGS solver.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -split <0/1>, to select the split form in the problem definition (enabled by default)\n";

/*
   Solve T(lambda)x=0 using NLEIGS solver
      with T(lambda) = -D+sqrt(lambda)*I
      where D is the Laplacian operator in 1 dimension
      and with the interpolation interval [.01,16]
*/

#include <slepcnep.h>

/*
   User-defined routines
*/
PetscErrorCode FormFunction(NEP,PetscScalar,Mat,Mat,void*);
PetscErrorCode FormJacobian(NEP,PetscScalar,Mat,void*);
PetscErrorCode ComputeSingularities(NEP,PetscInt*,PetscScalar*,void*);

int main(int argc,char **argv)
{
  NEP            nep;             /* nonlinear eigensolver context */
  Mat            F,J,A[2];
  NEPType        type;
  PetscInt       n=100,nev,Istart,Iend,i;
  PetscBool      terse,split=PETSC_TRUE;
  RG             rg;
  FN             f[2];
  PetscScalar    coeffs;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-split",&split,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSquare root eigenproblem, n=%" PetscInt_FMT "%s\n\n",n,split?" (in split form)":""));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear eigensolver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Select the NLEIGS solver and set required options for it
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,NULL));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,0,0));
#endif
  PetscCall(NEPSetTarget(nep,1.1));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define the nonlinear problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (split) {
    /*
       Create matrices for the split form
    */
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A[0]));
    PetscCall(MatSetSizes(A[0],PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(A[0]));
    PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {
      if (i>0) PetscCall(MatSetValue(A[0],i,i-1,1.0,INSERT_VALUES));
      if (i<n-1) PetscCall(MatSetValue(A[0],i,i+1,1.0,INSERT_VALUES));
      PetscCall(MatSetValue(A[0],i,i,-2.0,INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY));

    PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,1.0,&A[1]));

    /*
       Define functions for the split form
     */
    PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
    PetscCall(FNSetType(f[0],FNRATIONAL));
    coeffs = 1.0;
    PetscCall(FNRationalSetNumerator(f[0],1,&coeffs));
    PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
    PetscCall(FNSetType(f[1],FNSQRT));
    PetscCall(NEPSetSplitOperator(nep,2,A,f,SUBSET_NONZERO_PATTERN));

  } else {
    /*
       Callback form: create matrix and set Function evaluation routine
     */
    PetscCall(MatCreate(PETSC_COMM_WORLD,&F));
    PetscCall(MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(F));
    PetscCall(MatSeqAIJSetPreallocation(F,3,NULL));
    PetscCall(MatMPIAIJSetPreallocation(F,3,NULL,1,NULL));
    PetscCall(NEPSetFunction(nep,F,F,FormFunction,NULL));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&J));
    PetscCall(MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(J));
    PetscCall(MatSeqAIJSetPreallocation(J,1,NULL));
    PetscCall(MatMPIAIJSetPreallocation(J,1,NULL,1,NULL));
    PetscCall(NEPSetJacobian(nep,J,FormJacobian,NULL));
  }

  /*
     Set solver parameters at runtime
  */
  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetType(nep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type));
  PetscCall(NEPGetDimensions(nep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(NEPDestroy(&nep));
  if (split) {
    PetscCall(MatDestroy(&A[0]));
    PetscCall(MatDestroy(&A[1]));
    PetscCall(FNDestroy(&f[0]));
    PetscCall(FNDestroy(&f[1]));
  } else {
    PetscCall(MatDestroy(&F));
    PetscCall(MatDestroy(&J));
  }
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormFunction - Computes Function matrix  T(lambda)
*/
PetscErrorCode FormFunction(NEP nep,PetscScalar lambda,Mat fun,Mat B,void *ctx)
{
  PetscInt       i,n,col[3],Istart,Iend;
  PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;
  PetscScalar    value[3],t;

  PetscFunctionBeginUser;
  /*
     Compute Function entries and insert into matrix
  */
  t = PetscSqrtScalar(lambda);
  PetscCall(MatGetSize(fun,&n,NULL));
  PetscCall(MatGetOwnershipRange(fun,&Istart,&Iend));
  if (Istart==0) FirstBlock=PETSC_TRUE;
  if (Iend==n) LastBlock=PETSC_TRUE;
  value[0]=1.0; value[1]=t-2.0; value[2]=1.0;
  for (i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1;
    PetscCall(MatSetValues(fun,1,&i,3,col,value,INSERT_VALUES));
  }
  if (LastBlock) {
    i=n-1; col[0]=n-2; col[1]=n-1;
    PetscCall(MatSetValues(fun,1,&i,2,col,value,INSERT_VALUES));
  }
  if (FirstBlock) {
    i=0; col[0]=0; col[1]=1; value[0]=t-2.0; value[1]=1.0;
    PetscCall(MatSetValues(fun,1,&i,2,col,value,INSERT_VALUES));
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

/* ------------------------------------------------------------------- */
/*
   FormJacobian - Computes Jacobian matrix  T'(lambda)
*/
PetscErrorCode FormJacobian(NEP nep,PetscScalar lambda,Mat jac,void *ctx)
{
  Vec            d;

  PetscFunctionBeginUser;
  PetscCall(MatCreateVecs(jac,&d,NULL));
  PetscCall(VecSet(d,0.5/PetscSqrtScalar(lambda)));
  PetscCall(MatDiagonalSet(jac,d,INSERT_VALUES));
  PetscCall(VecDestroy(&d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
   ComputeSingularities - Computes maxnp points (at most) in the complex plane where
   the function T(.) is not analytic.

   In this case, we discretize the singularity region (-inf,0)~(-1e+6,-1e-5)
*/
PetscErrorCode ComputeSingularities(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *pt)
{
  PetscReal h;
  PetscInt  i;

  PetscFunctionBeginUser;
  h = 11.0/(*maxnp-1);
  xi[0] = -1e-5; xi[*maxnp-1] = -1e+6;
  for (i=1;i<*maxnp-1;i++) xi[i] = -PetscPowReal(10,-5+h*i);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -nep_nev 3 -terse
      output_file: output/ex27_1.out
      requires: !single
      filter: sed -e "s/[+-]0\.0*i//g"
      test:
         suffix: 1
         args: -nep_nleigs_interpolation_degree 90
      test:
         suffix: 3
         args: -nep_tol 1e-8 -nep_nleigs_rk_shifts 1.06,1.1,1.12,1.15 -nep_conv_norm -nep_nleigs_interpolation_degree 20
      test:
         suffix: 5_cuda
         args: -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 5_hip
         args: -mat_type aijhipsparse
         requires: hip

   testset:
      args: -split 0 -nep_nev 3 -terse
      output_file: output/ex27_2.out
      filter: sed -e "s/[+-]0\.0*i//g"
      test:
         suffix: 2
         args: -nep_nleigs_interpolation_degree 90
         requires: !single
      test:
         suffix: 4
         args: -nep_nleigs_rk_shifts 1.06,1.1,1.12,1.15 -nep_nleigs_interpolation_degree 20
         requires: double
      test:
         suffix: 6_cuda
         args: -mat_type aijcusparse
         requires: cuda !single
      test:
         suffix: 6_hip
         args: -mat_type aijhipsparse
         requires: hip !single

   testset:
      args: -split 0 -nep_type ciss -nep_ciss_extraction {{ritz hankel caa}} -rg_type ellipse -rg_ellipse_center 8 -rg_ellipse_radius .7 -nep_ciss_moments 4 -rg_ellipse_vscale 0.1 -terse
      requires: complex !single
      output_file: output/ex27_7.out
      timeoutfactor: 2
      test:
         suffix: 7
      test:
         suffix: 7_par
         nsize: 2
         args: -nep_ciss_partitions 2

   testset:
      args: -nep_type ciss -rg_type ellipse -rg_ellipse_center 8 -rg_ellipse_radius .7 -rg_ellipse_vscale 0.1 -terse
      requires: complex
      filter: sed -e "s/ (in split form)//" | sed -e "s/56925/56924/" | sed -e "s/60753/60754/" | sed -e "s/92630/92629/" | sed -e "s/24705/24706/"
      output_file: output/ex27_7.out
      timeoutfactor: 2
      test:
         suffix: 8
      test:
         suffix: 8_parallel
         nsize: 4
         args: -nep_ciss_partitions 4 -ds_parallel distributed
      test:
         suffix: 8_hpddm
         args: -nep_ciss_ksp_type hpddm
         requires: hpddm

   test:
      suffix: 9
      args: -nep_nev 4 -n 20 -terse
      requires: !single
      filter: sed -e "s/[+-]0\.0*i//g"

TEST*/
