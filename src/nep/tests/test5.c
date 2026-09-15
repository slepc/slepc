/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the INTERPOL solver with a user-provided PEP.\n\n"
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

int main(int argc,char **argv)
{
  NEP            nep;
  PEP            pep;
  Mat            Id,A,B;
  FN             f1,f2,f3;
  RG             rg;
  Mat            mats[3];
  FN             funs[3];
  PetscScalar    coeffs[2],b;
  PetscInt       n=128,nev,Istart,Iend,i,deg;
  PetscReal      tau=0.001,h,a=20,xi,tol;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau",&tau,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Delay Eigenproblem, n=%" PetscInt_FMT ", tau=%g\n\n",n,(double)tau));
  h = PETSC_PI/(PetscReal)(n+1);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create a standalone PEP and RG objects with appropriate settings
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PEPSetType(pep,PEPTOAR));
  PetscCall(PEPSetFromOptions(pep));

  PetscCall(RGCreate(PETSC_COMM_WORLD,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
  PetscCall(RGIntervalSetEndpoints(rg,5,20,-0.1,0.1));
  PetscCall(RGSetFromOptions(rg));

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

  /* Customize nonlinear solver; set runtime options */
  PetscCall(NEPSetType(nep,NEPINTERPOL));
  PetscCall(NEPSetRG(nep,rg));
  PetscCall(NEPInterpolSetPEP(nep,pep));
  PetscCall(NEPInterpolGetInterpolation(nep,&tol,&deg));
  PetscCall(NEPInterpolSetInterpolation(nep,tol,deg+2));  /* increase degree of interpolation */
  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetDimensions(nep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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
  PetscCall(PEPDestroy(&pep));
  PetscCall(RGDestroy(&rg));
  PetscCall(MatDestroy(&Id));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(FNDestroy(&f1));
  PetscCall(FNDestroy(&f2));
  PetscCall(FNDestroy(&f3));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -nep_nev 3 -nep_target 5 -terse
      output_file: output/test5_1.out
      filter: sed -e "s/[+-]0\.0*i//g"
      requires: !single
      test:
         suffix: 1
      test:
         suffix: 2_cuda
         args: -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 2_hip
         args: -mat_type aijhipsparse
         requires: hip
      test:
         suffix: 2_kokkos
         args: -mat_type aijkokkos
         requires: hip kokkos_kernels
      test:
         suffix: 3
         args: -nep_view_values draw
         requires: x

TEST*/
