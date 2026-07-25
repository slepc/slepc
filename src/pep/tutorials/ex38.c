/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum slicing on quadratic symmetric eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3]; /* problem matrices */
  PEP            pep;        /* polynomial eigenproblem solver context */
  ST             st;         /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PEPType        type;
  PetscBool      show=PETSC_FALSE,terse;
  PetscInt       n=100,Istart,Iend,nev,i,*inertias,ns;
  PetscReal      mu=1,tau=10,kappa=5,inta,intb,*shifts;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-show_inertias",&show,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum slicing on PEP, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(K,i,i-1,-kappa,INSERT_VALUES));
    PetscCall(MatSetValue(K,i,i,kappa*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(K,i,i+1,-kappa,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

  /* C is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(C,i,i-1,-tau,INSERT_VALUES));
    PetscCall(MatSetValue(C,i,i,tau*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(C,i,i+1,-tau,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  /* M is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&M));
  PetscCall(MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(M));
  PetscCall(MatGetOwnershipRange(M,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(M,i,i,mu,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));

  /*
     Set operators and set problem type
  */
  A[0] = K; A[1] = C; A[2] = M;
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPSetProblemType(pep,PEP_HYPERBOLIC));

  /*
     Set interval for spectrum slicing
  */
  inta = -11.3;
  intb = -9.5;
  PetscCall(PEPSetInterval(pep,inta,intb));
  PetscCall(PEPSetWhichEigenpairs(pep,PEP_ALL));

  /*
     Spectrum slicing requires STOAR
  */
  PetscCall(PEPSetType(pep,PEPSTOAR));

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */
  PetscCall(PEPGetST(pep,&st));
  PetscCall(STSetType(st,STSINVERT));

  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCCHOLESKY));

  /*
     Use MUMPS if available.
     Note that in complex scalars we cannot use MUMPS for spectrum slicing,
     because MatGetInertia() is not available in that case.
  */
#if PetscDefined(HAVE_MUMPS) && !PetscDefined(USE_COMPLEX)
  PetscCall(PEPSTOARSetDetectZeros(pep,PETSC_TRUE));  /* enforce zero detection */
  PetscCall(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));
  /*
     Add several MUMPS options (see ex43.c for a better way of setting them in program):
     '-st_mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
     '-st_mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
     '-st_mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

     Note: depending on the interval, it may be necessary also to increase the workspace:
     '-st_mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
  */
  PetscCall(PetscOptionsInsertString(NULL,"-st_mat_mumps_icntl_13 1 -st_mat_mumps_icntl_24 1 -st_mat_mumps_cntl_3 1e-12"));
#endif

  /*
     Set solver parameters at runtime
  */
  PetscCall(PEPSetFromOptions(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPSetUp(pep));
  if (show) {
    PetscCall(PEPSTOARGetInertias(pep,&ns,&shifts,&inertias));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n"));
    for (i=0;i<ns;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %" PetscInt_FMT " \n",(double)shifts[i],inertias[i]));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree(shifts));
    PetscCall(PetscFree(inertias));
  }
  PetscCall(PEPSolve(pep));
  if (show && !terse) {
    PetscCall(PEPSTOARGetInertias(pep,&ns,&shifts,&inertias));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"All shifts (after solve):\n"));
    for (i=0;i<ns;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %" PetscInt_FMT " \n",(double)shifts[i],inertias[i]));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree(shifts));
    PetscCall(PetscFree(inertias));
  }

  /*
     Show eigenvalues in interval and print solution
  */
  PetscCall(PEPGetType(pep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(PEPGetDimensions(pep,&nev,NULL,NULL));
  PetscCall(PEPGetInterval(pep,&inta,&intb));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," %" PetscInt_FMT " eigenvalues found in [%g, %g]\n",nev,(double)inta,(double)intb));

  /*
     Show detailed info unless -terse option is given by user
   */
  if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&K));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      requires: !single
      args: -show_inertias -terse
      output_file: output/ex38_1.out
      test:
         suffix: 1
         requires: !complex
      test:
         suffix: 1_complex
         requires: complex !mumps

   test:
      suffix: 2
      args: -pep_interval 0,2

TEST*/
