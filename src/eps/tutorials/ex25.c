/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum slicing on generalized symmetric eigenproblem.\n\n"
  "The problem is similar to ex13.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A,B;         /* matrices */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  EPSType        type;
  PetscInt       N,n=10,m,Istart,Iend,II,nev,i,j,*inertias,ns;
  PetscReal      inta,intb,*shifts;
  PetscBool      flag,show=PETSC_FALSE,terse;
#if defined(PETSC_HAVE_MUMPS) && !defined(PETSC_USE_COMPLEX)
  Mat            F;
#endif

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-show_inertias",&show,NULL));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum slicing on GHEP, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(B));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,II,II,4.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators and set problem type
  */
  PetscCall(EPSSetOperators(eps,A,B));
  PetscCall(EPSSetProblemType(eps,EPS_GHEP));

  /*
     Set interval for spectrum slicing
  */
  inta = 0.1;
  intb = 0.2;
  PetscCall(EPSSetInterval(eps,inta,intb));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_ALL));

  /*
     Spectrum slicing requires Krylov-Schur
  */
  PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */

  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSINVERT));
  PetscCall(EPSKrylovSchurGetKSP(eps,&ksp));
  PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCCHOLESKY));

  /*
     Use MUMPS if available.
     Note that in complex scalars we cannot use MUMPS for spectrum slicing,
     because MatGetInertia() is not available in that case.
  */
#if defined(PETSC_HAVE_MUMPS) && !defined(PETSC_USE_COMPLEX)
  PetscCall(EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE));  /* enforce zero detection */
  PetscCall(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));
  PetscCall(PCFactorSetUpMatSolverType(pc));
  /*
     Set several MUMPS options, the corresponding command-line options are:
     '-st_mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
     '-st_mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
     '-st_mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

     Note: depending on the interval, it may be necessary also to increase the workspace:
     '-st_mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
  */
  PetscCall(PCFactorGetMatrix(pc,&F));
  PetscCall(MatMumpsSetIcntl(F,13,1));
  PetscCall(MatMumpsSetIcntl(F,24,1));
  PetscCall(MatMumpsSetCntl(F,3,1e-12));
#endif

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSetUp(eps));
  if (show) {
    PetscCall(EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n"));
    for (i=0;i<ns;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %" PetscInt_FMT " \n",(double)shifts[i],inertias[i]));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree(shifts));
    PetscCall(PetscFree(inertias));
  }
  PetscCall(EPSSolve(eps));
  if (show) {
    PetscCall(EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"All shifts (after solve):\n"));
    for (i=0;i<ns;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %" PetscInt_FMT " \n",(double)shifts[i],inertias[i]));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree(shifts));
    PetscCall(PetscFree(inertias));
  }

  /*
     Show eigenvalues in interval and print solution
  */
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(EPSGetInterval(eps,&inta,&intb));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," %" PetscInt_FMT " eigenvalues found in [%g, %g]\n",nev,(double)inta,(double)intb));

  /*
     Show detailed info unless -terse option is given by user
   */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -terse -show_inertias
      output_file: output/ex25_1.out
      test:
         suffix: 1
         requires: !mumps
      test:
         suffix: 1_mumps
         requires: mumps !complex

TEST*/
