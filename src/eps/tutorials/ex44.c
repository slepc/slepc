/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Compute rightmost eigenvalues with Lyapunov inverse iteration.\n\n"
  "Loads matrix from a file or builds the same problem as ex36.c (with fixed parameters).\n\n"
  "The command line options are:\n"
  "  -file <filename>, where <filename> = matrix file in PETSc binary form.\n"
  "  -shift <sigma>, shift to make the matrix stable.\n"
  "  -n <n>, block dimension of the 2x2 block matrix (if matrix is generated).\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscScalar    alpha,beta,tau1,tau2,delta1,delta2,L,h,sigma=0.0;
  PetscInt       n=30,i,Istart,Iend,nev;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      flg,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetString(NULL,NULL,"-file",filename,sizeof(filename),&flg));
  if (flg) {

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Load the matrix from file
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nEigenproblem stored in file.\n\n"));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrix from a binary file...\n"));
#else
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrix from a binary file...\n"));
#endif
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatLoad(A,viewer));
    PetscCall(PetscViewerDestroy(&viewer));

  } else {

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Generate Brusselator matrix
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%" PetscInt_FMT "\n\n",n));

    alpha  = 2.0;
    beta   = 5.45;
    delta1 = 0.008;
    delta2 = 0.004;
    L      = 0.51302;

    h = 1.0 / (PetscReal)(n+1);
    tau1 = delta1 / ((h*L)*(h*L));
    tau2 = delta2 / ((h*L)*(h*L));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,2*n,2*n));
    PetscCall(MatSetFromOptions(A));

    PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {
      if (i<n) {  /* upper blocks */
        if (i>0) PetscCall(MatSetValue(A,i,i-1,tau1,INSERT_VALUES));
        if (i<n-1) PetscCall(MatSetValue(A,i,i+1,tau1,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i,-2.0*tau1+beta-1.0,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i+n,alpha*alpha,INSERT_VALUES));
      } else {  /* lower blocks */
        if (i>n) PetscCall(MatSetValue(A,i,i-1,tau2,INSERT_VALUES));
        if (i<2*n-1) PetscCall(MatSetValue(A,i,i+1,tau2,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i,-2.0*tau2-alpha*alpha,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i-n,-beta,INSERT_VALUES));
      }
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  }

  /* Shift the matrix to make it stable, A-sigma*I */
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-shift",&sigma,NULL));
  PetscCall(MatShift(A,-sigma));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 6 -shift 0.1 -eps_type {{krylovschur lyapii}} -eps_tol 1e-7 -terse
      requires: double
      filter: grep -v method | sed -e "s/-0.09981-2.13938i, -0.09981+2.13938i/-0.09981+2.13938i, -0.09981-2.13938i/" | sed -e "s/-0.77192-2.52712i, -0.77192+2.52712i/-0.77192+2.52712i, -0.77192-2.52712i/" | sed -e "s/-1.88445-3.02666i, -1.88445+3.02666i/-1.88445+3.02666i, -1.88445-3.02666i/"
      output_file: output/ex44_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -eps_lyapii_ranks 8,20 -options_left no

   test:
      suffix: 3
      args: -file ${SLEPC_DIR}/share/slepc/datafiles/matrices/rdb200.petsc -shift 35 -eps_type {{krylovschur lyapii}} -terse
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      filter: grep -v method

TEST*/
