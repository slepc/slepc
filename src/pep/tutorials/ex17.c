/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a polynomial eigenproblem P(l)x = 0 with matrices loaded from a file.\n\n"
  "The command line options are:\n"
  "-A <filename1,filename2, ...> , where <filename1,.. > = matrices A0 ... files in PETSc binary form.\n\n";

#include <slepcpep.h>

#define MAX_MATRICES 40

int main(int argc,char **argv)
{
  Mat            A[MAX_MATRICES]; /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscReal      tol;
  PetscInt       nev,maxit,its,nmat=MAX_MATRICES,i;
  char*          filenames[MAX_MATRICES];
  PetscViewer    viewer;
  PetscBool      flg,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the polynomial eigenproblem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nPolynomial eigenproblem stored in file.\n\n"));
#if PetscDefined(USE_COMPLEX)
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...\n"));
#else
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...\n"));
#endif
  PetscCall(PetscOptionsGetStringArray(NULL,NULL,"-A",filenames,&nmat,&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Must indicate a comma-separated list of file names with the -A option");
  for (i=0;i<nmat;i++) {
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filenames[i],FILE_MODE_READ,&viewer));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A[i]));
    PetscCall(MatSetFromOptions(A[i]));
    PetscCall(MatLoad(A[i],viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));

  /*
     Set matrices
  */
  PetscCall(PEPSetOperators(pep,nmat,A));
  /*
     Set solver parameters at runtime
  */
  PetscCall(PEPSetFromOptions(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPSolve(pep));
  PetscCall(PEPGetIterationNumber(pep,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(PEPGetDimensions(pep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(PEPGetTolerances(pep,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PEPDestroy(&pep));
  for (i=0;i<nmat;i++) {
    PetscCall(MatDestroy(&A[i]));
    PetscCall(PetscFree(filenames[i]));
  }
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      args: -A ${wSLEPC_DIR}/share/slepc/datafiles/matrices/speaker107k.petsc,${wSLEPC_DIR}/share/slepc/datafiles/matrices/speaker107c.petsc,${wSLEPC_DIR}/share/slepc/datafiles/matrices/speaker107m.petsc -pep_type {{toar qarnoldi linear}} -pep_nev 4 -pep_ncv 20 -pep_scale scalar -terse
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)

TEST*/
