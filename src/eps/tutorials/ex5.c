/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem associated with a Markov model of a random walk on a triangular grid. "
  "It is a standard nonsymmetric eigenproblem with real eigenvalues and the rightmost eigenvalue is known to be 1.\n"
  "This example illustrates how the user can set the initial vector.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);

int main(int argc,char **argv)
{
  Vec            v0;              /* initial vector */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  EPSStop        stop;
  PetscReal      thres;
  PetscInt       N,m=15,nev;
  PetscMPIInt    rank;
  PetscBool      terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  N = m*(m+1)/2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n\n",N,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMarkovModel(m,A));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /*
     Set the initial vector. This is optional, if not done the initial
     vector is set to random values
  */
  PetscCall(MatCreateVecs(A,&v0,NULL));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  if (!rank) {
    PetscCall(VecSetValue(v0,0,1.0,INSERT_VALUES));
    PetscCall(VecSetValue(v0,1,1.0,INSERT_VALUES));
    PetscCall(VecSetValue(v0,2,1.0,INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(v0));
  PetscCall(VecAssemblyEnd(v0));
  PetscCall(EPSSetInitialSpace(eps,1,&v0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetStoppingTest(eps,&stop));
  if (stop!=EPS_STOP_THRESHOLD) {
    PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  } else {
    PetscCall(EPSGetThreshold(eps,&thres,NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Using threshold: %.4g\n",(double)thres));
  }

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
  PetscCall(VecDestroy(&v0));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
    Matrix generator for a Markov model of a random walk on a triangular grid.

    This subroutine generates a test matrix that models a random walk on a
    triangular grid. This test example was used by G. W. Stewart ["{SRRIT} - a
    FORTRAN subroutine to calculate the dominant invariant subspaces of a real
    matrix", Tech. report. TR-514, University of Maryland (1978).] and in a few
    papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34, pp. 269-295
    (1980) ]. These matrices provide reasonably easy test problems for eigenvalue
    algorithms. The transpose of the matrix  is stochastic and so it is known
    that one is an exact eigenvalue. One seeks the eigenvector of the transpose
    associated with the eigenvalue unity. The problem is to calculate the steady
    state probability distribution of the system, which is the eigevector
    associated with the eigenvalue one and scaled in such a way that the sum all
    the components is equal to one.

    Note: the code will actually compute the transpose of the stochastic matrix
    that contains the transition probabilities.
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A)
{
  const PetscReal cst = 0.5/(PetscReal)(m-1);
  PetscReal       pd,pu;
  PetscInt        Istart,Iend,i,j,jmax,ix=0;

  PetscFunctionBeginUser;
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=1;i<=m;i++) {
    jmax = m-i+1;
    for (j=1;j<=jmax;j++) {
      ix = ix + 1;
      if (ix-1<Istart || ix>Iend) continue;  /* compute only owned rows */
      if (j!=jmax) {
        pd = cst*(PetscReal)(i+j-1);
        /* north */
        if (i==1) PetscCall(MatSetValue(A,ix-1,ix,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix,pd,INSERT_VALUES));
        /* east */
        if (j==1) PetscCall(MatSetValue(A,ix-1,ix+jmax-1,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix+jmax-1,pd,INSERT_VALUES));
      }
      /* south */
      pu = 0.5 - cst*(PetscReal)(i+j-3);
      if (j>1) PetscCall(MatSetValue(A,ix-1,ix-2,pu,INSERT_VALUES));
      /* west */
      if (i>1) PetscCall(MatSetValue(A,ix-1,ix-jmax-2,pu,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      nsize: 2
      args: -eps_largest_real -eps_nev 4 -eps_two_sided {{0 1}} -eps_krylovschur_locking {{0 1}} -ds_parallel synchronized -terse
      filter: sed -e "s/90424/90423/" | sed -e "s/85715/85714/"

   test:
      suffix: 2
      args: -eps_threshold_relative .9 -eps_ncv 10 -terse
      filter: sed -e "s/-0.85714/0.85714/" | sed -e "s/90424/90423/" | sed -e "s/-1.00000, 1.00000/1.00000, -1.00000/" | sed -e "s/-0.97137, 0.97137/0.97137, -0.97137/" | sed -e "s/-0.90423, 0.90423/0.90423, -0.90423/"

TEST*/
