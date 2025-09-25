/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a problem associated to the Brusselator wave model in chemical reactions, illustrating the use of shell matrices.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>

/*
   This example computes the eigenvalues with largest real part of the
   following matrix

        A = [ tau1*T+(beta-1)*I     alpha^2*I
                  -beta*I        tau2*T-alpha^2*I ],

   where

        T = tridiag{1,-2,1}
        h = 1/(n+1)
        tau1 = delta1/(h*L)^2
        tau2 = delta2/(h*L)^2
*/

/*
   Matrix operations
*/
PetscErrorCode MatMult_Brussel(Mat,Vec,Vec);
PetscErrorCode MatMultTranspose_Brussel(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_Brussel(Mat,Vec);

typedef struct {
  Mat         T;
  Vec         x1,x2,y1,y2;
  PetscScalar alpha,beta,tau1,tau2,sigma;
} CTX_BRUSSEL;

int main(int argc,char **argv)
{
  Mat            A;               /* eigenvalue problem matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscScalar    delta1,delta2,L,h;
  PetscInt       N=30,n,i,Istart,Iend,nev;
  CTX_BRUSSEL    *ctx;
  PetscBool      terse;
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%" PetscInt_FMT "\n\n",N));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create shell matrix context and set default parameters
  */
  PetscCall(PetscNew(&ctx));
  ctx->alpha = 2.0;
  ctx->beta  = 5.45;
  delta1     = 0.008;
  delta2     = 0.004;
  L          = 0.51302;

  /*
     Look the command line for user-provided parameters
  */
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-alpha",&ctx->alpha,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-beta",&ctx->beta,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL));

  /*
     Create matrix T
  */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&ctx->T));
  PetscCall(MatSetSizes(ctx->T,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(ctx->T));

  PetscCall(MatGetOwnershipRange(ctx->T,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(ctx->T,i,i-1,1.0,INSERT_VALUES));
    if (i<N-1) PetscCall(MatSetValue(ctx->T,i,i+1,1.0,INSERT_VALUES));
    PetscCall(MatSetValue(ctx->T,i,i,-2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(ctx->T,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(ctx->T,MAT_FINAL_ASSEMBLY));
  PetscCall(MatGetLocalSize(ctx->T,&n,NULL));

  /*
     Fill the remaining information in the shell matrix context
     and create auxiliary vectors
  */
  h = 1.0 / (PetscReal)(N+1);
  ctx->tau1 = delta1 / ((h*L)*(h*L));
  ctx->tau2 = delta2 / ((h*L)*(h*L));
  ctx->sigma = 0.0;
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x1));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x2));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y1));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y2));

  /*
     Create the shell matrix
  */
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,2*n,2*n,2*N,2*N,(void*)ctx,&A));
  PetscCall(MatShellSetOperation(A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Brussel));
  PetscCall(MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_Brussel));
  PetscCall(MatShellSetOperation(A,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Brussel));

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
     Ask for the rightmost eigenvalues
  */
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));

  /*
     Set other solver options at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
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
    PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
    PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,viewer));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&ctx->T));
  PetscCall(VecDestroy(&ctx->x1));
  PetscCall(VecDestroy(&ctx->x2));
  PetscCall(VecDestroy(&ctx->y1));
  PetscCall(VecDestroy(&ctx->y2));
  PetscCall(PetscFree(ctx));
  PetscCall(SlepcFinalize());
  return 0;
}

PetscErrorCode MatMult_Brussel(Mat A,Vec x,Vec y)
{
  PetscInt          n;
  const PetscScalar *px;
  PetscScalar       *py;
  CTX_BRUSSEL       *ctx;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatGetLocalSize(ctx->T,&n,NULL));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(VecPlaceArray(ctx->x1,px));
  PetscCall(VecPlaceArray(ctx->x2,px+n));
  PetscCall(VecPlaceArray(ctx->y1,py));
  PetscCall(VecPlaceArray(ctx->y2,py+n));

  PetscCall(MatMult(ctx->T,ctx->x1,ctx->y1));
  PetscCall(VecScale(ctx->y1,ctx->tau1));
  PetscCall(VecAXPY(ctx->y1,ctx->beta-1.0+ctx->sigma,ctx->x1));
  PetscCall(VecAXPY(ctx->y1,ctx->alpha*ctx->alpha,ctx->x2));

  PetscCall(MatMult(ctx->T,ctx->x2,ctx->y2));
  PetscCall(VecScale(ctx->y2,ctx->tau2));
  PetscCall(VecAXPY(ctx->y2,-ctx->beta,ctx->x1));
  PetscCall(VecAXPY(ctx->y2,-ctx->alpha*ctx->alpha+ctx->sigma,ctx->x2));

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscCall(VecResetArray(ctx->x1));
  PetscCall(VecResetArray(ctx->x2));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecResetArray(ctx->y2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMultTranspose_Brussel(Mat A,Vec x,Vec y)
{
  PetscInt          n;
  const PetscScalar *px;
  PetscScalar       *py;
  CTX_BRUSSEL       *ctx;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(MatGetLocalSize(ctx->T,&n,NULL));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(VecPlaceArray(ctx->x1,px));
  PetscCall(VecPlaceArray(ctx->x2,px+n));
  PetscCall(VecPlaceArray(ctx->y1,py));
  PetscCall(VecPlaceArray(ctx->y2,py+n));

  PetscCall(MatMultTranspose(ctx->T,ctx->x1,ctx->y1));
  PetscCall(VecScale(ctx->y1,ctx->tau1));
  PetscCall(VecAXPY(ctx->y1,ctx->beta-1.0+ctx->sigma,ctx->x1));
  PetscCall(VecAXPY(ctx->y1,-ctx->beta,ctx->x2));

  PetscCall(MatMultTranspose(ctx->T,ctx->x2,ctx->y2));
  PetscCall(VecScale(ctx->y2,ctx->tau2));
  PetscCall(VecAXPY(ctx->y2,ctx->alpha*ctx->alpha,ctx->x1));
  PetscCall(VecAXPY(ctx->y2,-ctx->alpha*ctx->alpha+ctx->sigma,ctx->x2));

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscCall(VecResetArray(ctx->x1));
  PetscCall(VecResetArray(ctx->x2));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecResetArray(ctx->y2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_Brussel(Mat A,Vec diag)
{
  Vec            d1,d2;
  PetscInt       n;
  PetscScalar    *pd;
  MPI_Comm       comm;
  CTX_BRUSSEL    *ctx;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCall(MatGetLocalSize(ctx->T,&n,NULL));
  PetscCall(VecGetArray(diag,&pd));
  PetscCall(VecCreateMPIWithArray(comm,1,n,PETSC_DECIDE,pd,&d1));
  PetscCall(VecCreateMPIWithArray(comm,1,n,PETSC_DECIDE,pd+n,&d2));

  PetscCall(VecSet(d1,-2.0*ctx->tau1 + ctx->beta - 1.0 + ctx->sigma));
  PetscCall(VecSet(d2,-2.0*ctx->tau2 - ctx->alpha*ctx->alpha + ctx->sigma));

  PetscCall(VecDestroy(&d1));
  PetscCall(VecDestroy(&d2));
  PetscCall(VecRestoreArray(diag,&pd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -n 50 -eps_nev 4 -eps_two_sided {{0 1}} -eps_type {{krylovschur lapack}} -terse
      requires: !single
      filter: grep -v method

   test:
      suffix: 2
      args: -eps_nev 8 -eps_max_it 300 -eps_target -28 -rg_type interval -rg_interval_endpoints -40,-20,-.1,.1 -terse
      requires: !single

   test:
      suffix: 3
      args: -n 50 -eps_nev 4 -eps_balance twoside -terse
      requires: double
      filter: grep -v method
      output_file: output/ex9_1.out

   test:
      suffix: 4
      args: -eps_smallest_imaginary -eps_ncv 24 -terse
      requires: !complex !single

   test:
      suffix: 4_complex
      args: -eps_smallest_imaginary -eps_ncv 24 -terse
      requires: complex !single

   test:
      suffix: 5
      args: -eps_nev 4 -eps_target_real -eps_target -3 -terse
      requires: !single

   test:
      suffix: 6
      args: -eps_nev 2 -eps_target_imaginary -eps_target 3i -terse
      requires: complex !single

   test:
      suffix: 7
      args: -n 40 -eps_nev 1 -eps_type arnoldi -eps_smallest_real -eps_refined -eps_ncv 42 -eps_max_it 300 -terse
      requires: double

   test:
      suffix: 8
      args: -eps_nev 2 -eps_target -30 -eps_type jd -st_matmode shell -eps_jd_fix 0.0001 -eps_jd_const_correction_tol 0 -terse
      requires: !single
      filter: sed -e "s/[+-]0\.0*i//g"
      timeoutfactor: 2

   test:
      suffix: 9
      args: -eps_largest_imaginary -eps_ncv 24 -terse
      requires: !single

TEST*/
