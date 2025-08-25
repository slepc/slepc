/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the use of a region for filtering; the number of wanted eigenvalues is not known a priori.\n\n"
  "The problem is the Brusselator wave model as in ex9.c.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>

/*
   This example tries to compute all eigenvalues lying outside the real axis.
   This could be achieved by computing LARGEST_IMAGINARY eigenvalues, but
   here we take a different route: define a region of the complex plane where
   eigenvalues must be emphasized (eigenvalues outside the region are filtered
   out). In this case, we select the region as the complement of a thin stripe
   around the real axis.
 */

PetscErrorCode MatMult_Brussel(Mat,Vec,Vec);
PetscErrorCode MyStoppingTest(EPS,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*);

typedef struct {
  Mat         T;
  Vec         x1,x2,y1,y2;
  PetscScalar alpha,beta,tau1,tau2,sigma;
  PetscInt    lastnconv;      /* last value of nconv; used in stopping test */
  PetscInt    nreps;          /* number of repetitions of nconv; used in stopping test */
} CTX_BRUSSEL;

int main(int argc,char **argv)
{
  Mat            A;               /* eigenvalue problem matrix */
  EPS            eps;             /* eigenproblem solver context */
  RG             rg;              /* region object */
  PetscScalar    delta1,delta2,L,h;
  PetscInt       N=30,n,i,Istart,Iend,mpd;
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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and configure the region
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));

  /*
     Define the region containing the eigenvalues of interest
  */
  PetscCall(EPSGetRG(eps,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
  PetscCall(RGIntervalSetEndpoints(rg,-PETSC_INFINITY,PETSC_INFINITY,-0.01,0.01));
  PetscCall(RGSetComplement(rg,PETSC_TRUE));
  /* sort eigenvalue approximations wrt a target, otherwise convergence will be erratic */
  PetscCall(EPSSetTarget(eps,0.0));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE));

  /*
     Set solver options. In particular, we must allocate sufficient
     storage for all eigenpairs that may converge (ncv). This is
     application-dependent.
  */
  mpd = 40;
  PetscCall(EPSSetDimensions(eps,2*mpd,3*mpd,mpd));
  PetscCall(EPSSetTolerances(eps,1e-7,2000));
  ctx->lastnconv = 0;
  ctx->nreps     = 0;
  PetscCall(EPSSetStoppingTestFunction(eps,MyStoppingTest,(void*)ctx,NULL));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(EPSConvergedReasonView(eps,viewer));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (!terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,viewer));
  PetscCall(PetscViewerPopFormat(viewer));

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
  PetscCall(VecAXPY(ctx->y1,ctx->beta - 1.0 + ctx->sigma,ctx->x1));
  PetscCall(VecAXPY(ctx->y1,ctx->alpha * ctx->alpha,ctx->x2));

  PetscCall(MatMult(ctx->T,ctx->x2,ctx->y2));
  PetscCall(VecScale(ctx->y2,ctx->tau2));
  PetscCall(VecAXPY(ctx->y2,-ctx->beta,ctx->x1));
  PetscCall(VecAXPY(ctx->y2,-ctx->alpha * ctx->alpha + ctx->sigma,ctx->x2));

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscCall(VecResetArray(ctx->x1));
  PetscCall(VecResetArray(ctx->x2));
  PetscCall(VecResetArray(ctx->y1));
  PetscCall(VecResetArray(ctx->y2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Function for user-defined stopping test.

    Ignores the value of nev. It only takes into account the number of
    eigenpairs that have converged in recent outer iterations (restarts);
    if no new eigenvalues have converged in the last few restarts,
    we stop the iteration, assuming that no more eigenvalues are present
    inside the region.
*/
PetscErrorCode MyStoppingTest(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ptr)
{
  CTX_BRUSSEL    *ctx = (CTX_BRUSSEL*)ptr;

  PetscFunctionBeginUser;
  /* check usual termination conditions, but ignoring the case nconv>=nev */
  PetscCall(EPSStoppingBasic(eps,its,max_it,nconv,PETSC_INT_MAX,reason,NULL));
  if (*reason==EPS_CONVERGED_ITERATING) {
    /* check if nconv is the same as before */
    if (nconv==ctx->lastnconv) ctx->nreps++;
    else {
      ctx->lastnconv = nconv;
      ctx->nreps     = 0;
    }
    /* check if no eigenvalues converged in last 10 restarts */
    if (nconv && ctx->nreps>10) *reason = EPS_CONVERGED_USER;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -n 100 -terse
      requires: !single

TEST*/
