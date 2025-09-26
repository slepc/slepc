/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates how to obtain invariant subspaces.\n\n"
  "Based on ex9.\n"
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

/* Matrix operations */
PetscErrorCode MatMult_Brussel(Mat,Vec,Vec);

typedef struct {
  Mat         T;
  Vec         x1,x2,y1,y2;
  PetscScalar alpha,beta,tau1,tau2,sigma;
} CTX_BRUSSEL;

int main(int argc,char **argv)
{
  EPS            eps;
  Mat            A;
  Vec            *Q,v;
  PetscScalar    delta1,delta2,L,h,kr,ki;
  PetscReal      errest,tol,re,im,lev;
  PetscInt       N=30,n,i,j,Istart,Iend,nev,nconv;
  CTX_BRUSSEL    *ctx;
  PetscBool      errok,trueres;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%" PetscInt_FMT "\n\n",N));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscNew(&ctx));
  ctx->alpha = 2.0;
  ctx->beta  = 5.45;
  delta1     = 0.008;
  delta2     = 0.004;
  L          = 0.51302;

  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-alpha",&ctx->alpha,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-beta",&ctx->beta,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL));

  /* Create matrix T */
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

  /* Fill the remaining information in the shell matrix context
     and create auxiliary vectors */
  h = 1.0 / (PetscReal)(N+1);
  ctx->tau1 = delta1 / ((h*L)*(h*L));
  ctx->tau2 = delta2 / ((h*L)*(h*L));
  ctx->sigma = 0.0;
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x1));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x2));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y1));
  PetscCall(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y2));

  /* Create the shell matrix */
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,2*n,2*n,2*N,2*N,(void*)ctx,&A));
  PetscCall(MatShellSetOperation(A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Brussel));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSSetTrueResidual(eps,PETSC_FALSE));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSSolve(eps));

  PetscCall(EPSGetTrueResidual(eps,&trueres));
  /*if (trueres) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Computing true residuals explicitly\n\n"));*/

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(EPSGetTolerances(eps,&tol,NULL));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv<nev) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem: less than %" PetscInt_FMT " eigenvalues converged\n\n",nev));
  else {
    /* Check that all converged eigenpairs satisfy the requested tolerance
       (in this example we use the solver's error estimate instead of computing
       the residual norm explicitly) */
    errok = PETSC_TRUE;
    for (i=0;i<nev;i++) {
      PetscCall(EPSGetErrorEstimate(eps,i,&errest));
      PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,NULL,NULL));
      errok = (errok && errest<5.0*SlepcAbsEigenvalue(kr,ki)*tol)? PETSC_TRUE: PETSC_FALSE;
    }
    if (!errok) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem: some of the first %" PetscInt_FMT " relative errors are higher than the tolerance\n\n",nev));
    else {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD," All requested eigenvalues computed up to the required tolerance:"));
      for (i=0;i<=(nev-1)/8;i++) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n     "));
        for (j=0;j<PetscMin(8,nev-8*i);j++) {
          PetscCall(EPSGetEigenpair(eps,8*i+j,&kr,&ki,NULL,NULL));
#if defined(PETSC_USE_COMPLEX)
          re = PetscRealPart(kr);
          im = PetscImaginaryPart(kr);
#else
          re = kr;
          im = ki;
#endif
          if (im!=0.0 && PetscAbs(re)/PetscAbs(im)<PETSC_SMALL) re = 0.0;
          if (re!=0.0 && PetscAbs(im)/PetscAbs(re)<PETSC_SMALL) im = 0.0;
          if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%.5f%+.5fi",(double)re,(double)im));
          else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%.5f",(double)re));
          if (8*i+j+1<nev) PetscCall(PetscPrintf(PETSC_COMM_WORLD,", "));
        }
      }
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n\n"));
    }
  }

  /* Get an orthogonal basis of the invariant subspace and check it is indeed
     orthogonal (note that eigenvectors are not orthogonal in this case) */
  if (nconv>1) {
    PetscCall(MatCreateVecs(A,&v,NULL));
    PetscCall(VecDuplicateVecs(v,nconv,&Q));
    PetscCall(EPSGetInvariantSubspace(eps,Q));
    PetscCall(VecCheckOrthonormality(Q,nconv,NULL,nconv,NULL,NULL,&lev));
    if (lev<10*tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality below the tolerance\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)lev));
    PetscCall(VecDestroyVecs(nconv,&Q));
    PetscCall(VecDestroy(&v));
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

/*TEST

   test:
      suffix: 1
      args: -eps_nev 4 -eps_true_residual {{0 1}}
      requires: !single

   test:
      suffix: 2
      args: -eps_nev 4 -eps_true_residual -eps_balance oneside -eps_tol 1e-7
      requires: !single

   test:
      suffix: 3
      args: -n 50 -eps_nev 4 -eps_ncv 16 -eps_type subspace -eps_largest_magnitude -bv_orthog_block {{gs tsqr chol tsqrchol svqb}}
      requires: !single

TEST*/
