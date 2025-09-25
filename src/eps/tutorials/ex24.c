/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum folding for a standard symmetric eigenproblem.\n\n"
  "The problem matrix is the 2-D Laplacian.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

/*
   User context for spectrum folding
*/
typedef struct {
  Mat       A;
  Vec       w;
  PetscReal target;
} CTX_FOLD;

/*
   Auxiliary routines
*/
PetscErrorCode MatMult_Fold(Mat,Vec,Vec);
PetscErrorCode RayleighQuotient(Mat,Vec,PetscScalar*);
PetscErrorCode ComputeResidualNorm(Mat,PetscScalar,Vec,PetscReal*);

int main(int argc,char **argv)
{
  Mat            A,M,P;       /* problem matrix, shell matrix and preconditioner */
  Vec            x;           /* eigenvector */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation */
  KSP            ksp;
  PC             pc;
  EPSType        type;
  CTX_FOLD       *ctx;
  PetscInt       nconv,N,n=10,m,nloc,mloc,Istart,Iend,II,i,j;
  PetscReal      *error,*evals,target=0.0,tol;
  PetscScalar    lambda;
  PetscBool      flag,terse,errok,hasmat;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-target",&target,NULL));
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum Folding, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid) target=%f\n\n",N,n,m,(double)target));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the 5-point stencil Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateVecs(A,&x,NULL));
  PetscCall(MatGetLocalSize(A,&nloc,&mloc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create shell matrix to perform spectrum folding
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PetscNew(&ctx));
  ctx->A = A;
  ctx->target = target;
  PetscCall(VecDuplicate(x,&ctx->w));

  PetscCall(MatCreateShell(PETSC_COMM_WORLD,nloc,mloc,N,N,ctx,&M));
  PetscCall(MatShellSetOperation(M,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Fold));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,M,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSetFromOptions(eps));

  PetscCall(PetscObjectTypeCompareAny((PetscObject)eps,&flag,EPSGD,EPSJD,EPSBLOPEX,EPSLOBPCG,EPSRQCG,""));
  if (flag) {
    /*
       Build preconditioner specific for this application (diagonal of A^2)
    */
    PetscCall(MatGetRowSum(A,x));
    PetscCall(VecScale(x,-1.0));
    PetscCall(VecShift(x,20.0));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&P));
    PetscCall(MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,N,N));
    PetscCall(MatSetFromOptions(P));
    PetscCall(MatDiagonalSet(P,x,INSERT_VALUES));
    /*
       Set diagonal preconditioner
    */
    PetscCall(EPSGetST(eps,&st));
    PetscCall(STSetType(st,STPRECOND));
    PetscCall(STSetPreconditionerMat(st,P));
    PetscCall(MatDestroy(&P));
    PetscCall(STGetKSP(st,&ksp));
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PCSetType(pc,PCJACOBI));
    PetscCall(STPrecondGetKSPHasMat(st,&hasmat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Preconditioned solver, hasmat=%s\n",hasmat?"true":"false"));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetTolerances(eps,&tol,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));
  if (nconv>0) {
    PetscCall(PetscMalloc2(nconv,&evals,nconv,&error));
    for (i=0;i<nconv;i++) {
      /*  Get i-th eigenvector, compute eigenvalue approximation from
          Rayleigh quotient and compute residual norm */
      PetscCall(EPSGetEigenpair(eps,i,NULL,NULL,x,NULL));
      PetscCall(RayleighQuotient(A,x,&lambda));
      PetscCall(ComputeResidualNorm(A,lambda,x,&error[i]));
#if defined(PETSC_USE_COMPLEX)
      evals[i] = PetscRealPart(lambda);
#else
      evals[i] = lambda;
#endif
    }
    PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
    if (!terse) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
           "           k              ||Ax-kx||\n"
           "   ----------------- ------------------\n"));
      for (i=0;i<nconv;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12.2g\n",(double)evals[i],(double)error[i]));
    } else {
      errok = PETSC_TRUE;
      for (i=0;i<nconv;i++) errok = (errok && error[i]<5.0*tol)? PETSC_TRUE: PETSC_FALSE;
      if (!errok) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Problem: some of the first %" PetscInt_FMT " relative errors are higher than the tolerance\n\n",nconv));
      else {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD," nconv=%" PetscInt_FMT " eigenvalues computed up to the required tolerance:",nconv));
        for (i=0;i<nconv;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %.5f",(double)evals[i]));
      }
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscFree2(evals,error));
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&M));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(VecDestroy(&x));
  PetscCall(PetscFree(ctx));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
    Matrix-vector product subroutine for the spectrum folding.
       y <-- (A-t*I)^2*x
 */
PetscErrorCode MatMult_Fold(Mat M,Vec x,Vec y)
{
  CTX_FOLD       *ctx;
  PetscScalar    sigma;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(M,&ctx));
  sigma = -ctx->target;
  PetscCall(MatMult(ctx->A,x,ctx->w));
  PetscCall(VecAXPY(ctx->w,sigma,x));
  PetscCall(MatMult(ctx->A,ctx->w,y));
  PetscCall(VecAXPY(y,sigma,ctx->w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Computes the Rayleigh quotient of a vector x
       r <-- x^T*A*x       (assumes x has unit norm)
 */
PetscErrorCode RayleighQuotient(Mat A,Vec x,PetscScalar *r)
{
  Vec            Ax;

  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(x,&Ax));
  PetscCall(MatMult(A,x,Ax));
  PetscCall(VecDot(Ax,x,r));
  PetscCall(VecDestroy(&Ax));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Computes the residual norm of an approximate eigenvector x, |A*x-lambda*x|
 */
PetscErrorCode ComputeResidualNorm(Mat A,PetscScalar lambda,Vec x,PetscReal *r)
{
  Vec            Ax;

  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(x,&Ax));
  PetscCall(MatMult(A,x,Ax));
  PetscCall(VecAXPY(Ax,-lambda,x));
  PetscCall(VecNorm(Ax,NORM_2,r));
  PetscCall(VecDestroy(&Ax));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -n 15 -eps_nev 1 -eps_ncv 12 -eps_max_it 1000 -eps_tol 1e-5 -terse
      filter: grep -v Solution
      test:
         suffix: 1
      test:
         suffix: 1_lobpcg
         args: -eps_type lobpcg
         requires: !single
      test:
         suffix: 1_gd
         args: -eps_type gd
         requires: !single

TEST*/
