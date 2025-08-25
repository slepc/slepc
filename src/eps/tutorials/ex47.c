/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Shows how to recover symmetry when solving a GHEP as non-symmetric.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

/*
   User context for shell matrix
*/
typedef struct {
  KSP       ksp;
  Mat       B;
  Vec       w;
} CTX_SHELL;

/*
    Matrix-vector product function for user matrix
       y <-- A^{-1}*B*x
    The matrix A^{-1}*B*x is not symmetric, but it is self-adjoint with respect
    to the B-inner product. Here we assume A is symmetric and B is SPD.
 */
PetscErrorCode MatMult_Sinvert0(Mat S,Vec x,Vec y)
{
  CTX_SHELL      *ctx;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(S,&ctx));
  PetscCall(MatMult(ctx->B,x,ctx->w));
  PetscCall(KSPSolve(ctx->ksp,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat               A,B,S;      /* matrices */
  EPS               eps;        /* eigenproblem solver context */
  BV                bv;
  Vec               *X,v;
  PetscReal         lev=0.0,tol=1000*PETSC_MACHINE_EPSILON;
  PetscInt          N,n=45,m,Istart,Iend,II,i,j,nconv;
  PetscBool         flag;
  CTX_SHELL         *ctx;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m=n;
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized Symmetric Eigenproblem, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,m));

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
    PetscCall(MatSetValue(B,II,II,2.0/PetscLogScalar(II+2),INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatCreateVecs(B,&v,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Create a shell matrix S = A^{-1}*B
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PetscNew(&ctx));
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ctx->ksp));
  PetscCall(KSPSetOperators(ctx->ksp,A,A));
  PetscCall(KSPSetTolerances(ctx->ksp,tol,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  PetscCall(KSPSetFromOptions(ctx->ksp));
  ctx->B = B;
  PetscCall(MatCreateVecs(A,&ctx->w,NULL));
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,(void*)ctx,&S));
  PetscCall(MatShellSetOperation(S,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Sinvert0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,S,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));  /* even though S is not symmetric */
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSetUp(eps));   /* explicitly call setup */
  PetscCall(EPSGetBV(eps,&bv));
  PetscCall(BVSetMatrix(bv,B,PETSC_FALSE));  /* set inner product matrix to recover symmetry */
  PetscCall(EPSSolve(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Display solution and check B-orthogonality
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSGetTolerances(eps,&tol,NULL));
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>1) {
    PetscCall(VecDuplicateVecs(v,nconv,&X));
    for (i=0;i<nconv;i++) PetscCall(EPSGetEigenvector(eps,i,X[i],NULL));
    PetscCall(VecCheckOrthonormality(X,nconv,NULL,nconv,B,NULL,&lev));
    if (lev<10*tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality below the tolerance\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)lev));
    PetscCall(VecDestroyVecs(nconv,&X));
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&v));
  PetscCall(KSPDestroy(&ctx->ksp));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(PetscFree(ctx));
  PetscCall(MatDestroy(&S));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      args: -n 18 -eps_nev 4 -eps_max_it 1500
      requires: !single

TEST*/
