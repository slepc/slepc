/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem with Bethe-Salpeter structure using shell matrices.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the blocks.\n\n";

#include <slepceps.h>

/*
   This example computes eigenvalues of a matrix

        H = [  R    C
              -C^H -R^T ],

   where R is Hermitian and C is complex symmetric. In particular, R and C have the
   following Toeplitz structure:

        R = pentadiag{a,b,c,conj(b),conj(a)}
        C = tridiag{b,d,b}

   where a,b,d are complex scalars, and c is real.
*/

/*
   User-defined routines
*/
PetscErrorCode MatMult_R(Mat R,Vec x,Vec y);
PetscErrorCode MatMultTranspose_R(Mat R,Vec x,Vec y);
PetscErrorCode MatMult_C(Mat C,Vec x,Vec y);
PetscErrorCode MatMultHermitianTranspose_C(Mat C,Vec x,Vec y);

/*
   User context for shell matrices
*/
typedef struct {
  PetscScalar a,b,c,d;
} CTX_SHELL;

int main(int argc,char **argv)
{
  Mat            H,R,C;      /* problem matrices */
  EPS            eps;        /* eigenproblem solver context */
  PetscReal      lev;
  PetscInt       n=24,i,nconv;
  PetscMPIInt    size;
  PetscBool      terse,checkorthog;
  Vec            t,*x,*y;
  CTX_SHELL      *ctx;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size==1,PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only");

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nShell Bethe-Salpeter eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               Generate the shell problem matrices R and C
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscNew(&ctx));
#if defined(PETSC_USE_COMPLEX)
  ctx->a = PetscCMPLX(-0.1,0.2);
  ctx->b = PetscCMPLX(1.0,0.5);
  ctx->d = PetscCMPLX(2.0,0.2);
#else
  ctx->a = -0.1;
  ctx->b = 1.0;
  ctx->d = 2.0;
#endif
  ctx->c = 4.5;

  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,(void*)ctx,&R));
  PetscCall(MatShellSetOperation(R,MATOP_MULT,(PetscErrorCodeFn*)MatMult_R));
  PetscCall(MatShellSetOperation(R,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMultTranspose_R));
  PetscCall(MatSetOption(R,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,(void*)ctx,&C));
  PetscCall(MatShellSetOperation(C,MATOP_MULT,(PetscErrorCodeFn*)MatMult_C));
  PetscCall(MatShellSetOperation(C,MATOP_MULT_HERMITIAN_TRANSPOSE,(PetscErrorCodeFn*)MatMultHermitianTranspose_C));
  PetscCall(MatSetOption(C,MAT_SYMMETRIC,PETSC_TRUE));

  PetscCall(MatCreateBSE(R,C,&H));
  PetscCall(MatDestroy(&R));
  PetscCall(MatDestroy(&C));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,H,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_BSE));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }

  /* check bi-orthogonality */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-checkorthog",&checkorthog));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (checkorthog && nconv>0) {
    PetscCall(MatCreateVecs(H,&t,NULL));
    PetscCall(VecDuplicateVecs(t,nconv,&x));
    PetscCall(VecDuplicateVecs(t,nconv,&y));
    for (i=0;i<nconv;i++) {
      PetscCall(EPSGetEigenvector(eps,i,x[i],NULL));
      PetscCall(EPSGetLeftEigenvector(eps,i,y[i],NULL));
    }
    PetscCall(VecCheckOrthogonality(x,nconv,y,nconv,NULL,NULL,&lev));
    if (lev<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD," Level of bi-orthogonality of eigenvectors < 100*eps\n\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD," Level of bi-orthogonality of eigenvectors: %g\n\n",(double)lev));
    PetscCall(VecDestroy(&t));
    PetscCall(VecDestroyVecs(nconv,&x));
    PetscCall(VecDestroyVecs(nconv,&y));
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&H));
  PetscCall(PetscFree(ctx));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
    Matrix-vector y = R*x.

    R = pentadiag{a,b,c,conj(b),conj(a)}
 */
PetscErrorCode MatMult_R(Mat R,Vec x,Vec y)
{
  CTX_SHELL         *ctx;
  PetscInt          n,i;
  const PetscScalar *px;
  PetscScalar       *py;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(R,&ctx));
  PetscCall(MatGetSize(R,NULL,&n));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  for (i=0;i<n;i++) {
    py[i] = ctx->c*px[i];
    if (i>1) py[i] += ctx->a*px[i-2];
    if (i>0) py[i] += ctx->b*px[i-1];
    if (i<n-1) py[i] += PetscConj(ctx->b)*px[i+1];
    if (i<n-2) py[i] += PetscConj(ctx->a)*px[i+2];
  }
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Matrix-vector y = R^T*x.

    Only needed to compute the residuals.
 */
PetscErrorCode MatMultTranspose_R(Mat R,Vec x,Vec y)
{
  Vec w;

  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(x,&w));
  PetscCall(VecCopy(x,w));
  PetscCall(VecConjugate(w));
  PetscCall(MatMult_R(R,w,y));
  PetscCall(VecConjugate(y));
  PetscCall(VecDestroy(&w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Matrix-vector y = C*x.

    C = tridiag{b,d,b}
 */
PetscErrorCode MatMult_C(Mat C,Vec x,Vec y)
{
  CTX_SHELL         *ctx;
  PetscInt          n,i;
  const PetscScalar *px;
  PetscScalar       *py;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(C,&ctx));
  PetscCall(MatGetSize(C,NULL,&n));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  for (i=0;i<n;i++) {
    py[i] = ctx->d*px[i];
    if (i>0) py[i] += ctx->b*px[i-1];
    if (i<n-1) py[i] += ctx->b*px[i+1];
  }
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Matrix-vector y = C^H*x.

    Only needed to compute the residuals.
 */
PetscErrorCode MatMultHermitianTranspose_C(Mat C,Vec x,Vec y)
{
  Vec w;

  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(x,&w));
  PetscCall(VecCopy(x,w));
  PetscCall(VecConjugate(w));
  PetscCall(MatMult_C(C,w,y));
  PetscCall(VecConjugate(y));
  PetscCall(VecDestroy(&w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -eps_nev 4 -eps_ncv 16 -eps_krylovschur_bse_type {{shao gruning}} -terse -checkorthog
      filter: sed -e "s/17496/17495/g" | sed -e "s/32172/32173/g" | sed -e "s/38566/38567/g"
      test:
         suffix: 1
         requires: complex
      test:
         suffix: 1_real
         requires: !complex

TEST*/
