/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the NLEIGS solver with shell matrices.\n\n"
  "This is based on ex27.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -split <0/1>, to select the split form in the problem definition (enabled by default).\n";

/*
   Solve T(lambda)x=0 using NLEIGS solver
      with T(lambda) = -D+sqrt(lambda)*I
      where D is the Laplacian operator in 1 dimension
      and with the interpolation interval [.01,16]
*/

#include <slepcnep.h>

/* User-defined routines */
PetscErrorCode FormFunction(NEP,PetscScalar,Mat,Mat,void*);
PetscErrorCode ComputeSingularities(NEP,PetscInt*,PetscScalar*,void*);
PetscErrorCode MatMult_A0(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_A0(Mat,Vec);
PetscErrorCode MatDuplicate_A0(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatMult_A1(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_A1(Mat,Vec);
PetscErrorCode MatDuplicate_A1(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatMult_F(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_F(Mat,Vec);
PetscErrorCode MatDuplicate_F(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatDestroy_F(Mat);

typedef struct {
  PetscScalar t;  /* square root of lambda */
} MatCtx;

int main(int argc,char **argv)
{
  NEP            nep;
  KSP            *ksp;
  PC             pc;
  Mat            F,A[2];
  NEPType        type;
  PetscInt       i,n=100,nev,its,nsolve;
  PetscReal      keep,tol=PETSC_SQRT_MACHINE_EPSILON/10;
  RG             rg;
  FN             f[2];
  PetscBool      terse,flg,lock,split=PETSC_TRUE;
  PetscScalar    coeffs;
  MatCtx         *ctx;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-split",&split,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSquare root eigenproblem, n=%" PetscInt_FMT "%s\n\n",n,split?" (in split form)":""));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create NEP context, configure NLEIGS with appropriate options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));
  PetscCall(NEPSetType(nep,NEPNLEIGS));
  PetscCall(NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,NULL));
  PetscCall(NEPGetRG(nep,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,-0.001,0.001));
#else
  PetscCall(RGIntervalSetEndpoints(rg,0.01,16.0,0,0));
#endif
  PetscCall(NEPSetTarget(nep,1.1));
  PetscCall(NEPNLEIGSGetKSPs(nep,&nsolve,&ksp));
  for (i=0;i<nsolve;i++) {
   PetscCall(KSPSetType(ksp[i],KSPBICG));
   PetscCall(KSPGetPC(ksp[i],&pc));
   PetscCall(PCSetType(pc,PCJACOBI));
   PetscCall(KSPSetTolerances(ksp[i],tol,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define the nonlinear problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (split) {
    /* Create matrix A0 (tridiagonal) */
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,NULL,&A[0]));
    PetscCall(MatShellSetOperation(A[0],MATOP_MULT,(void(*)(void))MatMult_A0));
    PetscCall(MatShellSetOperation(A[0],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A0));
    PetscCall(MatShellSetOperation(A[0],MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A0));
    PetscCall(MatShellSetOperation(A[0],MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A0));

    /* Create matrix A0 (identity) */
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,NULL,&A[1]));
    PetscCall(MatShellSetOperation(A[1],MATOP_MULT,(void(*)(void))MatMult_A1));
    PetscCall(MatShellSetOperation(A[1],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A1));
    PetscCall(MatShellSetOperation(A[1],MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A1));
    PetscCall(MatShellSetOperation(A[1],MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A1));

    /* Define functions for the split form */
    PetscCall(FNCreate(PETSC_COMM_WORLD,&f[0]));
    PetscCall(FNSetType(f[0],FNRATIONAL));
    coeffs = 1.0;
    PetscCall(FNRationalSetNumerator(f[0],1,&coeffs));
    PetscCall(FNCreate(PETSC_COMM_WORLD,&f[1]));
    PetscCall(FNSetType(f[1],FNSQRT));
    PetscCall(NEPSetSplitOperator(nep,2,A,f,SUBSET_NONZERO_PATTERN));
  } else {
    /* Callback form: create shell matrix for F=A0+sqrt(lambda)*A1  */
    PetscCall(PetscNew(&ctx));
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,(void*)ctx,&F));
    PetscCall(MatShellSetOperation(F,MATOP_MULT,(void(*)(void))MatMult_F));
    PetscCall(MatShellSetOperation(F,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_F));
    PetscCall(MatShellSetOperation(F,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_F));
    PetscCall(MatShellSetOperation(F,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_F));
    PetscCall(MatShellSetOperation(F,MATOP_DESTROY,(void(*)(void))MatDestroy_F));
    /* Set Function evaluation routine */
    PetscCall(NEPSetFunction(nep,F,F,FormFunction,NULL));
  }

  /* Set solver parameters at runtime */
  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetType(nep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type));
  PetscCall(NEPGetDimensions(nep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(PetscObjectTypeCompare((PetscObject)nep,NEPNLEIGS,&flg));
  if (flg) {
    PetscCall(NEPNLEIGSGetRestart(nep,&keep));
    PetscCall(NEPNLEIGSGetLocking(nep,&lock));
    PetscCall(NEPNLEIGSGetInterpolation(nep,&tol,&its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Restart factor is %3.2f",(double)keep));
    if (lock) PetscCall(PetscPrintf(PETSC_COMM_WORLD," (locking activated)"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Divided diferences with tol=%6.2g maxit=%" PetscInt_FMT "\n",(double)tol,its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(NEPConvergedReasonView(nep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(NEPDestroy(&nep));
  if (split) {
    PetscCall(MatDestroy(&A[0]));
    PetscCall(MatDestroy(&A[1]));
    PetscCall(FNDestroy(&f[0]));
    PetscCall(FNDestroy(&f[1]));
  } else PetscCall(MatDestroy(&F));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
   FormFunction - Computes Function matrix  T(lambda)
*/
PetscErrorCode FormFunction(NEP nep,PetscScalar lambda,Mat fun,Mat B,void *ctx)
{
  MatCtx         *ctxF;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(fun,&ctxF));
  ctxF->t = PetscSqrtScalar(lambda);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   ComputeSingularities - Computes maxnp points (at most) in the complex plane where
   the function T(.) is not analytic.

   In this case, we discretize the singularity region (-inf,0)~(-1e+6,-1e-6)
*/
PetscErrorCode ComputeSingularities(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *pt)
{
  PetscReal h;
  PetscInt  i;

  PetscFunctionBeginUser;
  h = 12.0/(*maxnp-1);
  xi[0] = -1e-6; xi[*maxnp-1] = -1e+6;
  for (i=1;i<*maxnp-1;i++) xi[i] = -PetscPowReal(10,-6+h*i);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* -------------------------------- A0 ----------------------------------- */

PetscErrorCode MatMult_A0(Mat A,Vec x,Vec y)
{
  PetscInt          i,n;
  PetscMPIInt       rank,size,next,prev;
  const PetscScalar *px;
  PetscScalar       *py,upper=0.0,lower=0.0;
  MPI_Comm          comm;

  PetscFunctionBeginUser;
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCallMPI(MPI_Comm_size(comm,&size));
  PetscCallMPI(MPI_Comm_rank(comm,&rank));
  next = rank==size-1? MPI_PROC_NULL: rank+1;
  prev = rank==0? MPI_PROC_NULL: rank-1;

  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(VecGetLocalSize(x,&n));

  PetscCallMPI(MPI_Sendrecv(px,1,MPIU_SCALAR,prev,0,&lower,1,MPIU_SCALAR,next,0,comm,MPI_STATUS_IGNORE));
  PetscCallMPI(MPI_Sendrecv(px+n-1,1,MPIU_SCALAR,next,0,&upper,1,MPIU_SCALAR,prev,0,comm,MPI_STATUS_IGNORE));

  py[0] = upper-2.0*px[0]+px[1];
  for (i=1;i<n-1;i++) py[i] = px[i-1]-2.0*px[i]+px[i+1];
  py[n-1] = px[n-2]-2.0*px[n-1]+lower;
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_A0(Mat A,Vec diag)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(diag,-2.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDuplicate_A0(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt       m,n,M,N;
  MPI_Comm       comm;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCall(MatCreateShell(comm,m,n,M,N,NULL,B));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_A0));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A0));
  PetscCall(MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A0));
  PetscCall(MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* -------------------------------- A1 ----------------------------------- */

PetscErrorCode MatMult_A1(Mat A,Vec x,Vec y)
{
  PetscFunctionBeginUser;
  PetscCall(VecCopy(x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_A1(Mat A,Vec diag)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(diag,1.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDuplicate_A1(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt       m,n,M,N;
  MPI_Comm       comm;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCall(MatCreateShell(comm,m,n,M,N,NULL,B));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_A1));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A1));
  PetscCall(MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A1));
  PetscCall(MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* -------------------------------- F ----------------------------------- */

PetscErrorCode MatMult_F(Mat A,Vec x,Vec y)
{
  PetscInt          i,n;
  PetscMPIInt       rank,size,next,prev;
  const PetscScalar *px;
  PetscScalar       *py,d,upper=0.0,lower=0.0;
  MatCtx            *ctx;
  MPI_Comm          comm;

  PetscFunctionBeginUser;
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCallMPI(MPI_Comm_size(comm,&size));
  PetscCallMPI(MPI_Comm_rank(comm,&rank));
  next = rank==size-1? MPI_PROC_NULL: rank+1;
  prev = rank==0? MPI_PROC_NULL: rank-1;

  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));
  PetscCall(VecGetLocalSize(x,&n));

  PetscCallMPI(MPI_Sendrecv(px,1,MPIU_SCALAR,prev,0,&lower,1,MPIU_SCALAR,next,0,comm,MPI_STATUS_IGNORE));
  PetscCallMPI(MPI_Sendrecv(px+n-1,1,MPIU_SCALAR,next,0,&upper,1,MPIU_SCALAR,prev,0,comm,MPI_STATUS_IGNORE));

  d = -2.0+ctx->t;
  py[0] = upper+d*px[0]+px[1];
  for (i=1;i<n-1;i++) py[i] = px[i-1]+d*px[i]+px[i+1];
  py[n-1] = px[n-2]+d*px[n-1]+lower;
  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_F(Mat A,Vec diag)
{
  MatCtx         *ctx;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(VecSet(diag,-2.0+ctx->t));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDuplicate_F(Mat A,MatDuplicateOption op,Mat *B)
{
  MatCtx         *actx,*bctx;
  PetscInt       m,n,M,N;
  MPI_Comm       comm;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&actx));
  PetscCall(MatGetSize(A,&M,&N));
  PetscCall(MatGetLocalSize(A,&m,&n));
  PetscCall(PetscNew(&bctx));
  bctx->t = actx->t;
  PetscCall(PetscObjectGetComm((PetscObject)A,&comm));
  PetscCall(MatCreateShell(comm,m,n,M,N,(void*)bctx,B));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_F));
  PetscCall(MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_F));
  PetscCall(MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_F));
  PetscCall(MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_F));
  PetscCall(MatShellSetOperation(*B,MATOP_DESTROY,(void(*)(void))MatDestroy_F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatDestroy_F(Mat A)
{
  MatCtx         *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&ctx));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      nsize: {{1 2}}
      args: -nep_nev 3 -nep_tol 1e-8 -terse
      filter: sed -e "s/[+-]0\.0*i//g"
      requires: !single
      test:
         suffix: 1
         args: -nep_nleigs_locking 0 -nep_nleigs_interpolation_degree 90 -nep_nleigs_interpolation_tol 1e-8 -nep_nleigs_restart 0.4
      test:
         suffix: 2
         args: -split 0

TEST*/
