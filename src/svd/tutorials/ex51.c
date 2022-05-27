/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Computes a partial GSVD of two matrices from IR Tools example.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x and y dimensions.\n\n";

#include <slepcsvd.h>

/* LookUp: returns an index i such that X(i) <= y < X(i+1), where X = linspace(0,1,N).
   Only elements start..end-1 are considered */
PetscErrorCode LookUp(PetscInt N,PetscInt start,PetscInt end,PetscReal y,PetscInt *i)
{
  PetscInt  n=end-start,j=n/2;
  PetscReal h=1.0/(N-1);

  PetscFunctionBeginUser;
  if (y<(start+j)*h) PetscCall(LookUp(N,start,start+j,y,i));
  else if (y<(start+j+1)*h) *i = start+j;
  else PetscCall(LookUp(N,start+j,end,y,i));
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat            A,B;             /* operator matrices */
  SVD            svd;             /* singular value problem solver context */
  KSP            ksp;
  PetscInt       n=32,N,i,i2,j,k,xidx,yidx,bl,Istart,Iend,col[3];
  PetscScalar    vals[] = { 1, -2, 1 },X,Y;
  PetscBool      flg,terse;
  PetscRandom    rctx;

  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg));
  N = n*n;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGSVD of inverse interpolation problem, (%" PetscInt_FMT "+%" PetscInt_FMT ")x%" PetscInt_FMT "\n\n",N,2*N,N));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                          Build the matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&rctx));
  PetscCall(PetscRandomSetInterval(rctx,0,1));
  PetscCall(PetscRandomSetFromOptions(rctx));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));

  /* make sure that the matrix is the same irrespective of the number of MPI processes */
  PetscCall(PetscRandomSetSeed(rctx,0x12345678));
  PetscCall(PetscRandomSeed(rctx));
  for (k=0;k<Istart;k++) {
    PetscCall(PetscRandomGetValue(rctx,&X));
    PetscCall(PetscRandomGetValue(rctx,&Y));
  }

  for (k=0;k<Iend-Istart;k++) {
    PetscCall(PetscRandomGetValue(rctx,&X));
    PetscCall(LookUp(n,0,n,PetscRealPart(X),&xidx));
    X = X*(n-1)-xidx;   /* scale value to a 1-spaced grid */
    PetscCall(PetscRandomGetValue(rctx,&Y));
    PetscCall(LookUp(n,0,n,PetscRealPart(Y),&yidx));
    Y = Y*(n-1)-yidx;   /* scale value to a 1-spaced grid */
    for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {
        if (i<n-1 && j<n-1 && xidx==j && yidx==i) PetscCall(MatSetValue(A,Istart+k,i+j*n,1.0-X-Y+X*Y,ADD_VALUES));
        if (i<n-1 && j>0 && xidx==j-1 && yidx==i) PetscCall(MatSetValue(A,Istart+k,i+j*n,X-X*Y,ADD_VALUES));
        if (i>0 && j<n-1 && xidx==j && yidx==i-1) PetscCall(MatSetValue(A,Istart+k,i+j*n,Y-X*Y,ADD_VALUES));
        if (i>0 && j>0 && xidx==j-1 && yidx==i-1) PetscCall(MatSetValue(A,Istart+k,i+j*n,X*Y,ADD_VALUES));
      }
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(PetscRandomDestroy(&rctx));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,2*N,N));
  PetscCall(MatSetFromOptions(B));
  PetscCall(MatSetUp(B));

  for (i=Istart;i<Iend;i++) {
    /* upper block: kron(speye(n),T1) where T1 is tridiagonal */
    i2 = i+Istart;
    if (i%n==0) PetscCall(MatSetValue(B,i2,i,1.0,INSERT_VALUES));
    else if (i%n==n-1) {
      PetscCall(MatSetValue(B,i2,i-1,-1.0,INSERT_VALUES));
      PetscCall(MatSetValue(B,i2,i,1.0,INSERT_VALUES));
    } else {
      col[0]=i-1; col[1]=i; col[2]=i+1;
      PetscCall(MatSetValues(B,1,&i2,3,col,vals,INSERT_VALUES));
    }
    /* lower block: kron(T2,speye(n)) where T2 is tridiagonal */
    i2 = i+Iend;
    bl = i/n;  /* index of block */
    j = i-bl*n; /* index within block */
    if (bl==0 || bl==n-1) PetscCall(MatSetValue(B,i2,i,1.0,INSERT_VALUES));
    else {
      col[0]=i-n; col[1]=i; col[2]=i+n;
      PetscCall(MatSetValues(B,1,&i2,3,col,vals,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A,B));
  PetscCall(SVDSetProblemType(svd,SVD_GENERALIZED));

  PetscCall(SVDSetType(svd,SVDTRLANCZOS));
  PetscCall(SVDSetDimensions(svd,6,PETSC_DEFAULT,PETSC_DEFAULT));
  PetscCall(SVDTRLanczosSetExplicitMatrix(svd,PETSC_TRUE));
  PetscCall(SVDTRLanczosSetScale(svd,-10));
  PetscCall(SVDTRLanczosGetKSP(svd,&ksp));
  PetscCall(KSPSetTolerances(ksp,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT));

  PetscCall(SVDSetFromOptions(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the problem and print solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SVDSolve(svd));

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(SVDConvergedReasonView(svd,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      args: -svd_trlanczos_gbidiag {{upper lower}} -svd_trlanczos_oneside -terse
      requires: double

TEST*/
