/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "A quadratic eigenproblem defined using shell matrices.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x and y dimensions.\n\n";

#include <slepcpep.h>

/*
   User-defined routines
*/
PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);
PetscErrorCode MatMult_Zero(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Zero(Mat A,Vec diag);
PetscErrorCode MatMult_Identity(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Identity(Mat A,Vec diag);

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PEPType        type;
  PetscInt       N,n=10,nev;
  PetscMPIInt    size;
  PetscBool      terse;
  ST             st;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size==1,PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only");

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  N = n*n;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nQuadratic Eigenproblem with shell matrices, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is the 2-D Laplacian */
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,&n,&K));
  PetscCall(MatShellSetOperation(K,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Laplacian2D));
  PetscCall(MatShellSetOperation(K,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMult_Laplacian2D));
  PetscCall(MatShellSetOperation(K,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Laplacian2D));

  /* C is the zero matrix */
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,NULL,&C));
  PetscCall(MatShellSetOperation(C,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Zero));
  PetscCall(MatShellSetOperation(C,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMult_Zero));
  PetscCall(MatShellSetOperation(C,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Zero));

  /* M is the identity matrix */
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,NULL,&M));
  PetscCall(MatShellSetOperation(M,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Identity));
  PetscCall(MatShellSetOperation(M,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMult_Identity));
  PetscCall(MatShellSetOperation(M,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Identity));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));

  /*
     Set matrices and problem type
  */
  A[0] = K; A[1] = C; A[2] = M;
  PetscCall(PEPSetOperators(pep,3,A));
  PetscCall(PEPGetST(pep,&st));
  PetscCall(STSetMatMode(st,ST_MATMODE_SHELL));

  /*
     Set solver parameters at runtime
  */
  PetscCall(PEPSetFromOptions(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPSolve(pep));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(PEPGetType(pep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(PEPGetDimensions(pep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PEPDestroy(&pep));
  PetscCall(MatDestroy(&M));
  PetscCall(MatDestroy(&C));
  PetscCall(MatDestroy(&K));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
    Compute the matrix vector multiplication y<---T*x where T is a nx by nx
    tridiagonal matrix with DD on the diagonal, DL on the subdiagonal, and
    DU on the superdiagonal.
 */
static void tv(int nx,const PetscScalar *x,PetscScalar *y)
{
  PetscScalar dd,dl,du;
  int         j;

  dd  = 4.0;
  dl  = -1.0;
  du  = -1.0;

  y[0] =  dd*x[0] + du*x[1];
  for (j=1;j<nx-1;j++)
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  y[nx-1] = dl*x[nx-2] + dd*x[nx-1];
}

/*
    Matrix-vector product subroutine for the 2D Laplacian.

    The matrix used is the 2 dimensional discrete Laplacian on unit square with
    zero Dirichlet boundary condition.

    Computes y <-- A*x, where A is the block tridiagonal matrix

                 | T -I          |
                 |-I  T -I       |
             A = |   -I  T       |
                 |        ...  -I|
                 |           -I T|

    The subroutine TV is called to compute y<--T*x.
 */
PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y)
{
  void              *ctx;
  int               nx,lo,i,j;
  const PetscScalar *px;
  PetscScalar       *py;

  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&ctx));
  nx = *(int*)ctx;
  PetscCall(VecGetArrayRead(x,&px));
  PetscCall(VecGetArray(y,&py));

  tv(nx,&px[0],&py[0]);
  for (i=0;i<nx;i++) py[i] -= px[nx+i];

  for (j=2;j<nx;j++) {
    lo = (j-1)*nx;
    tv(nx,&px[lo],&py[lo]);
    for (i=0;i<nx;i++) py[lo+i] -= px[lo-nx+i] + px[lo+nx+i];
  }

  lo = (nx-1)*nx;
  tv(nx,&px[lo],&py[lo]);
  for (i=0;i<nx;i++) py[lo+i] -= px[lo-nx+i];

  PetscCall(VecRestoreArrayRead(x,&px));
  PetscCall(VecRestoreArray(y,&py));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(diag,4.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Matrix-vector product subroutine for the Null matrix.
 */
PetscErrorCode MatMult_Zero(Mat A,Vec x,Vec y)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(y,0.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_Zero(Mat A,Vec diag)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(diag,0.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Matrix-vector product subroutine for the Identity matrix.
 */
PetscErrorCode MatMult_Identity(Mat A,Vec x,Vec y)
{
  PetscFunctionBeginUser;
  PetscCall(VecCopy(x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatGetDiagonal_Identity(Mat A,Vec diag)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(diag,1.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      suffix: 1
      args: -pep_type {{toar qarnoldi linear}} -pep_nev 4 -terse
      filter: grep -v Solution | sed -e "s/2.7996[1-8]i/2.79964i/g" | sed -e "s/2.7570[5-9]i/2.75708i/g" | sed -e "s/0.00000-2.79964i, 0.00000+2.79964i/0.00000+2.79964i, 0.00000-2.79964i/" | sed -e "s/0.00000-2.75708i, 0.00000+2.75708i/0.00000+2.75708i, 0.00000-2.75708i/"

TEST*/
