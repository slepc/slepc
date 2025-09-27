/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves the same eigenproblem as in example ex2, but using a shell matrix.\n\n"
  "The problem is a standard symmetric eigenproblem corresponding to the 2-D Laplacian operator.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in both x and y dimensions.\n\n";

#include <slepceps.h>
#include <petscblaslapack.h>

/*
   User-defined routines
*/
PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON;
  PetscMPIInt    size;
  PetscInt       N,n=10,nev;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size==1,PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only");

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  N = n*n;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n2-D Laplacian Eigenproblem (matrix-free version), N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid)\n\n",N,n,n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,&n,&A));
  PetscCall(MatShellSetOperation(A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_Laplacian2D));
  PetscCall(MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(PetscErrorCodeFn*)MatMult_Laplacian2D));
  PetscCall(MatShellSetOperation(A,MATOP_GET_DIAGONAL,(PetscErrorCodeFn*)MatGetDiagonal_Laplacian2D));

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
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
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

/*TEST

   testset:
      args: -n 20 -eps_nev 4 -eps_ncv 11 -eps_max_it 40000
      requires: !single
      output_file: output/test8_1.out
      test:
         suffix: 1
         args: -eps_type {{power subspace arnoldi}}
      test:
         suffix: 1_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog local
      test:
         suffix: 1_lapack
         args: -eps_type lapack
         timeoutfactor: 2
      test:
         suffix: 1_elemental
         args: -eps_type elemental
         requires: elemental
      test:
         suffix: 1_krylovschur_vecs
         args: -bv_type vecs -bv_orthog_refine always -eps_ncv 10 -vec_mdot_use_gemv 0
      test:
         suffix: 1_jd
         args: -eps_type jd -eps_jd_blocksize 3
         requires: !__float128
      test:
         suffix: 1_gd
         args: -eps_type gd -eps_gd_blocksize 3 -eps_tol 1e-8
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 1_primme
         args: -eps_type primme -eps_conv_abs -eps_largest_magnitude
         requires: primme

   testset:
      args: -eps_nev 4 -eps_smallest_real -eps_max_it 600
      output_file: output/test8_2.out
      test:
         suffix: 2
         args: -eps_type {{rqcg lobpcg}}
      test:
         suffix: 2_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog local
      test:
         suffix: 2_arpack
         args: -eps_type arpack -eps_ncv 6
         requires: arpack !single
      test:
         suffix: 2_primme
         args: -eps_type primme -eps_conv_abs -eps_primme_method lobpcg_orthobasisw -eps_ncv 24
         requires: primme

   testset:
      args: -eps_nev 12 -eps_mpd 9 -eps_smallest_real -eps_max_it 1000
      output_file: output/test8_3.out
      test:
         suffix: 3_rqcg
         args: -eps_type rqcg
      test:
         suffix: 3_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog local
      test:
         suffix: 3_lobpcg
         args: -eps_type lobpcg -eps_lobpcg_blocksize 3 -eps_lobpcg_locking 0 -st_ksp_type preonly -st_pc_type jacobi -eps_lobpcg_checkprecond
         requires: !__float128
      test:
         suffix: 3_lobpcg_quad
         args: -eps_type lobpcg -eps_lobpcg_blocksize 3 -eps_lobpcg_locking 0 -st_ksp_type preonly -st_pc_type jacobi -eps_tol 1e-25
         requires: __float128
TEST*/
