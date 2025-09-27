/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Estimates the 2-norm condition number of a matrix A, that is, the ratio of the largest to the smallest singular values of A.\n\n"
  "The matrix is a Grcar matrix.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n\n";

#include <slepcsvd.h>

/*
   This example computes the singular values of an nxn Grcar matrix,
   which is a nonsymmetric Toeplitz matrix:

              |  1  1  1  1               |
              | -1  1  1  1  1            |
              |    -1  1  1  1  1         |
              |       .  .  .  .  .       |
          A = |          .  .  .  .  .    |
              |            -1  1  1  1  1 |
              |               -1  1  1  1 |
              |                  -1  1  1 |
              |                     -1  1 |

 */

int main(int argc,char **argv)
{
  Mat            A;               /* Grcar matrix */
  SVD            svd;             /* singular value solver context */
  PetscInt       N=30,Istart,Iend,i,col[5],nconv1,nconv2;
  PetscScalar    value[] = { -1, 1, 1, 1, 1 };
  PetscReal      sigma_1,sigma_n;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nEstimate the condition number of a Grcar matrix, n=%" PetscInt_FMT "\n\n",N));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) PetscCall(MatSetValues(A,1,&i,PetscMin(4,N-i),col+1,value+1,INSERT_VALUES));
    else PetscCall(MatSetValues(A,1,&i,PetscMin(5,N-i+1),col,value,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create the singular value solver and set the solution method
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create singular value context
  */
  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));

  /*
     Set operator
  */
  PetscCall(SVDSetOperators(svd,A,NULL));

  /*
     Set solver parameters at runtime
  */
  PetscCall(SVDSetFromOptions(svd));
  PetscCall(SVDSetDimensions(svd,1,PETSC_DETERMINE,PETSC_DETERMINE));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the singular value problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     First request a singular value from one end of the spectrum
  */
  PetscCall(SVDSetWhichSingularTriplets(svd,SVD_LARGEST));
  PetscCall(SVDSolve(svd));
  /*
     Get number of converged singular values
  */
  PetscCall(SVDGetConverged(svd,&nconv1));
  /*
     Get converged singular values: largest singular value is stored in sigma_1.
     In this example, we are not interested in the singular vectors
  */
  if (nconv1 > 0) PetscCall(SVDGetSingularTriplet(svd,0,&sigma_1,NULL,NULL));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD," Unable to compute large singular value!\n\n"));

  /*
     Request a singular value from the other end of the spectrum
  */
  PetscCall(SVDSetWhichSingularTriplets(svd,SVD_SMALLEST));
  PetscCall(SVDSolve(svd));
  /*
     Get number of converged singular triplets
  */
  PetscCall(SVDGetConverged(svd,&nconv2));
  /*
     Get converged singular values: smallest singular value is stored in sigma_n.
     As before, we are not interested in the singular vectors
  */
  if (nconv2 > 0) PetscCall(SVDGetSingularTriplet(svd,0,&sigma_n,NULL,NULL));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD," Unable to compute small singular value!\n\n"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (nconv1 > 0 && nconv2 > 0) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Computed singular values: sigma_1=%.4f, sigma_n=%.4f\n",(double)sigma_1,(double)sigma_n));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Estimated condition number: sigma_1/sigma_n=%.4f\n\n",(double)(sigma_1/sigma_n)));
  }

  /*
     Free work space
  */
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1

TEST*/
