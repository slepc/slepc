/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix function evaluation via diagonalization.\n\n";

#include <slepcfn.h>

int main(int argc,char **argv)
{
  FN             fn;
  Mat            A,F,G;
  PetscInt       i,j,n=10;
  PetscReal      nrm;
  PetscScalar    *As,alpha,beta;
  PetscViewer    viewer;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix function of symmetric/Hermitian matrix, n=%" PetscInt_FMT ".\n",n));

  /* Create function object */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNEXP));   /* default to exponential */
#if PetscDefined(USE_COMPLEX)
  alpha = PetscCMPLX(0.3,0.8);
  beta  = PetscCMPLX(1.1,-0.1);
#else
  alpha = 0.3;
  beta  = 1.1;
#endif
  PetscCall(FNSetScale(fn,alpha,beta));
  PetscCall(FNSetFromOptions(fn));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Create a symmetric/Hermitian Toeplitz matrix */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));
  PetscCall(MatDenseGetArray(A,&As));
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
#if PetscDefined(USE_COMPLEX)
      As[i+(i+j)*n]=PetscCMPLX(1.0,0.1); As[(i+j)+i*n]=PetscCMPLX(1.0,-0.1);
#else
      As[i+(i+j)*n]=0.5; As[(i+j)+i*n]=0.5;
#endif
    }
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    PetscCall(MatView(A,viewer));
  }

  /* compute matrix function */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&F));
  PetscCall(PetscObjectSetName((PetscObject)F,"F"));
  PetscCall(FNEvaluateFunctionMat(fn,A,F));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) - - - - - - -\n"));
    PetscCall(MatView(F,viewer));
  }

  /* Repeat with MAT_HERMITIAN flag set */
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&G));
  PetscCall(PetscObjectSetName((PetscObject)G,"G"));
  PetscCall(FNEvaluateFunctionMat(fn,A,G));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) symm - - - - - - -\n"));
    PetscCall(MatView(G,viewer));
  }

  /* compare the two results */
  PetscCall(MatAXPY(F,-1.0,G,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(F,NORM_FROBENIUS,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of F-G is %g\n",(double)nrm));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed results match.\n"));

  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&F));
  PetscCall(MatDestroy(&G));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -fn_type {{exp sqrt}shared output}
      output_file: output/test12_1.out

   test:
      suffix: 1_rational
      nsize: 1
      args: -fn_type rational -fn_rational_numerator 2,-1.5 -fn_rational_denominator 1,0.8
      output_file: output/test12_1.out

TEST*/
