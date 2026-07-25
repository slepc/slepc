/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BVNormalize().\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  BV             X,Y,Z;
  Mat            B;
  Vec            v,t;
  PetscInt       i,j,n=20,k=8,l=3,Istart,Iend;
  PetscViewer    view;
  PetscBool      verbose;
  PetscReal      norm,error;
  PetscScalar    alpha;
#if !PetscDefined(USE_COMPLEX)
  PetscScalar    *eigi;
  PetscRandom    rand;
  PetscReal      normr,normi;
  Vec            vi;
#endif

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV normalization with %" PetscInt_FMT " columns of length %" PetscInt_FMT ".\n",k,n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=0;i<=n/2;i++) {
      if (i+j<n) {
        alpha = (3.0*i+j-2)/(2*(i+j+1));
        PetscCall(VecSetValue(v,i+j,alpha,INSERT_VALUES));
      }
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Create copies on Y and Z */
  PetscCall(BVDuplicate(X,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVCopy(X,Y));
  PetscCall(BVDuplicate(X,&Z));
  PetscCall(PetscObjectSetName((PetscObject)Z,"Z"));
  PetscCall(BVCopy(X,Z));
  PetscCall(BVSetActiveColumns(X,l,k));
  PetscCall(BVSetActiveColumns(Y,l,k));
  PetscCall(BVSetActiveColumns(Z,l,k));

  /* Test BVNormalize */
  PetscCall(BVNormalize(X,NULL));
  if (verbose) PetscCall(BVView(X,view));

  /* Check unit norm of columns */
  error = 0.0;
  for (j=l;j<k;j++) {
    PetscCall(BVNormColumn(X,j,NORM_2,&norm));
    error = PetscMax(error,PetscAbsReal(norm-PetscRealConstant(1.0)));
  }
  if (error<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from normalized vectors < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from normalized vectors: %g\n",(double)norm));

  /* Create inner product matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(B));
  PetscCall(PetscObjectSetName((PetscObject)B,"B"));

  PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(B,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(B,i,i+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(B,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  if (verbose) PetscCall(MatView(B,view));

  /* Test BVNormalize with B-norm */
  PetscCall(BVSetMatrix(Y,B,PETSC_FALSE));
  PetscCall(BVNormalize(Y,NULL));
  if (verbose) PetscCall(BVView(Y,view));

  /* Check unit B-norm of columns */
  error = 0.0;
  for (j=l;j<k;j++) {
    PetscCall(BVNormColumn(Y,j,NORM_2,&norm));
    error = PetscMax(error,PetscAbsReal(norm-PetscRealConstant(1.0)));
  }
  if (error<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from B-normalized vectors < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from B-normalized vectors: %g\n",(double)norm));

#if !PetscDefined(USE_COMPLEX)
  /* fill imaginary parts */
  PetscCall(PetscCalloc1(k,&eigi));
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&rand));
  PetscCall(PetscRandomSetFromOptions(rand));
  for (j=l+1;j<k-1;j+=5) {
    PetscCall(PetscRandomGetValue(rand,&alpha));
    eigi[j]   =  alpha;
    eigi[j+1] = -alpha;
  }
  PetscCall(PetscRandomDestroy(&rand));
  if (verbose) {
    PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF,1,k,eigi,&v));
    PetscCall(VecView(v,view));
    PetscCall(VecDestroy(&v));
  }

  /* Test BVNormalize with complex conjugate columns */
  PetscCall(BVNormalize(Z,eigi));
  if (verbose) PetscCall(BVView(Z,view));

  /* Check unit norm of (complex conjugate) columns */
  error = 0.0;
  for (j=l;j<k;j++) {
    if (eigi[j]) {
      PetscCall(BVGetColumn(Z,j,&v));
      PetscCall(BVGetColumn(Z,j+1,&vi));
      PetscCall(VecNormBegin(v,NORM_2,&normr));
      PetscCall(VecNormBegin(vi,NORM_2,&normi));
      PetscCall(VecNormEnd(v,NORM_2,&normr));
      PetscCall(VecNormEnd(vi,NORM_2,&normi));
      PetscCall(BVRestoreColumn(Z,j+1,&vi));
      PetscCall(BVRestoreColumn(Z,j,&v));
      norm = SlepcAbsEigenvalue(normr,normi);
      j++;
    } else {
      PetscCall(BVGetColumn(Z,j,&v));
      PetscCall(VecNorm(v,NORM_2,&norm));
      PetscCall(BVRestoreColumn(Z,j,&v));
    }
    error = PetscMax(error,PetscAbsReal(norm-PetscRealConstant(1.0)));
  }
  if (error<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from normalized conjugate vectors < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Deviation from normalized conjugate vectors: %g\n",(double)norm));
  PetscCall(PetscFree(eigi));
#endif

  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(BVDestroy(&Z));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -n 250 -l 6 -k 15
      nsize: {{1 2}}
      requires: !complex
      output_file: output/test18_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}}
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip

   testset:
      args: -n 250 -l 6 -k 15
      nsize: {{1 2}}
      requires: complex
      output_file: output/test18_1_complex.out
      test:
         suffix: 1_complex
         args: -bv_type {{vecs contiguous svec mat}}
      test:
         suffix: 1_cuda_complex
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 1_hip_complex
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip

TEST*/
