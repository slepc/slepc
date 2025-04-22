/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV orthogonalization with selected columns.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  BV             X;
  Vec            v,t,z;
  PetscInt       i,j,n=20,k=8;
  PetscViewer    view;
  PetscBool      verbose,*which;
  PetscReal      norm;
  PetscScalar    alpha,*pz;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV orthogonalization with selected columns of length %" PetscInt_FMT ".\n",n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetOrthogonalization(X,BV_ORTHOG_MGS,BV_ORTHOG_REFINE_IFNEEDED,PETSC_CURRENT,BV_ORTHOG_BLOCK_GS));
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

  /* Orthonormalize first k-1 columns */
  for (j=0;j<k-1;j++) {
    PetscCall(BVOrthogonalizeColumn(X,j,NULL,&norm,NULL));
    alpha = 1.0/norm;
    PetscCall(BVScaleColumn(X,j,alpha));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Select odd columns and orthogonalize last column against those only */
  PetscCall(PetscMalloc1(k,&which));
  for (i=0;i<k;i++) which[i] = (i%2)? PETSC_TRUE: PETSC_FALSE;
  PetscCall(BVOrthogonalizeSomeColumn(X,k-1,which,NULL,NULL,NULL));
  PetscCall(PetscFree(which));
  if (verbose) PetscCall(BVView(X,view));

  /* Check orthogonality */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Orthogonalization coefficients:\n"));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF,k-1,&z));
  PetscCall(PetscObjectSetName((PetscObject)z,"z"));
  PetscCall(VecGetArray(z,&pz));
  PetscCall(BVDotColumn(X,k-1,pz));
  for (i=0;i<k-1;i++) {
    if (PetscAbsScalar(pz[i])<5.0*PETSC_MACHINE_EPSILON) pz[i]=0.0;
  }
  PetscCall(VecRestoreArray(z,&pz));
  PetscCall(VecView(z,view));
  PetscCall(VecDestroy(&z));

  PetscCall(BVDestroy(&X));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test8_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip
      test:
         suffix: 2
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_refine never
         requires: !single
      test:
         suffix: 3
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_refine always

TEST*/
