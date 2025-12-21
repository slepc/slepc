/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test block orthogonalization on a rank-deficient BV.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  BV             X,Z;
  Mat            M,R;
  Vec            v,w,t;
  PetscInt       i,j,n=20,k=8;
  PetscViewer    view;
  PetscBool      verbose;
  PetscReal      norm;
  PetscScalar    alpha;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV block orthogonalization (length %" PetscInt_FMT ", k=%" PetscInt_FMT ").\n",n,k));
  PetscCheck(k>5,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"k must be at least 6");

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

  /* Fill X entries (first half) */
  for (j=0;j<k/2;j++) {
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

  /* make middle column linearly dependent wrt columns 0 and 1 */
  PetscCall(BVCopyColumn(X,0,j));
  PetscCall(BVGetColumn(X,j,&v));
  PetscCall(BVGetColumn(X,1,&w));
  PetscCall(VecAXPY(v,0.5,w));
  PetscCall(BVRestoreColumn(X,1,&w));
  PetscCall(BVRestoreColumn(X,j,&v));
  j++;

  /* Fill X entries (second half) */
  for (;j<k-1;j++) {
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

  /* make middle column linearly dependent wrt columns 1 and k/2+1 */
  PetscCall(BVCopyColumn(X,1,j));
  PetscCall(BVGetColumn(X,j,&v));
  PetscCall(BVGetColumn(X,k/2+1,&w));
  PetscCall(VecAXPY(v,-1.2,w));
  PetscCall(BVRestoreColumn(X,k/2+1,&w));
  PetscCall(BVRestoreColumn(X,j,&v));

  if (verbose) PetscCall(BVView(X,view));

  /* Create a copy on Z */
  PetscCall(BVDuplicate(X,&Z));
  PetscCall(PetscObjectSetName((PetscObject)Z,"Z"));
  PetscCall(BVCopy(X,Z));

  /* Test BVOrthogonalize */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&R));
  PetscCall(PetscObjectSetName((PetscObject)R,"R"));
  PetscCall(BVOrthogonalize(X,R));
  if (verbose) {
    PetscCall(BVView(X,view));
    PetscCall(MatView(R,view));
  }

  /* Check orthogonality */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));
  PetscCall(MatShift(M,1.0));   /* set leading part to identity */
  PetscCall(BVDot(X,X,M));
  PetscCall(MatShift(M,-1.0));
  PetscCall(MatNorm(M,NORM_1,&norm));
  if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)norm));

  /* Check residual */
  PetscCall(BVMult(Z,-1.0,1.0,X,R));
  PetscCall(BVNorm(Z,NORM_FROBENIUS,&norm));
  if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Residual ||X-QR|| < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Residual ||X-QR||: %g\n",(double)norm));

  PetscCall(MatDestroy(&R));
  PetscCall(MatDestroy(&M));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Z));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: 1
      args: -bv_orthog_block gs
      output_file: output/test12_1.out
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

TEST*/
