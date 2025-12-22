/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BVGetSplitRows().\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Mat            M,T,D,block[4];
  Vec            t,v,v1,v2,w,*C;
  VecType        vtype;
  BV             X,U,L;
  IS             is[2];
  PetscInt       i,j,n=10,k=5,l=3,nc=0,rstart,rend;
  PetscViewer    view;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nc",&nc,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BVGetSplitRows (length %" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT ", nc=%" PetscInt_FMT ").\n",n,l,k,nc));

  /* Create Nest matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&T));
  PetscCall(MatSetSizes(T,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(T));
  PetscCall(MatGetOwnershipRange(T,&rstart,&rend));
  for (i=rstart;i<rend;i++) {
    if (i>0) PetscCall(MatSetValue(T,i,i-1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(T,i,i,2.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(T,i,i+1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY));
  PetscCall(MatGetVecType(T,&vtype));
  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,1.0,&D));
  PetscCall(MatSetVecType(D,vtype));

  block[0] = T;
  block[1] = block[2] = NULL;
  block[3] = D;
  PetscCall(MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,block,&M));
  PetscCall(MatDestroy(&T));
  PetscCall(MatDestroy(&D));
  PetscCall(MatNestGetISs(M,is,NULL));

  PetscCall(MatView(M,NULL));

  /* Create template vector */
  PetscCall(MatCreateVecs(M,&t,NULL));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));

  /* Generate constraints and attach them to X */
  if (nc>0) {
    PetscCall(VecDuplicateVecs(t,nc,&C));
    for (j=0;j<nc;j++) {
      for (i=0;i<=j;i++) PetscCall(VecSetValue(C[j],nc-i+1,1.0,INSERT_VALUES));
      PetscCall(VecAssemblyBegin(C[j]));
      PetscCall(VecAssemblyEnd(C[j]));
    }
    PetscCall(BVInsertConstraints(X,&nc,C));
    PetscCall(VecDestroyVecs(nc,&C));
  }

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,0.0));
    PetscCall(VecGetSubVector(v,is[0],&v1));
    PetscCall(VecGetSubVector(v,is[1],&v2));
    for (i=0;i<4;i++) {
      if (i+j>=rstart && i+j<rend) {
        PetscCall(VecSetValue(v1,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
        PetscCall(VecSetValue(v2,i+j,(PetscScalar)(-i+2*(j+1)),INSERT_VALUES));
      }
    }
    PetscCall(VecAssemblyBegin(v1));
    PetscCall(VecAssemblyBegin(v2));
    PetscCall(VecAssemblyEnd(v1));
    PetscCall(VecAssemblyEnd(v2));
    PetscCall(VecRestoreSubVector(v,is[0],&v1));
    PetscCall(VecRestoreSubVector(v,is[1],&v2));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Get split BVs */
  PetscCall(BVGetSplitRows(X,is[0],is[1],&U,&L));
  PetscCall(PetscObjectSetName((PetscObject)U,"U"));
  PetscCall(PetscObjectSetName((PetscObject)L,"L"));

  if (verbose) {
    PetscCall(BVView(U,view));
    PetscCall(BVView(L,view));
  }

  /* Copy l-th column of U to first column of L */
  PetscCall(BVGetColumn(U,l,&v));
  PetscCall(BVGetColumn(L,0,&w));
  PetscCall(VecCopy(v,w));
  PetscCall(BVRestoreColumn(U,l,&v));
  PetscCall(BVRestoreColumn(L,0,&w));

  /* Finished using the split BVs */
  PetscCall(BVRestoreSplitRows(X,is[0],is[1],&U,&L));
  if (verbose) PetscCall(BVView(X,view));

  /* Check: print bottom part of first column */
  PetscCall(BVGetColumn(X,0,&v));
  PetscCall(VecGetSubVector(v,is[1],&v2));
  PetscCall(VecView(v2,NULL));
  PetscCall(VecRestoreSubVector(v,is[1],&v2));
  PetscCall(BVRestoreColumn(X,0,&v));

  PetscCall(BVDestroy(&X));
  PetscCall(VecDestroy(&t));
  PetscCall(MatDestroy(&M));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: {{1 2}}
      output_file: output/test19_1.out
      filter: grep -v Process | grep -v Object | sed -e 's/mpi/seq/' | sed -e 's/seqcuda/seq/' | sed -e 's/seqaijcusparse/seqaij/' | sed -e 's/seqhip/seq/' | sed -e 's/seqaijhipsparse/seqaij/' | sed -e 's/nc=2/nc=0/'
      test:
         suffix: 1
         args: -nc {{0 2}} -bv_type {{svec mat}}
      test:
         suffix: 1_cuda
         args: -nc {{0 2}} -bv_type {{svec mat}} -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 1_hip
         args: -nc {{0 2}} -bv_type {{svec mat}} -mat_type aijhipsparse
         requires: hip

TEST*/
