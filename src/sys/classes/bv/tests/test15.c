/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BVGetSplit().\n\n";

#include <slepcbv.h>

/*
   Print the first row of a BV
 */
PetscErrorCode PrintFirstRow(BV X)
{
  PetscMPIInt       rank;
  PetscInt          i,ld,k,nc;
  const PetscScalar *pX;
  const char        *name;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)X),&rank));
  if (!rank) {
    PetscCall(BVGetActiveColumns(X,NULL,&k));
    PetscCall(BVGetLeadingDimension(X,&ld));
    PetscCall(BVGetNumConstraints(X,&nc));
    PetscCall(PetscObjectGetName((PetscObject)X,&name));
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)X),"First row of %s =\n",name));
    PetscCall(BVGetArrayRead(X,&pX));
    for (i=0;i<nc+k;i++) PetscCall(PetscPrintf(PetscObjectComm((PetscObject)X),"%g ",(double)PetscRealPart(pX[i*ld])));
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)X),"\n"));
    PetscCall(BVRestoreArrayRead(X,&pX));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Vec            t,v,*C;
  BV             X,L,R;
  PetscInt       i,j,n=10,k=5,l=3,nc=0;
  PetscReal      norm;
  PetscScalar    alpha;
  PetscViewer    view;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nc",&nc,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BVGetSplit (length %" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT ", nc=%" PetscInt_FMT ").\n",n,l,k,nc));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

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
    for (i=0;i<4;i++) {
      if (i+j<n) PetscCall(VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Get split BVs */
  PetscCall(BVSetActiveColumns(X,l,k));
  PetscCall(BVGetSplit(X,&L,&R));
  PetscCall(PetscObjectSetName((PetscObject)L,"L"));
  PetscCall(PetscObjectSetName((PetscObject)R,"R"));

  if (verbose) {
    PetscCall(BVView(L,view));
    PetscCall(BVView(R,view));
  }

  /* Modify first column of R */
  PetscCall(BVGetColumn(R,0,&v));
  PetscCall(VecSet(v,-1.0));
  PetscCall(BVRestoreColumn(R,0,&v));

  /* Finished using the split BVs */
  PetscCall(BVRestoreSplit(X,&L,&R));
  PetscCall(PrintFirstRow(X));
  if (verbose) PetscCall(BVView(X,view));

  /* Get the left split BV only */
  PetscCall(BVGetSplit(X,&L,NULL));
  for (j=0;j<l;j++) {
    PetscCall(BVOrthogonalizeColumn(L,j,NULL,&norm,NULL));
    alpha = 1.0/norm;
    PetscCall(BVScaleColumn(L,j,alpha));
  }
  PetscCall(BVRestoreSplit(X,&L,NULL));
  PetscCall(PrintFirstRow(X));
  if (verbose) PetscCall(BVView(X,view));

  /* Now get the right split BV after changing the number of leading columns */
  PetscCall(BVSetActiveColumns(X,l-1,k));
  PetscCall(BVGetSplit(X,NULL,&R));
  PetscCall(BVGetColumn(R,0,&v));
  PetscCall(BVInsertVec(X,0,v));
  PetscCall(BVRestoreColumn(R,0,&v));
  PetscCall(BVRestoreSplit(X,NULL,&R));
  PetscCall(PrintFirstRow(X));
  if (verbose) PetscCall(BVView(X,view));

  PetscCall(BVDestroy(&X));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: 2
      output_file: output/test15_1.out
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
         suffix: 1_kokkos
         args: -bv_type mat -vec_type kokkos
         requires: hip kokkos_kernels

   testset:
      nsize: 2
      output_file: output/test15_2.out
      test:
         suffix: 2
         args: -nc 2 -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 2_cuda
         args: -nc 2 -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 2_hip
         args: -nc 2 -bv_type {{svec mat}} -vec_type hip
         requires: hip
      test:
         suffix: 2_kokkos
         args: -nc 2 -bv_type mat -vec_type kokkos
         requires: hip kokkos_kernels

TEST*/
