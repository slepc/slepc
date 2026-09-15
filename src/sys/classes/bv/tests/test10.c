/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test split reductions in BV.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v,w,y,z,zsplit;
  BV             X;
  PetscInt       i,j,n=10,k=5;
  PetscScalar    *zarray;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCheck(k>2,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Should specify at least k=3 columns");
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"BV split ops (%" PetscInt_FMT " columns of dimension %" PetscInt_FMT ").\n",k,n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));
  PetscCall(VecDuplicate(t,&v));
  PetscCall(VecSet(v,1.0));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetFromOptions(X));

  /* Fill X entries */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&w));
    PetscCall(VecSet(w,0.0));
    for (i=0;i<4;i++) {
      if (i+j<n) PetscCall(VecSetValue(w,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(w));
    PetscCall(VecAssemblyEnd(w));
    PetscCall(BVRestoreColumn(X,j,&w));
  }

  /* Use regular operations */
  PetscCall(VecCreateSeq(PETSC_COMM_SELF,k+6,&z));
  PetscCall(VecGetArray(z,&zarray));
  PetscCall(BVGetColumn(X,0,&w));
  PetscCall(VecDot(w,v,zarray));
  PetscCall(BVRestoreColumn(X,0,&w));
  PetscCall(BVDotVec(X,v,zarray+1));
  PetscCall(BVDotColumn(X,2,zarray+1+k));

  PetscCall(BVGetColumn(X,1,&y));
  PetscCall(VecNorm(y,NORM_2,&nrm));
  PetscCall(BVRestoreColumn(X,1,&y));
  zarray[k+3] = nrm;
  PetscCall(BVNormVec(X,v,NORM_2,&nrm));
  zarray[k+4] = nrm;
  PetscCall(BVNormColumn(X,0,NORM_2,&nrm));
  zarray[k+5] = nrm;
  PetscCall(VecRestoreArray(z,&zarray));

  /* Use split operations */
  PetscCall(VecCreateSeq(PETSC_COMM_SELF,k+6,&zsplit));
  PetscCall(VecGetArray(zsplit,&zarray));
  PetscCall(BVGetColumn(X,0,&w));
  PetscCall(VecDotBegin(w,v,zarray));
  PetscCall(BVDotVecBegin(X,v,zarray+1));
  PetscCall(BVDotColumnBegin(X,2,zarray+1+k));
  PetscCall(VecDotEnd(w,v,zarray));
  PetscCall(BVRestoreColumn(X,0,&w));
  PetscCall(BVDotVecEnd(X,v,zarray+1));
  PetscCall(BVDotColumnEnd(X,2,zarray+1+k));

  PetscCall(BVGetColumn(X,1,&y));
  PetscCall(VecNormBegin(y,NORM_2,&nrm));
  PetscCall(BVNormVecBegin(X,v,NORM_2,&nrm));
  PetscCall(BVNormColumnBegin(X,0,NORM_2,&nrm));
  PetscCall(VecNormEnd(y,NORM_2,&nrm));
  PetscCall(BVRestoreColumn(X,1,&y));
  zarray[k+3] = nrm;
  PetscCall(BVNormVecEnd(X,v,NORM_2,&nrm));
  zarray[k+4] = nrm;
  PetscCall(BVNormColumnEnd(X,0,NORM_2,&nrm));
  zarray[k+5] = nrm;
  PetscCall(VecRestoreArray(zsplit,&zarray));

  /* Show difference */
  PetscCall(VecAXPY(z,-1.0,zsplit));
  PetscCall(VecNorm(z,NORM_1,&nrm));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%g\n",(double)nrm));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,NULL));

  PetscCall(BVDestroy(&X));
  PetscCall(VecDestroy(&t));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&z));
  PetscCall(VecDestroy(&zsplit));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: 2
      output_file: output/test10_1.out
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

TEST*/
