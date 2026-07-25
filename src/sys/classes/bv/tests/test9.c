/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV matrix projection.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v;
  Mat            B,G,H0,H1;
  BV             X,Y,Z;
  PetscInt       i,j,n=20,kx=6,lx=3,ky=5,ly=2,Istart,Iend,col[5];
  PetscScalar    alpha,value[] = { -1, 1, 1, 1, 1 };
  PetscViewer    view;
  PetscReal      norm;
  PetscBool      verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-kx",&kx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-lx",&lx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ky",&ky,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ly",&ly,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV projection (n=%" PetscInt_FMT ").\n",n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X has %" PetscInt_FMT " active columns (%" PetscInt_FMT " leading columns).\n",kx,lx));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Y has %" PetscInt_FMT " active columns (%" PetscInt_FMT " leading columns).\n",ky,ly));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Create non-symmetric matrix G (Toeplitz) */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&G));
  PetscCall(MatSetSizes(G,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(G));
  PetscCall(PetscObjectSetName((PetscObject)G,"G"));

  PetscCall(MatGetOwnershipRange(G,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) PetscCall(MatSetValues(G,1,&i,PetscMin(4,n-i),col+1,value+1,INSERT_VALUES));
    else PetscCall(MatSetValues(G,1,&i,PetscMin(5,n-i+1),col,value,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(G,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(G,MAT_FINAL_ASSEMBLY));
  if (verbose) PetscCall(MatView(G,view));

  /* Create symmetric matrix B (1-D Laplacian) */
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
  PetscCall(MatCreateVecs(B,&t,NULL));
  if (verbose) PetscCall(MatView(B,view));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,kx+2));  /* two extra columns to test active columns */
  PetscCall(BVSetFromOptions(X));

  /* Fill X entries */
  for (j=0;j<kx+2;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=0;i<4;i++) {
      if (i+j<n) {
#if PetscDefined(USE_COMPLEX)
        alpha = PetscCMPLX((PetscReal)(3*i+j-2),(PetscReal)(2*i));
#else
        alpha = (PetscReal)(3*i+j-2);
#endif
        PetscCall(VecSetValue(v,i+j,alpha,INSERT_VALUES));
      }
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) PetscCall(BVView(X,view));

  /* Duplicate BV object and store Z=G*X */
  PetscCall(BVDuplicate(X,&Z));
  PetscCall(PetscObjectSetName((PetscObject)Z,"Z"));
  PetscCall(BVSetActiveColumns(X,0,kx));
  PetscCall(BVSetActiveColumns(Z,0,kx));
  PetscCall(BVMatMult(X,G,Z));
  PetscCall(BVSetActiveColumns(X,lx,kx));
  PetscCall(BVSetActiveColumns(Z,lx,kx));

  /* Create BV object Y */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVSetSizesFromVec(Y,t,ky+1));
  PetscCall(BVSetFromOptions(Y));
  PetscCall(BVSetActiveColumns(Y,ly,ky));

  /* Fill Y entries */
  for (j=0;j<ky+1;j++) {
    PetscCall(BVGetColumn(Y,j,&v));
#if PetscDefined(USE_COMPLEX)
    alpha = PetscCMPLX((PetscReal)(j+1)/4.0,-(PetscReal)j);
#else
    alpha = (PetscReal)(j+1)/4.0;
#endif
    PetscCall(VecSet(v,(PetscScalar)(j+1)/4.0));
    PetscCall(BVRestoreColumn(Y,j,&v));
  }
  if (verbose) PetscCall(BVView(Y,view));

  /* Test BVMatProject for non-symmetric matrix G */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H0));
  PetscCall(PetscObjectSetName((PetscObject)H0,"H0"));
  PetscCall(BVMatProject(X,G,Y,H0));
  if (verbose) PetscCall(MatView(H0,view));

  /* Test BVMatProject with previously stored G*X */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H1));
  PetscCall(PetscObjectSetName((PetscObject)H1,"H1"));
  PetscCall(BVMatProject(Z,NULL,Y,H1));
  if (verbose) PetscCall(MatView(H1,view));

  /* Check that H0 and H1 are equal */
  PetscCall(MatAXPY(H0,-1.0,H1,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(H0,NORM_1,&norm));
  if (norm<10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||H0-H1|| < 10*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||H0-H1||=%g\n",(double)norm));
  PetscCall(MatDestroy(&H0));
  PetscCall(MatDestroy(&H1));

  /* Test BVMatProject for symmetric matrix B with orthogonal projection */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,kx,kx,NULL,&H0));
  PetscCall(PetscObjectSetName((PetscObject)H0,"H0"));
  PetscCall(BVMatProject(X,B,X,H0));
  if (verbose) PetscCall(MatView(H0,view));

  /* Repeat previous test with symmetry flag set */
  PetscCall(MatSetOption(B,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,kx,kx,NULL,&H1));
  PetscCall(PetscObjectSetName((PetscObject)H1,"H1"));
  PetscCall(BVMatProject(X,B,X,H1));
  if (verbose) PetscCall(MatView(H1,view));

  /* Check that H0 and H1 are equal */
  PetscCall(MatAXPY(H0,-1.0,H1,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(H0,NORM_1,&norm));
  if (norm<10*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||H0-H1|| < 10*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||H0-H1||=%g\n",(double)norm));
  PetscCall(MatDestroy(&H0));
  PetscCall(MatDestroy(&H1));

  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(BVDestroy(&Z));
  PetscCall(MatDestroy(&B));
  PetscCall(MatDestroy(&G));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test9_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 1_svec_vecs
         args: -bv_type svec -bv_matmult vecs
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -mat_type aijhipsparse
         requires: hip
      test:
         suffix: 2
         nsize: 2
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 2_svec_vecs
         nsize: 2
         args: -bv_type svec -bv_matmult vecs

TEST*/
