/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test multiplication of a Mat times a BV.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,r,v;
  Mat            B,Ymat;
  BV             X,Y,Z=NULL,Zcopy=NULL;
  PetscInt       i,j,m=10,n,k=5,rep=1,Istart,Iend,ld;
  PetscScalar    *pZ;
  PetscReal      norm;
  PetscViewer    view;
  PetscBool      flg,verbose,fromfile;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  BVMatMultType  vmm;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-rep",&rep,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-file",filename,sizeof(filename),&fromfile));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
  PetscCall(PetscObjectSetName((PetscObject)B,"B"));
  if (fromfile) {
#if PetscDefined(USE_COMPLEX)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrix from a binary file...\n"));
#else
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrix from a binary file...\n"));
#endif
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    PetscCall(MatSetFromOptions(B));
    PetscCall(MatLoad(B,viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(MatGetSize(B,&m,&n));
    PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
  } else {
    /* Create 1-D Laplacian matrix */
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg));
    if (!flg) n = m;
    PetscCall(MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,m,n));
    PetscCall(MatSetFromOptions(B));
    PetscCall(MatGetOwnershipRange(B,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {
      if (i>0 && i-1<n) PetscCall(MatSetValue(B,i,i-1,-1.0,INSERT_VALUES));
      if (i+1<n) PetscCall(MatSetValue(B,i,i+1,-1.0,INSERT_VALUES));
      if (i<n) PetscCall(MatSetValue(B,i,i,2.0,INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  }

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BVMatMult (m=%" PetscInt_FMT ", n=%" PetscInt_FMT ", k=%" PetscInt_FMT ").\n",m,n,k));
  PetscCall(MatCreateVecs(B,&t,&r));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,k));
  PetscCall(BVSetMatMultMethod(X,BV_MATMULT_VECS));
  PetscCall(BVSetFromOptions(X));
  PetscCall(BVGetMatMultMethod(X,&vmm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Using method: %s\n",BVMatMultTypes[vmm]));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  for (j=0;j<k;j++) {
    PetscCall(BVGetColumn(X,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=Istart;i<PetscMin(j+1,Iend);i++) PetscCall(VecSetValue(v,i,1.0,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(X,j,&v));
  }
  if (verbose) {
    PetscCall(MatView(B,view));
    PetscCall(BVView(X,view));
  }

  /* Create BV object Y */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVSetSizesFromVec(Y,r,k+4));
  PetscCall(BVSetMatMultMethod(Y,BV_MATMULT_VECS));
  PetscCall(BVSetFromOptions(Y));
  PetscCall(BVSetActiveColumns(Y,2,k+2));

  /* Test BVMatMult */
  for (i=0;i<rep;i++) PetscCall(BVMatMult(X,B,Y));
  if (verbose) PetscCall(BVView(Y,view));

  if (fromfile) {
    /* Test BVMatMultTranspose */
    PetscCall(BVDuplicate(X,&Z));
    PetscCall(BVSetRandom(Z));
    for (i=0;i<rep;i++) PetscCall(BVMatMultTranspose(Z,B,Y));
    if (verbose) {
      PetscCall(BVView(Z,view));
      PetscCall(BVView(Y,view));
    }
    PetscCall(BVDestroy(&Z));
    PetscCall(BVMatMultTransposeColumn(Y,B,2));
    if (verbose) PetscCall(BVView(Y,view));
  }

  /* Test BVGetMat/RestoreMat */
  PetscCall(BVGetMat(Y,&Ymat));
  PetscCall(PetscObjectSetName((PetscObject)Ymat,"Ymat"));
  if (verbose) PetscCall(MatView(Ymat,view));
  PetscCall(BVRestoreMat(Y,&Ymat));

  if (!fromfile) {
    /* Create BV object Z */
    PetscCall(BVDuplicateResize(Y,k,&Z));
    PetscCall(PetscObjectSetName((PetscObject)Z,"Z"));

    /* Fill Z entries */
    for (j=0;j<k;j++) {
      PetscCall(BVGetColumn(Z,j,&v));
      PetscCall(VecSet(v,0.0));
      if (!Istart) PetscCall(VecSetValue(v,0,1.0,ADD_VALUES));
      if (j<n && j>=Istart && j<Iend) PetscCall(VecSetValue(v,j,1.0,ADD_VALUES));
      if (j+1<n && j>=Istart && j<Iend) PetscCall(VecSetValue(v,j+1,-1.0,ADD_VALUES));
      PetscCall(VecAssemblyBegin(v));
      PetscCall(VecAssemblyEnd(v));
      PetscCall(BVRestoreColumn(Z,j,&v));
    }
    if (verbose) PetscCall(BVView(Z,view));

    /* Save a copy of Z */
    PetscCall(BVDuplicate(Z,&Zcopy));
    PetscCall(BVCopy(Z,Zcopy));

    /* Test BVMult, check result of previous operations */
    PetscCall(BVMult(Z,-1.0,1.0,Y,NULL));
    PetscCall(BVNorm(Z,NORM_FROBENIUS,&norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of error: %g\n",(double)norm));
  }

  /* Test BVMatMultColumn, multiply Y(:,2), result in Y(:,3) */
  if (m==n) {
    PetscCall(BVMatMultColumn(Y,B,2));
    if (verbose) PetscCall(BVView(Y,view));

    if (!fromfile) {
      /* Test BVGetArray, modify Z to match Y */
      PetscCall(BVCopy(Zcopy,Z));
      PetscCall(BVGetArray(Z,&pZ));
      if (Istart==0) {
        PetscCheck(Iend>2,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"First process must have at least 3 rows");
        PetscCall(BVGetLeadingDimension(Z,&ld));
        pZ[ld]   = 5.0;   /* modify 3 first entries of second column */
        pZ[ld+1] = -4.0;
        pZ[ld+2] = 1.0;
      }
      PetscCall(BVRestoreArray(Z,&pZ));
      if (verbose) PetscCall(BVView(Z,view));

      /* Check result again with BVMult */
      PetscCall(BVMult(Z,-1.0,1.0,Y,NULL));
      PetscCall(BVNorm(Z,NORM_FROBENIUS,&norm));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm of error: %g\n",(double)norm));
    }
  }

  PetscCall(BVDestroy(&Z));
  PetscCall(BVDestroy(&Zcopy));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&t));
  PetscCall(VecDestroy(&r));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test7_1.out
      filter: grep -v "Using method"
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_matmult vecs
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -mat_type aijcusparse -bv_matmult vecs
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -mat_type aijhipsparse -bv_matmult vecs
         requires: hip
      test:
         suffix: 1_mat
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_matmult mat

   testset:
      output_file: output/test7_2.out
      filter: grep -v "Using method" | sed -e "s/error: -0\./error: 0./"
      args: -m 34 -n 38 -k 9
      nsize: 2
      test:
         suffix: 2
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_matmult vecs
      test:
         suffix: 2_cuda
         args: -bv_type {{svec mat}} -mat_type aijcusparse -bv_matmult vecs
         requires: cuda
      test:
         suffix: 2_hip
         args: -bv_type {{svec mat}} -mat_type aijhipsparse -bv_matmult vecs
         requires: hip
      test:
         suffix: 2_mat
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_matmult mat

   testset:
      output_file: output/test7_3.out
      filter: grep -v "Using method"
      args: -file ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -bv_reproducible_random
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      nsize: 2
      test:
         suffix: 3
         args: -bv_type {{vecs contiguous svec mat}shared output} -bv_matmult {{vecs mat}}

TEST*/
