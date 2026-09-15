/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test tensor BV.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec               t,v;
  Mat               S,M,Q;
  BV                U,V,UU;
  PetscInt          i,ii,j,jj,n=10,k=6,l=3,d=3,deg,id,lds;
  PetscScalar       *pS,*q;
  PetscReal         norm;
  PetscViewer       view;
  PetscBool         verbose;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-d",&d,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test tensor BV of degree %" PetscInt_FMT " with %" PetscInt_FMT " columns of dimension %" PetscInt_FMT "*d.\n",d,k,n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

  /* Create BV object U */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&U));
  PetscCall(PetscObjectSetName((PetscObject)U,"U"));
  PetscCall(BVSetSizesFromVec(U,t,k+d-1));
  PetscCall(BVSetFromOptions(U));
  PetscCall(PetscObjectSetName((PetscObject)U,"U"));

  /* Fill first d columns of U */
  for (j=0;j<d;j++) {
    PetscCall(BVGetColumn(U,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=0;i<4;i++) {
      if (i+j<n) PetscCall(VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(U,j,&v));
  }

  /* Create tensor BV */
  PetscCall(BVCreateTensor(U,d,&V));
  PetscCall(PetscObjectSetName((PetscObject)V,"V"));
  PetscCall(BVTensorGetDegree(V,&deg));
  PetscCheck(deg==d,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Wrong degree");

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(BVView(V,view));
  PetscCall(PetscViewerPopFormat(view));
  if (verbose) {
    PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(BVView(V,view));
  }

  /* Build first column from previously introduced coefficients */
  PetscCall(BVTensorBuildFirstColumn(V,d));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After building the first column - - - - -\n"));
    PetscCall(BVView(V,view));
  }

  /* Test orthogonalization */
  PetscCall(BVTensorGetFactors(V,&UU,&S));
  PetscCall(BVGetActiveColumns(UU,NULL,&j));
  PetscCall(BVGetSizes(UU,NULL,NULL,&id));
  PetscCheck(id==k+d-1,PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Wrong dimensions");
  lds = id*d;
  for (jj=1;jj<k;jj++) {
    /* set new orthogonal column in U */
    PetscCall(BVGetColumn(UU,j,&v));
    PetscCall(VecSet(v,0.0));
    for (i=0;i<4;i++) {
      if (i+j<n) PetscCall(VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(v));
    PetscCall(VecAssemblyEnd(v));
    PetscCall(BVRestoreColumn(UU,j,&v));
    PetscCall(BVOrthonormalizeColumn(UU,j,PETSC_TRUE,NULL,NULL));
    j++;
    PetscCall(BVSetActiveColumns(UU,0,j));
    /* set new column of S */
    PetscCall(MatDenseGetArray(S,&pS));
    for (ii=0;ii<d;ii++) {
      for (i=0;i<ii+jj+1;i++) {
        pS[i+ii*id+jj*lds] = (PetscScalar)(2*ii+i+0.5*jj);
      }
    }
    PetscCall(MatDenseRestoreArray(S,&pS));
    PetscCall(BVOrthonormalizeColumn(V,jj,PETSC_TRUE,NULL,NULL));
  }
  PetscCall(BVTensorRestoreFactors(V,&UU,&S));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After orthogonalization - - - - -\n"));
    PetscCall(BVView(V,view));
  }

  /* Check orthogonality */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));
  PetscCall(BVDot(V,V,M));
  PetscCall(MatShift(M,-1.0));
  PetscCall(MatNorm(M,NORM_1,&norm));
  if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)norm));

  /* Test BVTensorCompress */
  PetscCall(BVSetActiveColumns(V,0,l));
  PetscCall(BVTensorCompress(V,0));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVTensorCompress - - - - -\n"));
    PetscCall(BVView(V,view));
  }

  /* Create Mat */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,l,NULL,&Q));
  PetscCall(PetscObjectSetName((PetscObject)Q,"Q"));
  PetscCall(MatDenseGetArray(Q,&q));
  for (i=0;i<k;i++)
    for (j=0;j<l;j++)
      q[i+j*k] = (i<j)? 2.0: -0.5;
  PetscCall(MatDenseRestoreArray(Q,&q));
  if (verbose) PetscCall(MatView(Q,NULL));

  /* Test BVMultInPlace */
  PetscCall(BVMultInPlace(V,Q,1,l));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultInPlace - - - - -\n"));
    PetscCall(BVView(V,view));
  }

  /* Test BVNorm */
  PetscCall(BVNorm(V,NORM_1,&norm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Norm: %g\n",(double)norm));

  PetscCall(BVDestroy(&U));
  PetscCall(BVDestroy(&V));
  PetscCall(MatDestroy(&Q));
  PetscCall(MatDestroy(&M));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      nsize: 2
      output_file: output/test16_1.out
      filter: grep -v "doing matmult"
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}}
      test:
         suffix: 1_cuda
         args: -bv_type {{vecs svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{vecs svec mat}} -vec_type hip
         requires: hip
      test:
         suffix: 1_kokkos
         args: -bv_type mat -vec_type kokkos
         requires: hip kokkos_kernels

TEST*/
