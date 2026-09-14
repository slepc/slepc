/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec               t,v;
  Mat               Q=NULL,M=NULL;
  BV                X,Y;
  PetscInt          i,j,n=10,k=5,l=3,ldx,lda;
  PetscMPIInt       rank;
  PetscScalar       *q,*z;
  const PetscScalar *pX;
  PetscReal         nrm;
  PetscViewer       view;
  PetscBool         verbose,matcuda,mathip,testlda=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-mathip",&mathip));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-testlda",&testlda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV with %" PetscInt_FMT " columns of dimension %" PetscInt_FMT ".\n",k,n));

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
  PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_INFO_DETAIL));
  PetscCall(BVView(X,view));
  PetscCall(PetscViewerPopFormat(view));

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

  /* Create BV object Y */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVSetSizesFromVec(Y,t,l));
  PetscCall(BVSetFromOptions(Y));

  /* Fill Y entries */
  for (j=0;j<l;j++) {
    PetscCall(BVGetColumn(Y,j,&v));
    PetscCall(VecSet(v,(PetscScalar)(j+1)/4.0));
    PetscCall(BVRestoreColumn(Y,j,&v));
  }
  if (verbose) PetscCall(BVView(Y,view));

  /* Create Mat */
  PetscCall(MatCreate(PETSC_COMM_SELF,&Q));
  if (matcuda && PetscDefined(HAVE_CUDA)) PetscCall(MatSetType(Q,MATSEQDENSECUDA));
  else if (mathip && PetscDefined(HAVE_HIP)) PetscCall(MatSetType(Q,MATSEQDENSEHIP));
  else PetscCall(MatSetType(Q,MATSEQDENSE));
  PetscCall(MatSetSizes(Q,k,l,k,l));
  if (testlda) PetscCall(MatDenseSetLDA(Q,k+2));
  PetscCall(MatSeqDenseSetPreallocation(Q,NULL));
  PetscCall(PetscObjectSetName((PetscObject)Q,"Q"));
  PetscCall(MatDenseGetArrayWrite(Q,&q));
  PetscCall(MatDenseGetLDA(Q,&lda));
  for (i=0;i<k;i++)
    for (j=0;j<l;j++)
      q[i+j*lda] = (i<j)? 2.0: -0.5;
  PetscCall(MatDenseRestoreArrayWrite(Q,&q));
  if (verbose) PetscCall(MatView(Q,NULL));

  /* Test BVMult */
  PetscCall(BVMult(Y,2.0,1.0,X,Q));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMult - - - - - - - - -\n"));
    PetscCall(BVView(Y,view));
  }

  /* Test BVMultVec */
  PetscCall(BVGetColumn(Y,0,&v));
  PetscCall(PetscMalloc1(k,&z));
  z[0] = 2.0;
  for (i=1;i<k;i++) z[i] = -0.5*z[i-1];
  PetscCall(BVMultVec(X,-1.0,1.0,v,z));
  PetscCall(PetscFree(z));
  PetscCall(BVRestoreColumn(Y,0,&v));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultVec - - - - - - -\n"));
    PetscCall(BVView(Y,view));
  }

  /* Test BVDot */
  PetscCall(MatCreate(PETSC_COMM_SELF,&M));
  if (matcuda && PetscDefined(HAVE_CUDA)) PetscCall(MatSetType(M,MATSEQDENSECUDA));
  else if (mathip && PetscDefined(HAVE_HIP)) PetscCall(MatSetType(M,MATSEQDENSEHIP));
  else PetscCall(MatSetType(M,MATSEQDENSE));
  PetscCall(MatSetSizes(M,l,k,l,k));
  if (testlda) PetscCall(MatDenseSetLDA(M,l+2));
  PetscCall(MatSeqDenseSetPreallocation(M,NULL));
  PetscCall(PetscObjectSetName((PetscObject)M,"M"));
  PetscCall(BVDot(X,Y,M));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVDot - - - - - - - - -\n"));
    PetscCall(MatView(M,NULL));
  }

  /* Test BVDotVec */
  PetscCall(BVGetColumn(Y,0,&v));
  PetscCall(PetscMalloc1(k,&z));
  PetscCall(BVDotVec(X,v,z));
  PetscCall(BVRestoreColumn(Y,0,&v));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVDotVec - - - - - - -\n"));
    PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF,1,k,z,&v));
    PetscCall(PetscObjectSetName((PetscObject)v,"z"));
    PetscCall(VecView(v,view));
    PetscCall(VecDestroy(&v));
  }
  PetscCall(PetscFree(z));

  /* Test BVMultInPlace and BVScale */
  PetscCall(BVMultInPlace(X,Q,1,l));
  PetscCall(BVScale(X,2.0));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultInPlace - - - - -\n"));
    PetscCall(BVView(X,view));
  }

  /* Test BVNorm */
  PetscCall(BVNormColumn(X,0,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"2-Norm of X[0] = %g\n",(double)nrm));
  PetscCall(BVNorm(X,NORM_FROBENIUS,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Frobenius Norm of X = %g\n",(double)nrm));

  /* Test BVGetArrayRead */
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  if (!rank) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"First row of X =\n"));
    PetscCall(BVGetLeadingDimension(X,&ldx));
    PetscCall(BVGetArrayRead(X,&pX));
    for (i=0;i<k;i++) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%g ",(double)PetscRealPart(pX[i*ldx])));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(BVRestoreArrayRead(X,&pX));
  }

  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(MatDestroy(&Q));
  PetscCall(MatDestroy(&M));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   test:
      args: -bv_type {{vecs contiguous svec mat}separate output} -verbose
      suffix: 1
      filter: sed -e 's/-0[.]/0./g'

   testset:
      args: -bv_type svec -vec_type cuda -verbose
      requires: cuda
      output_file: output/test1_1_svec_gpu.out
      test:
         suffix: 1_svec_cuda
      test:
         suffix: 1_svec_cuda_mat
         args: -matcuda
         filter: sed -e "s/seqdensecuda/seqdense/"

   testset:
      args: -bv_type mat -vec_type cuda -verbose
      requires: cuda
      output_file: output/test1_1_mat_gpu.out
      filter: sed -e "s/seqdensecuda/seqdense/"
      test:
         suffix: 1_mat_cuda
      test:
         suffix: 1_mat_cuda_mat
         args: -matcuda

   testset:
      args: -bv_type svec -vec_type hip -verbose
      requires: hip
      output_file: output/test1_1_svec_gpu.out
      filter: sed -e "s/seqhip/seqcuda/" | sed -e "s/HIP/CUDA/"
      test:
         suffix: 1_svec_hip
      test:
         suffix: 1_svec_hip_mat
         args: -mathip
         filter: sed -e "s/seqhip/seqcuda/" | sed -e "s/seqdensehip/seqdense/" | sed -e "s/HIP/CUDA/"

   testset:
      args: -bv_type mat -vec_type hip -verbose
      requires: hip
      output_file: output/test1_1_mat_gpu.out
      filter: sed -e "s/seqdensehip/seqdense/" | sed -e "s/HIP/CUDA/"
      test:
         suffix: 1_mat_hip
      test:
         suffix: 1_mat_hip_mat
         args: -mathip

   test:
      args: -bv_type mat -vec_type kokkos -verbose
      suffix: 1_mat_kokkos
      requires: kokkos_kernels !cuda !hip
      output_file: output/test1_1_bv_type-mat.out
      filter: sed -e "s/-0[.]/0./g" | sed -e "s/ (Kokkos)//"

   test:
      args: -bv_type mat -vec_type kokkos
      suffix: 1_mat_kokkos_np2
      nsize: 2
      requires: kokkos_kernels !cuda !hip

   test:
      args: -bv_type {{vecs contiguous svec mat}separate output} -verbose -testlda
      suffix: 2
      filter: sed -e 's/-0[.]/0./g'

   testset:
      args: -bv_type svec -vec_type cuda -verbose -testlda
      requires: cuda
      output_file: output/test1_1_svec_gpu.out
      test:
         suffix: 2_svec_cuda
      test:
         suffix: 2_svec_cuda_mat
         args: -matcuda
         filter: sed -e "s/seqdensecuda/seqdense/"

   testset:
      args: -bv_type mat -vec_type cuda -verbose -testlda
      requires: cuda
      output_file: output/test1_1_mat_gpu.out
      filter: sed -e "s/seqdensecuda/seqdense/"
      test:
         suffix: 2_mat_cuda
      test:
         suffix: 2_mat_cuda_mat
         args: -matcuda

   testset:
      args: -bv_type svec -vec_type hip -verbose -testlda
      requires: hip
      output_file: output/test1_1_svec_gpu.out
      test:
         suffix: 2_svec_hip
         filter: sed -e "s/seqhip/seqcuda/" | sed -e "s/HIP/CUDA/"
      test:
         suffix: 2_svec_hip_mat
         args: -mathip
         filter: sed -e "s/seqhip/seqcuda/" | sed -e "s/seqdensehip/seqdense/" | sed -e "s/HIP/CUDA/"

   testset:
      args: -bv_type mat -vec_type hip -verbose -testlda
      requires: hip
      output_file: output/test1_1_mat_gpu.out
      filter: sed -e "s/seqdensehip/seqdense/" | sed -e "s/HIP/CUDA/"
      test:
         suffix: 2_mat_hip
      test:
         suffix: 2_mat_hip_mat
         args: -mathip

TEST*/
