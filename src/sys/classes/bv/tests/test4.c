/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations, changing the number of active columns.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  Vec            t,v;
  Mat            Q,M;
  BV             X,Y;
  PetscInt       i,j,n=10,kx=6,lx=3,ky=5,ly=2;
  PetscScalar    *q,*z;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose,trans;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-kx",&kx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-lx",&lx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ky",&ky,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ly",&ly,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"First BV with %" PetscInt_FMT " active columns (%" PetscInt_FMT " leading columns) of dimension %" PetscInt_FMT ".\n",kx,lx,n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Second BV with %" PetscInt_FMT " active columns (%" PetscInt_FMT " leading columns) of dimension %" PetscInt_FMT ".\n",ky,ly,n));

  /* Create template vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&t));
  PetscCall(VecSetSizes(t,PETSC_DECIDE,n));
  PetscCall(VecSetFromOptions(t));

  /* Create BV object X */
  PetscCall(BVCreate(PETSC_COMM_WORLD,&X));
  PetscCall(PetscObjectSetName((PetscObject)X,"X"));
  PetscCall(BVSetSizesFromVec(X,t,kx+2));  /* two extra columns to test active columns */
  PetscCall(BVSetFromOptions(X));
  PetscCall(BVSetActiveColumns(X,lx,kx));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view));
  if (verbose) PetscCall(PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB));

  /* Fill X entries */
  for (j=0;j<kx+2;j++) {
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
  PetscCall(BVSetSizesFromVec(Y,t,ky+1));
  PetscCall(BVSetFromOptions(Y));
  PetscCall(BVSetActiveColumns(Y,ly,ky));

  /* Fill Y entries */
  for (j=0;j<ky+1;j++) {
    PetscCall(BVGetColumn(Y,j,&v));
    PetscCall(VecSet(v,(PetscScalar)(j+1)/4.0));
    PetscCall(BVRestoreColumn(Y,j,&v));
  }
  if (verbose) PetscCall(BVView(Y,view));

  /* Create Mat */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,kx,ky,NULL,&Q));
  PetscCall(PetscObjectSetName((PetscObject)Q,"Q"));
  PetscCall(MatDenseGetArray(Q,&q));
  for (i=0;i<kx;i++)
    for (j=0;j<ky;j++)
      q[i+j*kx] = (i<j)? 2.0: -0.5;
  PetscCall(MatDenseRestoreArray(Q,&q));
  if (verbose) PetscCall(MatView(Q,NULL));

  /* Test BVResize */
  PetscCall(BVResize(X,kx+4,PETSC_TRUE));

  /* Test BVMult */
  PetscCall(BVMult(Y,2.0,0.5,X,Q));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMult - - - - - - - - -\n"));
    PetscCall(BVView(Y,view));
  }

  /* Test BVMultVec */
  PetscCall(BVGetColumn(Y,0,&v));
  PetscCall(PetscMalloc1(kx-lx,&z));
  z[0] = 2.0;
  for (i=1;i<kx-lx;i++) z[i] = -0.5*z[i-1];
  PetscCall(BVMultVec(X,-1.0,1.0,v,z));
  PetscCall(PetscFree(z));
  PetscCall(BVRestoreColumn(Y,0,&v));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultVec - - - - - - -\n"));
    PetscCall(BVView(Y,view));
  }

  /* Test BVDot */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&M));
  PetscCall(PetscObjectSetName((PetscObject)M,"M"));
  PetscCall(BVDot(X,Y,M));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVDot - - - - - - - - -\n"));
    PetscCall(MatView(M,NULL));
  }

  /* Test BVDotVec */
  PetscCall(BVGetColumn(Y,0,&v));
  PetscCall(PetscMalloc1(kx-lx,&z));
  PetscCall(BVDotVec(X,v,z));
  PetscCall(BVRestoreColumn(Y,0,&v));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVDotVec - - - - - - -\n"));
    PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF,1,kx-lx,z,&v));
    PetscCall(PetscObjectSetName((PetscObject)v,"z"));
    PetscCall(VecView(v,view));
    PetscCall(VecDestroy(&v));
  }
  PetscCall(PetscFree(z));

  /* Test BVMultInPlace and BVScale */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-trans",&trans));
  if (trans) {
    Mat Qt;
    PetscCall(MatTranspose(Q,MAT_INITIAL_MATRIX,&Qt));
    PetscCall(BVMultInPlaceHermitianTranspose(X,Qt,lx+1,ky));
    PetscCall(MatDestroy(&Qt));
  } else PetscCall(BVMultInPlace(X,Q,lx+1,ky));
  PetscCall(BVScale(X,2.0));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"After BVMultInPlace - - - - -\n"));
    PetscCall(BVView(X,view));
  }

  /* Test BVNorm */
  PetscCall(BVNormColumn(X,lx,NORM_2,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"2-Norm of X[%" PetscInt_FMT "] = %g\n",lx,(double)nrm));
  PetscCall(BVNorm(X,NORM_FROBENIUS,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Frobenius Norm of X = %g\n",(double)nrm));

  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(MatDestroy(&Q));
  PetscCall(MatDestroy(&M));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test4_1.out
      args: -n 18 -kx 12 -ky 8
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}}
      test:
         suffix: 1_vecs_vmip
         args: -bv_type vecs -bv_vecs_vmip 1
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
         args: -bv_type {{vecs contiguous svec mat}} -trans
      test:
         suffix: 2_cuda
         args: -bv_type {{svec mat}} -vec_type cuda -trans
         requires: cuda
      test:
         suffix: 2_hip
         args: -bv_type {{svec mat}} -vec_type hip -trans
         requires: hip

TEST*/
