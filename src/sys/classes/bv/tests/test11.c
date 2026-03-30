/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV block orthogonalization.\n\n";

#include <slepcbv.h>

/*
   Compute the Frobenius norm ||A(l:k,l:k)-diag||_F
 */
PetscErrorCode MyMatNorm(Mat A,PetscInt lda,PetscInt l,PetscInt k,PetscScalar diag,PetscReal *norm)
{
  PetscInt          i,j;
  const PetscScalar *pA;
  PetscReal         s,val;

  PetscFunctionBeginUser;
  PetscCall(MatDenseGetArrayRead(A,&pA));
  s = 0.0;
  for (i=l;i<k;i++) {
    for (j=l;j<k;j++) {
      val = (i==j)? PetscAbsScalar(pA[i+j*lda]-diag): PetscAbsScalar(pA[i+j*lda]);
      s += val*val;
    }
  }
  *norm = PetscSqrtReal(s);
  PetscCall(MatDenseRestoreArrayRead(A,&pA));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  BV             X,Y,Z,cached;
  Mat            B=NULL,M,R=NULL;
  Vec            v,t;
  PetscInt       i,j,n=20,l=2,k=8,Istart,Iend;
  PetscViewer    view;
  PetscBool      withb,resid,rand,verbose;
  PetscReal      norm;
  PetscScalar    alpha;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-withb",&withb));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-resid",&resid));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-rand",&rand));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV block orthogonalization (length %" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT ")%s.\n",n,l,k,withb?" with non-standard inner product":""));

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

  /* Fill X entries */
  if (rand) PetscCall(BVSetRandom(X));
  else {
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
  }
  if (verbose) PetscCall(BVView(X,view));

  if (withb) {
    /* Create inner product matrix */
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
    if (verbose) PetscCall(MatView(B,view));
    PetscCall(BVSetMatrix(X,B,PETSC_FALSE));
  }

  /* Create copy on Y */
  PetscCall(BVDuplicate(X,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  PetscCall(BVCopy(X,Y));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));

  if (resid) {
    /* Create matrix R to store triangular factor */
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&R));
    PetscCall(PetscObjectSetName((PetscObject)R,"R"));
  }

  if (l>0) {
    /* First orthogonalize leading columns */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Orthogonalizing leading columns\n"));
    PetscCall(BVSetActiveColumns(Y,0,l));
    PetscCall(BVSetActiveColumns(X,0,l));
    PetscCall(BVOrthogonalize(Y,R));
    if (verbose) {
      PetscCall(BVView(Y,view));
      if (resid) PetscCall(MatView(R,view));
    }

    if (withb) {
      /* Extract cached BV and check it is equal to B*X */
      PetscCall(BVGetCachedBV(Y,&cached));
      PetscCall(BVDuplicate(X,&Z));
      PetscCall(BVSetMatrix(Z,NULL,PETSC_FALSE));
      PetscCall(BVSetActiveColumns(Z,0,l));
      PetscCall(BVCopy(X,Z));
      PetscCall(BVMatMult(X,B,Z));
      PetscCall(BVMult(Z,-1.0,1.0,cached,NULL));
      PetscCall(BVNorm(Z,NORM_FROBENIUS,&norm));
      if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Difference ||cached-BX|| < 100*eps\n"));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Difference ||cached-BX||: %g\n",(double)norm));
      PetscCall(BVDestroy(&Z));
    }

    /* Check orthogonality */
    PetscCall(BVDot(Y,Y,M));
    PetscCall(MyMatNorm(M,k,0,l,1.0,&norm));
    if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q1 < 100*eps\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q1: %g\n",(double)norm));

    if (resid) {
      /* Check residual */
      PetscCall(BVDuplicate(X,&Z));
      PetscCall(BVSetMatrix(Z,NULL,PETSC_FALSE));
      PetscCall(BVSetActiveColumns(Z,0,l));
      PetscCall(BVCopy(X,Z));
      PetscCall(BVMult(Z,-1.0,1.0,Y,R));
      PetscCall(BVNorm(Z,NORM_FROBENIUS,&norm));
      if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X1-Q1*R11|| < 100*eps\n"));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X1-Q1*R11||: %g\n",(double)norm));
      PetscCall(BVDestroy(&Z));
    }

  }

  /* Now orthogonalize the rest of columns */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Orthogonalizing active columns\n"));
  PetscCall(BVSetActiveColumns(Y,l,k));
  PetscCall(BVSetActiveColumns(X,l,k));
  PetscCall(BVOrthogonalize(Y,R));
  if (verbose) {
    PetscCall(BVView(Y,view));
    if (resid) PetscCall(MatView(R,view));
  }

  if (l>0) {
    /* Check orthogonality */
    PetscCall(BVDot(Y,Y,M));
    PetscCall(MyMatNorm(M,k,l,k,1.0,&norm));
    if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q2 < 100*eps\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q2: %g\n",(double)norm));
  }

  /* Check the complete decomposition */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Overall decomposition\n"));
  PetscCall(BVSetActiveColumns(Y,0,k));
  PetscCall(BVSetActiveColumns(X,0,k));

  /* Check orthogonality */
  PetscCall(BVDot(Y,Y,M));
  PetscCall(MyMatNorm(M,k,0,k,1.0,&norm));
  if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q: %g\n",(double)norm));

  if (resid) {
    /* Check residual */
    PetscCall(BVMult(X,-1.0,1.0,Y,R));
    PetscCall(BVSetMatrix(X,NULL,PETSC_FALSE));
    PetscCall(BVNorm(X,NORM_FROBENIUS,&norm));
    if (norm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X-Q*R|| < 100*eps\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X-Q*R||: %g\n",(double)norm));
    PetscCall(MatDestroy(&R));
  }

  PetscCall(MatDestroy(&B));
  PetscCall(MatDestroy(&M));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -bv_orthog_block {{gs chol tsqr tsqrchol svqb}}
      nsize: 2
      output_file: output/test11_1.out
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

   testset:
      args: -withb -bv_orthog_block {{gs chol svqb}}
      nsize: 2
      output_file: output/test11_4.out
      test:
         suffix: 4
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 4_cuda
         args: -bv_type {{svec mat}} -vec_type cuda -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 4_hip
         args: -bv_type {{svec mat}} -vec_type hip -mat_type aijhipsparse
         requires: hip

   testset:
      args: -resid -bv_orthog_block {{gs chol tsqr tsqrchol svqb}}
      nsize: 2
      output_file: output/test11_6.out
      test:
         suffix: 6
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 6_cuda
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda
      test:
         suffix: 6_hip
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip

   testset:
      args: -resid -withb -bv_orthog_block {{gs chol svqb}}
      nsize: 2
      output_file: output/test11_9.out
      test:
         suffix: 9
         args: -bv_type {{vecs contiguous svec mat}shared output}
      test:
         suffix: 9_cuda
         args: -bv_type {{svec mat}} -vec_type cuda -mat_type aijcusparse
         requires: cuda
      test:
         suffix: 9_hip
         args: -bv_type {{svec mat}} -vec_type hip -mat_type aijhipsparse
         requires: hip

   testset:
      args: -bv_orthog_block tsqr
      nsize: 7
      output_file: output/test11_1.out
      test:
         suffix: 11
         args: -bv_type {{vecs contiguous svec mat}shared output}
         requires: !defined(PETSCTEST_VALGRIND)
      test:
         suffix: 11_cuda
         TODO: too many processes accessing the GPU
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda !defined(PETSCTEST_VALGRIND)
      test:
         suffix: 11_hip
         TODO: too many processes accessing the GPU
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip !defined(PETSCTEST_VALGRIND)

   testset:
      args: -resid -n 180 -l 0 -k 7 -bv_orthog_block tsqr
      nsize: 7
      output_file: output/test11_12.out
      test:
         suffix: 12
         args: -bv_type {{vecs contiguous svec mat}shared output}
         requires: double !defined(PETSCTEST_VALGRIND)
      test:
         suffix: 12_cuda
         TODO: too many processes accessing the GPU
         args: -bv_type {{svec mat}} -vec_type cuda
         requires: cuda !single !defined(PETSCTEST_VALGRIND)
      test:
         suffix: 12_hip
         TODO: too many processes accessing the GPU
         args: -bv_type {{svec mat}} -vec_type hip
         requires: hip !single !defined(PETSCTEST_VALGRIND)

TEST*/
