/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV bi-orthogonalization functions.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  BV             X,Y;
  Mat            M;
  Vec            v,t;
  PetscInt       i,j,n=20,k=7;
  PetscViewer    view;
  PetscBool      verbose;
  PetscReal      norm,condn=1.0;
  PetscScalar    alpha;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-condn",&condn,NULL));
  PetscCheck(condn>=1.0,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The condition number must be > 1");
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Test BV bi-orthogonalization with %" PetscInt_FMT " columns of length %" PetscInt_FMT ".\n",k,n));
  if (condn>1.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," - Using random BVs with condition number = %g\n",(double)condn));

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
  if (condn==1.0) {
    for (j=0;j<k;j++) {
      PetscCall(BVGetColumn(X,j,&v));
      PetscCall(VecSet(v,0.0));
      for (i=0;i<=n/2;i++) {
        if (i+j<n) {
          alpha = (3.0*i+j-2)/(2*(i+j+1));
#if PetscDefined(USE_COMPLEX)
          alpha += (1.2*i+j-2)/(0.1*(i+j+1))*PETSC_i;
#endif
          PetscCall(VecSetValue(v,i+j,alpha,INSERT_VALUES));
        }
      }
      PetscCall(VecAssemblyBegin(v));
      PetscCall(VecAssemblyEnd(v));
      PetscCall(BVRestoreColumn(X,j,&v));
    }
  } else PetscCall(BVSetRandomCond(X,condn));
  if (verbose) PetscCall(BVView(X,view));

  /* Create Y and fill its entries */
  PetscCall(BVDuplicate(X,&Y));
  PetscCall(PetscObjectSetName((PetscObject)Y,"Y"));
  if (condn==1.0) {
    for (j=0;j<k;j++) {
      PetscCall(BVGetColumn(Y,j,&v));
      PetscCall(VecSet(v,0.0));
      for (i=PetscMin(n,2+(2*j)%6);i<PetscMin(n,6+(3*j)%9);i++) {
        if (i%5 && i!=j) {
          alpha = (1.5*i+j)/(2.2*(i-j));
          PetscCall(VecSetValue(v,i+j,alpha,INSERT_VALUES));
        }
      }
      PetscCall(VecAssemblyBegin(v));
      PetscCall(VecAssemblyEnd(v));
      PetscCall(BVRestoreColumn(Y,j,&v));
    }
  } else PetscCall(BVSetRandomCond(Y,condn));
  if (verbose) PetscCall(BVView(Y,view));

  /* Test BVBiorthonormalizeColumn */
  for (j=0;j<k;j++) PetscCall(BVBiorthonormalizeColumn(X,Y,j,NULL));
  if (verbose) {
    PetscCall(BVView(X,view));
    PetscCall(BVView(Y,view));
  }

  /* Check orthogonality */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M));
  PetscCall(PetscObjectSetName((PetscObject)M,"M"));
  PetscCall(BVDot(X,Y,M));
  if (verbose) PetscCall(MatView(M,view));
  PetscCall(MatShift(M,-1.0));
  PetscCall(MatNorm(M,NORM_1,&norm));
  if (norm<200*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of bi-orthogonality < 200*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Level of bi-orthogonality: %g\n",(double)norm));

  PetscCall(MatDestroy(&M));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&Y));
  PetscCall(VecDestroy(&t));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test17_1.out
      test:
         suffix: 1
         args: -bv_type {{vecs contiguous svec mat}} -bv_orthog_type cgs
         requires: double
      test:
         suffix: 1_cuda
         args: -bv_type {{svec mat}} -vec_type cuda -bv_orthog_type cgs
         requires: cuda
      test:
         suffix: 1_hip
         args: -bv_type {{svec mat}} -vec_type hip -bv_orthog_type cgs
         requires: hip
      test:
         suffix: 2
         args: -bv_type {{vecs contiguous svec mat}} -bv_orthog_type mgs
      test:
         suffix: 2_cuda
         args: -bv_type {{svec mat}} -vec_type cuda -bv_orthog_type mgs
         requires: cuda
      test:
         suffix: 2_hip
         args: -bv_type {{svec mat}} -vec_type hip -bv_orthog_type mgs
         requires: hip

TEST*/
