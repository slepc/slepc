/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix square root.\n\n";

#include <slepcfn.h>

/*
   Compute matrix square root B = sqrtm(A)
   Check result as norm(B*B-A)
 */
PetscErrorCode TestMatSqrt(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscScalar    tau,eta;
  PetscReal      nrm;
  PetscBool      set,flg;
  PetscInt       n;
  Mat            S,R,Acopy;
  Vec            v,f0;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&S));
  PetscCall(PetscObjectSetName((PetscObject)S,"S"));
  PetscCall(FNGetScale(fn,&tau,&eta));
  /* compute square root */
  if (inplace) {
    PetscCall(MatCopy(A,S,SAME_NONZERO_PATTERN));
    PetscCall(MatIsHermitianKnown(A,&set,&flg));
    if (set && flg) PetscCall(MatSetOption(S,MAT_HERMITIAN,PETSC_TRUE));
    PetscCall(FNEvaluateFunctionMat(fn,S,NULL));
  } else {
    PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&Acopy));
    PetscCall(FNEvaluateFunctionMat(fn,A,S));
    /* check that A has not been modified */
    PetscCall(MatAXPY(Acopy,-1.0,A,SAME_NONZERO_PATTERN));
    PetscCall(MatNorm(Acopy,NORM_1,&nrm));
    if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the input matrix has changed by %g\n",(double)nrm));
    PetscCall(MatDestroy(&Acopy));
  }
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    PetscCall(MatView(A,viewer));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed sqrtm(A) - - - - - - -\n"));
    PetscCall(MatView(S,viewer));
  }
  /* check error ||S*S-A||_F */
  PetscCall(MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
  if (eta!=1.0) PetscCall(MatScale(R,1.0/(eta*eta)));
  PetscCall(MatAXPY(R,-tau,A,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(R,NORM_FROBENIUS,&nrm));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S-A||_F < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S-A||_F = %g\n",(double)nrm));
  /* check FNEvaluateFunctionMatVec() */
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(S,f0,0));
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&S));
  PetscCall(MatDestroy(&R));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  FN             fn;
  Mat            A=NULL;
  PetscInt       i,j,n=10;
  PetscScalar    *As;
  PetscViewer    viewer;
  PetscBool      verbose,inplace,matcuda;
  PetscRandom    myrand;
  PetscReal      v;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix square root, n=%" PetscInt_FMT ".\n",n));

  /* Create function object */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNSQRT));
  PetscCall(FNSetFromOptions(fn));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNView(fn,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Create matrix */
  if (matcuda) {
#if PetscDefined(HAVE_CUDA)
    PetscCall(MatCreateSeqDenseCUDA(PETSC_COMM_SELF,n,n,NULL,&A));
#endif
  } else PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  /* Compute square root of a symmetric matrix A */
  PetscCall(MatDenseGetArray(A,&As));
  for (i=0;i<n;i++) As[i+i*n]=2.5;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) { As[i+(i+j)*n]=1.0; As[(i+j)+i*n]=1.0; }
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(TestMatSqrt(fn,A,viewer,verbose,inplace));

  /* Repeat with upper triangular A */
  PetscCall(MatDenseGetArray(A,&As));
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) As[(i+j)+i*n]=0.0;
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(TestMatSqrt(fn,A,viewer,verbose,inplace));

  /* Repeat with non-symmetic A */
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&myrand));
  PetscCall(PetscRandomSetFromOptions(myrand));
  PetscCall(PetscRandomSetInterval(myrand,0.0,1.0));
  PetscCall(MatDenseGetArray(A,&As));
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
      PetscCall(PetscRandomGetValueReal(myrand,&v));
      As[(i+j)+i*n]=v;
    }
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(PetscRandomDestroy(&myrand));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(TestMatSqrt(fn,A,viewer,verbose,inplace));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -fn_scale .05,2 -n 100
      filter: grep -v "computing matrix functions"
      output_file: output/test7_1.out
      requires: !__float128
      timeoutfactor: 2
      test:
         suffix: 1
         args: -fn_method {{0 1 2}}
      test:
         suffix: 1_sadeghi
         args: -fn_method 3
         requires: !single
      test:
         suffix: 1_cuda
         args: -fn_method 2 -matcuda
         requires: cuda !single
      test:
         suffix: 1_magma
         args: -fn_method {{1 3}} -matcuda
         requires: cuda magma !single
      test:
         suffix: 2
         args: -inplace -fn_method {{0 1 2}}
      test:
         suffix: 2_sadeghi
         args: -inplace -fn_method 3
         requires: !single
      test:
         suffix: 2_cuda
         args: -inplace -fn_method 2 -matcuda
         requires: cuda !single
      test:
         suffix: 2_magma
         args: -inplace -fn_method {{1 3}} -matcuda
         requires: cuda magma !single

   testset:
      nsize: 3
      args: -fn_scale .05,2 -n 100 -fn_parallel synchronized
      filter: grep -v "computing matrix functions" | grep -v "SYNCHRONIZED" | sed -e "s/3 MPI processes/1 MPI process/g"
      requires: !__float128
      output_file: output/test7_1.out
      test:
         suffix: 3
      test:
         suffix: 3_inplace
         args: -inplace

TEST*/
