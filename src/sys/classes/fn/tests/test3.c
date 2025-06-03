/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix exponential.\n\n";

#include <slepcfn.h>

/*
   Compute matrix exponential B = expm(A)
 */
PetscErrorCode TestMatExp(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace,PetscBool checkerror)
{
  PetscScalar    tau,eta;
  PetscBool      set,flg;
  PetscInt       n;
  Mat            F,R,Finv,Acopy;
  Vec            v,f0;
  FN             finv;
  PetscReal      nrm,nrmf;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&F));
  PetscCall(PetscObjectSetName((PetscObject)F,"F"));
  /* compute matrix exponential */
  if (inplace) {
    PetscCall(MatCopy(A,F,SAME_NONZERO_PATTERN));
    PetscCall(MatIsHermitianKnown(A,&set,&flg));
    if (set && flg) PetscCall(MatSetOption(F,MAT_HERMITIAN,PETSC_TRUE));
    PetscCall(FNEvaluateFunctionMat(fn,F,NULL));
  } else {
    PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&Acopy));
    PetscCall(FNEvaluateFunctionMat(fn,A,F));
    /* check that A has not been modified */
    PetscCall(MatAXPY(Acopy,-1.0,A,SAME_NONZERO_PATTERN));
    PetscCall(MatNorm(Acopy,NORM_1,&nrm));
    if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the input matrix has changed by %g\n",(double)nrm));
    PetscCall(MatDestroy(&Acopy));
  }
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    PetscCall(MatView(A,viewer));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed expm(A) - - - - - - -\n"));
    PetscCall(MatView(F,viewer));
  }
  /* print matrix norm for checking */
  PetscCall(MatNorm(F,NORM_1,&nrmf));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The 1-norm of f(A) is %g\n",(double)nrmf));
  if (checkerror) {
    PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Finv));
    PetscCall(PetscObjectSetName((PetscObject)Finv,"Finv"));
    PetscCall(FNGetScale(fn,&tau,&eta));
    /* compute inverse exp(-tau*A)/eta */
    PetscCall(FNCreate(PETSC_COMM_WORLD,&finv));
    PetscCall(FNSetType(finv,FNEXP));
    PetscCall(FNSetFromOptions(finv));
    PetscCall(FNSetScale(finv,-tau,1.0/eta));
    if (inplace) {
      PetscCall(MatCopy(A,Finv,SAME_NONZERO_PATTERN));
      PetscCall(MatIsHermitianKnown(A,&set,&flg));
      if (set && flg) PetscCall(MatSetOption(Finv,MAT_HERMITIAN,PETSC_TRUE));
      PetscCall(FNEvaluateFunctionMat(finv,Finv,NULL));
    } else PetscCall(FNEvaluateFunctionMat(finv,A,Finv));
    PetscCall(FNDestroy(&finv));
    /* check error ||F*Finv-I||_F */
    PetscCall(MatMatMult(F,Finv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
    PetscCall(MatShift(R,-1.0));
    PetscCall(MatNorm(R,NORM_FROBENIUS,&nrm));
    if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||exp(A)*exp(-A)-I||_F < 100*eps\n"));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||exp(A)*exp(-A)-I||_F = %g\n",(double)nrm));
    PetscCall(MatDestroy(&R));
    PetscCall(MatDestroy(&Finv));
  }
  /* check FNEvaluateFunctionMatVec() */
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(F,f0,0));
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm/nrmf>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&F));
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
  PetscBool      verbose,inplace,checkerror,matcuda;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-checkerror",&checkerror));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix exponential, n=%" PetscInt_FMT ".\n",n));

  /* Create exponential function object */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNEXP));
  PetscCall(FNSetFromOptions(fn));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNView(fn,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Create matrices */
  if (matcuda) {
#if defined(PETSC_HAVE_CUDA)
    PetscCall(MatCreateSeqDenseCUDA(PETSC_COMM_SELF,n,n,NULL,&A));
#endif
  } else PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  /* Fill A with a symmetric Toeplitz matrix */
  PetscCall(MatDenseGetArray(A,&As));
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
      As[i+(i+j)*n]=1.0;
      As[(i+j)+i*n]=1.0;
    }
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(TestMatExp(fn,A,viewer,verbose,inplace,checkerror));

  /* Repeat with non-symmetric A */
  PetscCall(MatDenseGetArray(A,&As));
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) As[(i+j)+i*n]=-1.0;
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(TestMatExp(fn,A,viewer,verbose,inplace,checkerror));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: grep -v "computing matrix functions"
      output_file: output/test3_1.out
      test:
         suffix: 1
         args: -fn_method {{0 1}}
      test:
         suffix: 1_subdiagonalpade
         args: -fn_method {{2 3}}
         requires: c99_complex !single
      test:
         suffix: 1_cuda
         args: -fn_method 1 -matcuda
         requires: cuda !magma
      test:
         suffix: 1_magma
         args: -fn_method {{0 1 2 3}} -matcuda
         requires: cuda magma
      test:
         suffix: 2
         args: -inplace -fn_method{{0 1}}
      test:
         suffix: 2_subdiagonalpade
         args: -inplace -fn_method{{2 3}}
         requires: c99_complex !single
      test:
         suffix: 2_cuda
         args: -inplace -fn_method 1 -matcuda
         requires: cuda !magma
      test:
         suffix: 2_magma
         args: -inplace -fn_method {{0 1 2 3}} -matcuda
         requires: cuda magma

   testset:
      args: -fn_scale 0.1
      filter: grep -v "computing matrix functions"
      output_file: output/test3_3.out
      test:
         suffix: 3
         args: -fn_method {{0 1}}
      test:
        suffix: 3_subdiagonalpade
        args: -fn_method {{2 3}}
        requires: c99_complex !single

   testset:
      args: -n 120 -fn_scale 0.6,1.5
      filter: grep -v "computing matrix functions"
      output_file: output/test3_4.out
      test:
         suffix: 4
         args: -fn_method {{0 1}}
         requires: !single
      test:
        suffix: 4_subdiagonalpade
        args: -fn_method {{2 3}}
        requires: c99_complex !single

   test:
      suffix: 5
      args: -fn_scale 30 -fn_method {{2 3}}
      filter: grep -v "computing matrix functions"
      requires: c99_complex !single
      output_file: output/test3_5.out

   test:
      suffix: 6
      args: -fn_scale 1e-9 -fn_method {{2 3}}
      filter: grep -v "computing matrix functions"
      requires: c99_complex !single
      output_file: output/test3_6.out

TEST*/
