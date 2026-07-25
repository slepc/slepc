/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix inverse square root.\n\n";

#include <slepcfn.h>

/*
   Compute matrix inverse square root B = inv(sqrtm(A))
   Check result as norm(B*B*A-I)
 */
PetscErrorCode TestMatInvSqrt(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
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
  /* compute inverse square root */
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
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed inv(sqrtm(A)) - - - - - - -\n"));
    PetscCall(MatView(S,viewer));
  }
  /* check error ||S*S*A-I||_F */
  PetscCall(MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
  if (eta!=1.0) PetscCall(MatScale(R,1.0/(eta*eta)));
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(S,f0,0));
  PetscCall(MatCopy(R,S,SAME_NONZERO_PATTERN));
  PetscCall(MatDestroy(&R));
  if (tau!=1.0) PetscCall(MatScale(S,tau));
  PetscCall(MatMatMult(S,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R));
  PetscCall(MatShift(R,-1.0));
  PetscCall(MatNorm(R,NORM_FROBENIUS,&nrm));
  PetscCall(MatDestroy(&R));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||S*S*A-I||_F = %g\n",(double)nrm));
  /* check FNEvaluateFunctionMatVec() */
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&S));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  FN             fn;
  Mat            A=NULL;
  PetscInt       i,j,n=10;
  PetscScalar    x,y,yp,*As;
  PetscViewer    viewer;
  PetscBool      verbose,inplace,matcuda;
  PetscRandom    myrand;
  PetscReal      v;
  char           strx[50],str[50];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix inverse square root, n=%" PetscInt_FMT ".\n",n));

  /* Create function object */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNINVSQRT));
  PetscCall(FNSetFromOptions(fn));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNView(fn,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Scalar evaluation */
  x = 2.2;
  PetscCall(SlepcSNPrintfScalar(strx,sizeof(strx),x,PETSC_FALSE));
  PetscCall(FNEvaluateFunction(fn,x,&y));
  PetscCall(FNEvaluateDerivative(fn,x,&yp));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),y,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),yp,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str));

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
  PetscCall(TestMatInvSqrt(fn,A,viewer,verbose,inplace));

  /* Repeat with upper triangular A */
  PetscCall(MatDenseGetArray(A,&As));
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) As[(i+j)+i*n]=0.0;
  }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(TestMatInvSqrt(fn,A,viewer,verbose,inplace));

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
  PetscCall(TestMatInvSqrt(fn,A,viewer,verbose,inplace));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -fn_scale 0.9,0.5 -n 10
      filter: grep -v "computing matrix functions"
      requires: !__float128
      output_file: output/test8_1.out
      test:
         suffix: 1
         args: -fn_method {{0 1 2 3}}
      test:
         suffix: 1_cuda
         args: -fn_method 2 -matcuda
         requires: cuda
      test:
         suffix: 1_magma
         args: -fn_method {{1 3}} -matcuda
         requires: cuda magma
      test:
         suffix: 2
         args: -inplace -fn_method {{0 1 2 3}}
      test:
         suffix: 2_cuda
         args: -inplace -fn_method 2 -matcuda
         requires: cuda
      test:
         suffix: 2_magma
         args: -inplace -fn_method {{1 3}} -matcuda
         requires: cuda magma

TEST*/
