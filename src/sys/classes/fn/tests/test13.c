/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix logarithm.\n\n";

#include <slepcfn.h>

/*
   Compute matrix logarithm B = logm(A)
 */
PetscErrorCode TestMatLog(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscBool      set,flg;
  PetscScalar    tau,eta;
  PetscInt       n;
  Mat            F,R;
  Vec            v,f0;
  FN             fnexp;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&F));
  PetscCall(PetscObjectSetName((PetscObject)F,"F"));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&R));
  PetscCall(PetscObjectSetName((PetscObject)R,"R"));
  PetscCall(FNGetScale(fn,&tau,&eta));
  /* compute matrix logarithm */
  if (inplace) {
    PetscCall(MatCopy(A,F,SAME_NONZERO_PATTERN));
    PetscCall(MatIsHermitianKnown(A,&set,&flg));
    if (set && flg) PetscCall(MatSetOption(F,MAT_HERMITIAN,PETSC_TRUE));
    PetscCall(FNEvaluateFunctionMat(fn,F,NULL));
  } else PetscCall(FNEvaluateFunctionMat(fn,A,F));
  if (verbose) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n"));
    PetscCall(MatView(A,viewer));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed logm(A) - - - - - - -\n"));
    PetscCall(MatView(F,viewer));
  }
  /* check error ||expm(F)-A||_F */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fnexp));
  PetscCall(FNSetType(fnexp,FNEXP));
  PetscCall(MatCopy(F,R,SAME_NONZERO_PATTERN));
  if (eta!=1.0) PetscCall(MatScale(R,1.0/eta));
  PetscCall(FNEvaluateFunctionMat(fnexp,R,NULL));
  PetscCall(FNDestroy(&fnexp));
  PetscCall(MatAXPY(R,-tau,A,SAME_NONZERO_PATTERN));
  PetscCall(MatNorm(R,NORM_FROBENIUS,&nrm));
  if (nrm<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||expm(F)-A||_F < 100*eps\n"));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"||expm(F)-A||_F = %g\n",(double)nrm));
  /* check FNEvaluateFunctionMatVec() */
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(F,f0,0));
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&F));
  PetscCall(MatDestroy(&R));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  FN             fn;
  Mat            A;
  PetscInt       i,j,n=10;
  PetscScalar    *As;
  PetscViewer    viewer;
  PetscBool      verbose,inplace,random,triang;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-random",&random));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-triang",&triang));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Matrix logarithm, n=%" PetscInt_FMT ".\n",n));

  /* Create logarithm function object */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&fn));
  PetscCall(FNSetType(fn,FNLOG));
  PetscCall(FNSetFromOptions(fn));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNView(fn,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Create matrices */
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  if (random) PetscCall(MatSetRandom(A,NULL));
  else {
    /* Fill A with a non-symmetric Toeplitz matrix */
    PetscCall(MatDenseGetArray(A,&As));
    for (i=0;i<n;i++) As[i+i*n]=2.0;
    for (j=1;j<3;j++) {
      for (i=0;i<n-j;i++) {
        As[i+(i+j)*n]=1.0;
        if (!triang) As[(i+j)+i*n]=-1.0;
      }
    }
    As[(n-1)*n] = -5.0;
    As[0] = 2.01;
    PetscCall(MatDenseRestoreArray(A,&As));
  }
  PetscCall(TestMatLog(fn,A,viewer,verbose,inplace));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&fn));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: grep -v "computing matrix functions"
      output_file: output/test13_1.out
      test:
         suffix: 1
         args: -fn_scale .04,2 -n 75 -inplace {{0 1}}
         requires: c99_complex !__float128
      test:
         suffix: 1_triang
         args: -fn_scale .04,2 -n 75 -triang
         requires: c99_complex !__float128
      test:
         suffix: 1_random
         args: -fn_scale .02,2 -n 75 -random
         requires: complex !__float128
         filter_output: sed -e 's/04/02/'

TEST*/
