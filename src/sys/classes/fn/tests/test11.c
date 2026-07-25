/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Define the function

        f(x) = (exp(x)-1)/x    (the phi_1 function)

   with the following tree:

            f(x)                  f(x)              (combined by division)
           /    \                 p(x) = x          (polynomial)
        a(x)    p(x)              a(x)              (combined by addition)
       /    \                     e(x) = exp(x)     (exponential)
     e(x)   c(x)                  c(x) = -1         (constant)
*/

static char help[] = "Another test of a combined function.\n\n";

#include <slepcfn.h>

/*
   Compute matrix function B = A\(exp(A)-I)
 */
PetscErrorCode TestMatCombine(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscBool      set,flg;
  PetscInt       n;
  Mat            F,Acopy;
  Vec            v,f0;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&F));
  PetscCall(PetscObjectSetName((PetscObject)F,"F"));
  /* compute matrix function */
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
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) - - - - - - -\n"));
    PetscCall(MatView(F,viewer));
  }
  /* print matrix norm for checking */
  PetscCall(MatNorm(F,NORM_1,&nrm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"The 1-norm of f(A) is %6.3f\n",(double)nrm));
  /* check FNEvaluateFunctionMatVec() */
  PetscCall(MatCreateVecs(A,&v,&f0));
  PetscCall(MatGetColumnVector(F,f0,0));
  PetscCall(FNEvaluateFunctionMatVec(fn,A,v));
  PetscCall(VecAXPY(v,-1.0,f0));
  PetscCall(VecNorm(v,NORM_2,&nrm));
  if (nrm>100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm));
  PetscCall(MatDestroy(&F));
  PetscCall(VecDestroy(&v));
  PetscCall(VecDestroy(&f0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  FN             f,p,a,e,c,f1,f2;
  FNCombineType  ctype;
  Mat            A=NULL;
  PetscInt       i,j,n=10,np;
  PetscScalar    x,y,yp,*As,coeffs[10];
  char           strx[50],str[50];
  PetscViewer    viewer;
  PetscBool      verbose,inplace,matcuda;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-inplace",&inplace));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-matcuda",&matcuda));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Phi1 via a combined function, n=%" PetscInt_FMT ".\n",n));

  /* Create function */

  /* e(x) = exp(x) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&e));
  PetscCall(PetscObjectSetName((PetscObject)e,"e"));
  PetscCall(FNSetType(e,FNEXP));
  PetscCall(FNSetFromOptions(e));
  /* c(x) = -1 */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&c));
  PetscCall(PetscObjectSetName((PetscObject)c,"c"));
  PetscCall(FNSetType(c,FNRATIONAL));
  PetscCall(FNSetFromOptions(c));
  np = 1;
  coeffs[0] = -1.0;
  PetscCall(FNRationalSetNumerator(c,np,coeffs));
  /* a(x) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&a));
  PetscCall(PetscObjectSetName((PetscObject)a,"a"));
  PetscCall(FNSetType(a,FNCOMBINE));
  PetscCall(FNSetFromOptions(a));
  PetscCall(FNCombineSetChildren(a,FN_COMBINE_ADD,e,c));
  /* p(x) = x */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&p));
  PetscCall(PetscObjectSetName((PetscObject)p,"p"));
  PetscCall(FNSetType(p,FNRATIONAL));
  PetscCall(FNSetFromOptions(p));
  np = 2;
  coeffs[0] = 1.0; coeffs[1] = 0.0;
  PetscCall(FNRationalSetNumerator(p,np,coeffs));
  /* f(x) */
  PetscCall(FNCreate(PETSC_COMM_WORLD,&f));
  PetscCall(PetscObjectSetName((PetscObject)f,"f"));
  PetscCall(FNSetType(f,FNCOMBINE));
  PetscCall(FNSetFromOptions(f));
  PetscCall(FNCombineSetChildren(f,FN_COMBINE_DIVIDE,a,p));

  /* Set up viewer */
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(FNCombineGetChildren(f,&ctype,&f1,&f2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Two functions combined with division:\n"));
  PetscCall(FNView(f1,viewer));
  PetscCall(FNView(f2,viewer));
  if (verbose) PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));

  /* Scalar evaluation */
  x = 2.2;
  PetscCall(SlepcSNPrintfScalar(strx,sizeof(strx),x,PETSC_FALSE));
  PetscCall(FNEvaluateFunction(f,x,&y));
  PetscCall(FNEvaluateDerivative(f,x,&yp));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),y,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),yp,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str));

  /* Create matrices */
  if (matcuda) {
#if PetscDefined(HAVE_CUDA)
    PetscCall(MatCreateSeqDenseCUDA(PETSC_COMM_SELF,n,n,NULL,&A));
#endif
  } else PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A));
  PetscCall(PetscObjectSetName((PetscObject)A,"A"));

  /* Fill A with 1-D Laplacian matrix */
  PetscCall(MatDenseGetArray(A,&As));
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  j=1;
  for (i=0;i<n-j;i++) { As[i+(i+j)*n]=-1.0; As[(i+j)+i*n]=-1.0; }
  PetscCall(MatDenseRestoreArray(A,&As));
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE));
  PetscCall(TestMatCombine(f,A,viewer,verbose,inplace));

  /* Repeat with same matrix as non-symmetric */
  PetscCall(MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE));
  PetscCall(TestMatCombine(f,A,viewer,verbose,inplace));

  PetscCall(MatDestroy(&A));
  PetscCall(FNDestroy(&f));
  PetscCall(FNDestroy(&p));
  PetscCall(FNDestroy(&a));
  PetscCall(FNDestroy(&e));
  PetscCall(FNDestroy(&c));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      output_file: output/test11_1.out
      test:
         suffix: 1
      test:
         suffix: 1_cuda
         args: -matcuda
         requires: cuda
      test:
         suffix: 2
         args: -inplace
      test:
         suffix: 2_cuda
         args: -inplace -matcuda
         requires: cuda

TEST*/
