/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Shell spectral transformations with a non-injective mapping.\n\n"
  "Implements spectrum folding for the 2-D Laplacian, as in ex24.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

/* Context for spectrum folding spectral transformation */
typedef struct {
  Mat         A;
  Vec         w;
  PetscScalar target;
} FoldShellST;

/* Routines for shell spectral transformation */
PetscErrorCode STCreate_Fold(Mat,PetscScalar,FoldShellST**);
PetscErrorCode STApply_Fold(ST,Vec,Vec);
PetscErrorCode STDestroy_Fold(FoldShellST*);

int main (int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation context */
  FoldShellST    *fold;           /* user-defined spectral transform context */
  EPSType        type;
  PetscInt       N,n=10,m,i,j,II,Istart,Iend,nev;
  PetscBool      isShell,terse,flag;
  PetscScalar    target=1.1;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag));
  if (!flag) m = n;
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-target",&target,NULL));
  N = n*m;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum Folding via shell ST, N=%" PetscInt_FMT " (%" PetscInt_FMT "x%" PetscInt_FMT " grid) target=%3.2f\n\n",N,n,m,(double)PetscRealPart(target)));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the 5-point stencil Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) PetscCall(MatSetValue(A,II,II-n,-1.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A,II,II+n,-1.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A,II,II-1,-1.0,INSERT_VALUES));
    if (j<n-1) PetscCall(MatSetValue(A,II,II+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,II,II,4.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTarget(eps,target));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSHELL));
  PetscCall(EPSSetFromOptions(eps));

  /*
     Initialize shell spectral transformation
  */
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&isShell));
  if (isShell) {
    /* Change sorting criterion since this shell ST computes eigenvalues
       of the transformed operator closest to 0 */
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

    /* Create the context for the user-defined spectral transform */
    PetscCall(STCreate_Fold(A,target,&fold));
    PetscCall(STShellSetContext(st,fold));

    /* Set callback function for applying the operator (in this case we do not
       provide a back-transformation callback since the mapping is not one-to-one) */
    PetscCall(STShellSetApply(st,STApply_Fold));
    PetscCall(PetscObjectSetName((PetscObject)st,"STFOLD"));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  if (isShell) PetscCall(STDestroy_Fold(fold));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*
   STCreate_Fold - Creates the spectrum folding ST context.

   Input Parameters:
+  A - problem matrix
-  target - target value

   Output Parameter:
.  fold - user-defined spectral transformation context
*/
PetscErrorCode STCreate_Fold(Mat A,PetscScalar target,FoldShellST **fold)
{
  FoldShellST    *newctx;

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&newctx));
  newctx->A = A;
  PetscCall(PetscObjectReference((PetscObject)A));
  newctx->target = target;
  PetscCall(MatCreateVecs(A,&newctx->w,NULL));
  *fold = newctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   STApply_Fold - Applies the operator (A-target*I)^2 to a given vector.

   Input Parameters:
+  st - spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApply_Fold(ST st,Vec x,Vec y)
{
  FoldShellST    *fold;
  PetscScalar    sigma;

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&fold));
  sigma = -fold->target;
  PetscCall(MatMult(fold->A,x,fold->w));
  PetscCall(VecAXPY(fold->w,sigma,x));
  PetscCall(MatMult(fold->A,fold->w,y));
  PetscCall(VecAXPY(y,sigma,fold->w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   STDestroy_Fold - This routine destroys the shell ST context.

   Input Parameter:
.  fold - user-defined spectral transformation context
*/
PetscErrorCode STDestroy_Fold(FoldShellST *fold)
{
  PetscFunctionBeginUser;
  PetscCall(MatDestroy(&fold->A));
  PetscCall(VecDestroy(&fold->w));
  PetscCall(PetscFree(fold));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   test:
      args: -m 11 -eps_nev 4 -terse
      suffix: 1
      requires: !single

TEST*/
