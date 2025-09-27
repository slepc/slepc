/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the use of shell spectral transformations.\n\n"
  "The problem to be solved is the same as ex1.c and"
  "corresponds to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>

/* Define context for user-provided spectral transformation */
typedef struct {
  KSP        ksp;
} SampleShellST;

/* Declare routines for user-provided spectral transformation */
PetscErrorCode STCreate_User(SampleShellST**);
PetscErrorCode STSetUp_User(SampleShellST*,ST);
PetscErrorCode STApply_User(ST,Vec,Vec);
PetscErrorCode STApplyTranspose_User(ST,Vec,Vec);
#if defined(PETSC_USE_COMPLEX)
PetscErrorCode STApplyHermitianTranspose_User(ST,Vec,Vec);
#endif
PetscErrorCode STBackTransform_User(ST,PetscInt,PetscScalar*,PetscScalar*);
PetscErrorCode STDestroy_User(SampleShellST*);

int main (int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation context */
  SampleShellST  *shell;          /* user-defined spectral transform context */
  EPSType        type;
  PetscInt       n=30,i,Istart,Iend,nev;
  PetscBool      isShell,terse;
  PetscBool      set_ht=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem (shell-enabled), n=%" PetscInt_FMT "\n\n",n));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-set_ht",&set_ht,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /*
     Initialize shell spectral transformation if selected by user
  */
  PetscCall(EPSGetST(eps,&st));
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&isShell));
  if (isShell) {
    /* Change sorting criterion since this ST example computes values
       closest to 0 */
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

    /* (Required) Create a context for the user-defined spectral transform;
       this context can be defined to contain any application-specific data. */
    PetscCall(STCreate_User(&shell));
    PetscCall(STShellSetContext(st,shell));

    /* (Required) Set the user-defined routine for applying the operator */
    PetscCall(STShellSetApply(st,STApply_User));

    /* (Optional) Set the user-defined routine for applying the transposed operator */
    PetscCall(STShellSetApplyTranspose(st,STApplyTranspose_User));

    /* (Optional) Set the user-defined routine for applying the conjugate-transposed operator */
#if defined(PETSC_USE_COMPLEX)
    if (set_ht) PetscCall(STShellSetApplyHermitianTranspose(st,STApplyHermitianTranspose_User));
#endif

    /* (Optional) Set the user-defined routine for back-transformation */
    PetscCall(STShellSetBackTransform(st,STBackTransform_User));

    /* (Optional) Set a name for the transformation, used for STView() */
    PetscCall(PetscObjectSetName((PetscObject)st,"MyTransformation"));

    /* (Optional) Do any setup required for the new transformation */
    PetscCall(STSetUp_User(shell,st));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
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
  if (isShell) PetscCall(STDestroy_User(shell));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/***********************************************************************/
/*     Routines for a user-defined shell spectral transformation       */
/***********************************************************************/

/*
   STCreate_User - This routine creates a user-defined
   spectral transformation context.

   Output Parameter:
.  shell - user-defined spectral transformation context
*/
PetscErrorCode STCreate_User(SampleShellST **shell)
{
  SampleShellST  *newctx;

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&newctx));
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&newctx->ksp));
  PetscCall(KSPAppendOptionsPrefix(newctx->ksp,"st_"));
  *shell = newctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* ------------------------------------------------------------------- */
/*
   STSetUp_User - This routine sets up a user-defined
   spectral transformation context.

   Input Parameters:
+  shell - user-defined spectral transformation context
-  st    - spectral transformation context containing the operator matrices

   Output Parameter:
.  shell - fully set up user-defined transformation context

   Notes:
   In this example, the user-defined transformation is simply OP=A^-1.
   Therefore, the eigenpairs converge in reversed order. The KSP object
   used for the solution of linear systems with A is handled via the
   user-defined context SampleShellST.
*/
PetscErrorCode STSetUp_User(SampleShellST *shell,ST st)
{
  Mat            A;

  PetscFunctionBeginUser;
  PetscCall(STGetMatrix(st,0,&A));
  PetscCall(KSPSetOperators(shell->ksp,A,A));
  PetscCall(KSPSetFromOptions(shell->ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* ------------------------------------------------------------------- */
/*
   STApply_User - This routine demonstrates the use of a
   user-provided spectral transformation.

   Input Parameters:
+  st - spectral transformation context
-  x - input vector

   Output Parameter:
.  y - output vector

   Notes:
   The transformation implemented in this code is just OP=A^-1 and
   therefore it is of little use, merely as an example of working with
   a STSHELL.
*/
PetscErrorCode STApply_User(ST st,Vec x,Vec y)
{
  SampleShellST  *shell;

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&shell));
  PetscCall(KSPSolve(shell->ksp,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* ------------------------------------------------------------------- */
/*
   STApplyTranspose_User - This is not required unless using a two-sided
   eigensolver.

   Input Parameters:
+  st - spectral transformation context
-  x - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApplyTranspose_User(ST st,Vec x,Vec y)
{
  SampleShellST  *shell;

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&shell));
  PetscCall(KSPSolveTranspose(shell->ksp,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#if defined(PETSC_USE_COMPLEX)
/* ------------------------------------------------------------------- */
/*
   STApplyHermitianTranspose_User - This is not required unless using a two-sided
   eigensolver in complex scalars.

   Input Parameters:
+  st - spectral transformation context
-  x - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApplyHermitianTranspose_User(ST st,Vec x,Vec y)
{
  SampleShellST  *shell;
  Vec            w;

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&shell));
  PetscCall(VecDuplicate(x,&w));
  PetscCall(VecCopy(x,w));
  PetscCall(VecConjugate(w));
  PetscCall(KSPSolveTranspose(shell->ksp,w,y));
  PetscCall(VecConjugate(y));
  PetscCall(VecDestroy(&w));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
/* ------------------------------------------------------------------- */
/*
   STBackTransform_User - This routine demonstrates the use of a
   user-provided spectral transformation.

   Input Parameters:
+  st - spectral transformation context
-  n  - number of eigenvalues to transform

   Input/Output Parameters:
+  eigr - pointer to real part of eigenvalues
-  eigi - pointer to imaginary part of eigenvalues

   Notes:
   This code implements the back transformation of eigenvalues in
   order to retrieve the eigenvalues of the original problem. In this
   example, simply set k_i = 1/k_i.
*/
PetscErrorCode STBackTransform_User(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt j;

  PetscFunctionBeginUser;
  for (j=0;j<n;j++) {
    eigr[j] = 1.0 / eigr[j];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* ------------------------------------------------------------------- */
/*
   STDestroy_User - This routine destroys a user-defined
   spectral transformation context.

   Input Parameter:
.  shell - user-defined spectral transformation context
*/
PetscErrorCode STDestroy_User(SampleShellST *shell)
{
  PetscFunctionBeginUser;
  PetscCall(KSPDestroy(&shell->ksp));
  PetscCall(PetscFree(shell));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -eps_nev 5 -eps_non_hermitian -terse
      output_file: output/ex10_1.out
      test:
         suffix: 1_sinvert
         args: -st_type sinvert
         requires: !single
      test:
         suffix: 1_sinvert_twoside
         args: -st_type sinvert -eps_balance twoside -set_ht {{0 1}}
         requires: !single
      test:
         suffix: 1_shell
         args: -st_type shell
         requires: !single
      test:
         suffix: 1_shell_twoside
         args: -st_type shell -eps_balance twoside -set_ht {{0 1}}

TEST*/
