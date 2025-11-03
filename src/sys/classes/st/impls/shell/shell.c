/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This provides a simple shell interface for programmers to create
   their own spectral transformations without writing much interface code
*/

#include <slepc/private/stimpl.h>        /*I "slepcst.h" I*/

typedef struct {
  void           *ctx;                       /* user provided context */
  PetscErrorCode (*apply)(ST,Vec,Vec);
  PetscErrorCode (*applytrans)(ST,Vec,Vec);
  PetscErrorCode (*applyhermtrans)(ST,Vec,Vec);
  PetscErrorCode (*backtransform)(ST,PetscInt n,PetscScalar*,PetscScalar*);
} ST_SHELL;

/*@C
   STShellGetContext - Returns the user-provided context associated with an `STSHELL`.

   Not Collective

   Input Parameter:
.  st - spectral transformation context

   Output Parameter:
.  ctx - the user provided context

   Level: advanced

.seealso: [](ch:st), `STSHELL`, `STShellSetContext()`
@*/
PetscErrorCode STShellGetContext(ST st,void *ctx)
{
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(ctx,2);
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&flg));
  if (!flg) *(void**)ctx = NULL;
  else      *(void**)ctx = ((ST_SHELL*)st->data)->ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STShellSetContext - Sets the user-defined context for an `STSHELL`.

   Logically Collective

   Input Parameters:
+  st - the shell ST
-  ctx - the context

   Level: advanced

   Fortran Notes:
   To use this from Fortran you must write a Fortran interface definition
   for this function that tells Fortran the Fortran derived data type that
   you are passing in as the `ctx` argument.

.seealso: [](ch:st), `STSHELL`, `STShellGetContext()`
@*/
PetscErrorCode STShellSetContext(ST st,void *ctx)
{
  ST_SHELL       *shell = (ST_SHELL*)st->data;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&flg));
  if (flg) shell->ctx = ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STApply_Shell(ST st,Vec x,Vec y)
{
  ST_SHELL         *shell = (ST_SHELL*)st->data;
  PetscObjectState instate,outstate;

  PetscFunctionBegin;
  PetscCheck(shell->apply,PetscObjectComm((PetscObject)st),PETSC_ERR_USER,"No apply() routine provided to Shell ST");
  PetscCall(VecGetState(y,&instate));
  PetscCallBack("STSHELL user function apply()",(*shell->apply)(st,x,y));
  PetscCall(VecGetState(y,&outstate));
  if (instate == outstate) {
    /* user forgot to increase the state of the output vector */
    PetscCall(PetscObjectStateIncrease((PetscObject)y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STApplyTranspose_Shell(ST st,Vec x,Vec y)
{
  ST_SHELL         *shell = (ST_SHELL*)st->data;
  PetscObjectState instate,outstate;

  PetscFunctionBegin;
  PetscCheck(shell->applytrans,PetscObjectComm((PetscObject)st),PETSC_ERR_USER,"No applytrans() routine provided to Shell ST");
  PetscCall(VecGetState(y,&instate));
  PetscCallBack("STSHELL user function applytrans()",(*shell->applytrans)(st,x,y));
  PetscCall(VecGetState(y,&outstate));
  if (instate == outstate) {
    /* user forgot to increase the state of the output vector */
    PetscCall(PetscObjectStateIncrease((PetscObject)y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode STApplyHermitianTranspose_Shell(ST st,Vec x,Vec y)
{
  ST_SHELL         *shell = (ST_SHELL*)st->data;
  PetscObjectState instate,outstate;
  Vec              w;

  PetscFunctionBegin;
  if (shell->applyhermtrans) {
    PetscCall(VecGetState(y,&instate));
    PetscCallBack("STSHELL user function applyhermtrans()",(*shell->applyhermtrans)(st,x,y));
    PetscCall(VecGetState(y,&outstate));
    if (instate == outstate) {
      /* user forgot to increase the state of the output vector */
      PetscCall(PetscObjectStateIncrease((PetscObject)y));
    }
  } else {
    PetscCall(VecDuplicate(x,&w));
    PetscCall(VecCopy(x,w));
    PetscCall(VecConjugate(w));
    PetscCall(STApplyTranspose_Shell(st,w,y));
    PetscCall(VecDestroy(&w));
    PetscCall(VecConjugate(y));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode STBackTransform_Shell(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  ST_SHELL       *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  if (shell->backtransform) PetscCallBack("STSHELL user function backtransform()",(*shell->backtransform)(st,n,eigr,eigi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   STIsInjective_Shell - Check if the user has provided the backtransform operation.
*/
PetscErrorCode STIsInjective_Shell(ST st,PetscBool* is)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  *is = shell->backtransform? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Shell(ST st)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(st->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApply_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyTranspose_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyHermitianTranspose_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetBackTransform_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STShellSetApply_Shell(ST st,STShellApplyFn *apply)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->apply = apply;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   STShellSetApply - Sets routine to use as the application of the
   operator to a vector in the user-defined spectral transformation.

   Logically Collective

   Input Parameters:
+  st    - the spectral transformation context
-  apply - the application-provided transformation routine

   Level: advanced

.seealso: [](ch:st), `STSHELL`, `STShellSetBackTransform()`, `STShellSetApplyTranspose()`, `STShellSetApplyHermitianTranspose()`
@*/
PetscErrorCode STShellSetApply(ST st,STShellApplyFn *apply)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscTryMethod(st,"STShellSetApply_C",(ST,STShellApplyFn*),(st,apply));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STShellSetApplyTranspose_Shell(ST st,STShellApplyTransposeFn *applytrans)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->applytrans = applytrans;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   STShellSetApplyTranspose - Sets routine to use as the application of the
   transposed operator to a vector in the user-defined spectral transformation.

   Logically Collective

   Input Parameters:
+  st    - the spectral transformation context
-  applytrans - the application-provided transformation routine

   Level: advanced

.seealso: [](ch:st), `STSHELL`, `STShellSetApply()`, `STShellSetBackTransform()`
@*/
PetscErrorCode STShellSetApplyTranspose(ST st,STShellApplyTransposeFn *applytrans)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscTryMethod(st,"STShellSetApplyTranspose_C",(ST,STShellApplyTransposeFn*),(st,applytrans));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode STShellSetApplyHermitianTranspose_Shell(ST st,STShellApplyHermitianTransposeFn *applyhermtrans)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->applyhermtrans = applyhermtrans;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/*@C
   STShellSetApplyHermitianTranspose - Sets routine to use as the application of the
   conjugate-transposed operator to a vector in the user-defined spectral transformation.

   Logically Collective

   Input Parameters:
+  st    - the spectral transformation context
-  applyhermtrans - the application-provided transformation routine

   Note:
   If configured with real scalars, this function has the same effect as `STShellSetApplyTranspose()`,
   so no need to call both.

   Level: advanced

.seealso: [](ch:st), `STSHELL`, `STShellSetApply()`, `STShellSetApplyTranspose()`, `STShellSetBackTransform()`
@*/
PetscErrorCode STShellSetApplyHermitianTranspose(ST st,STShellApplyHermitianTransposeFn *applyhermtrans)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscTryMethod(st,"STShellSetApplyHermitianTranspose_C",(ST,STShellApplyHermitianTransposeFn*),(st,applyhermtrans));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STShellSetBackTransform_Shell(ST st,STShellBackTransformFn *backtr)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->backtransform = backtr;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   STShellSetBackTransform - Sets the routine to be called after the
   eigensolution process has finished in order to transform back the
   computed eigenvalues.

   Logically Collective

   Input Parameters:
+  st     - the spectral transformation context
-  backtr - the application-provided backtransform routine

   Level: advanced

.seealso: [](ch:st), `STSHELL`, `STShellSetApply()`, `STShellSetApplyTranspose()`
@*/
PetscErrorCode STShellSetBackTransform(ST st,STShellBackTransformFn *backtr)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscTryMethod(st,"STShellSetBackTransform_C",(ST,STShellBackTransformFn*),(st,backtr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   STSHELL - STSHELL = "shell" - User-defined spectral transformation via callback
   functions for the application of the operator to a vector and (optionally) the
   backtransform operation.

   Level: beginner

   Note:
   In order to define a `shell` spectral transformation, the user has to provide
   the `apply` operation via `STShellSetApply()` and related functions, and
   optionally a `backtransform` operation via `STShellSetBackTransform()`, and
   in some cases a user-defined context containing relevant data via
   `STShellSetContext()`.

.seealso: [](ch:st), `ST`, `STType`, `STSetType()`, `STShellSetApply()`, `STShellSetBackTransform()`, `STShellSetContext()`
M*/

SLEPC_EXTERN PetscErrorCode STCreate_Shell(ST st)
{
  ST_SHELL       *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  st->data = (void*)ctx;

  st->usesksp = PETSC_FALSE;

  st->ops->apply           = STApply_Shell;
  st->ops->applytrans      = STApplyTranspose_Shell;
#if defined(PETSC_USE_COMPLEX)
  st->ops->applyhermtrans  = STApplyHermitianTranspose_Shell;
#else
  st->ops->applyhermtrans  = STApplyTranspose_Shell;
#endif
  st->ops->backtransform   = STBackTransform_Shell;
  st->ops->destroy         = STDestroy_Shell;

  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApply_C",STShellSetApply_Shell));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyTranspose_C",STShellSetApplyTranspose_Shell));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyHermitianTranspose_C",STShellSetApplyHermitianTranspose_Shell));
#else
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyHermitianTranspose_C",STShellSetApplyTranspose_Shell));
#endif
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STShellSetBackTransform_C",STShellSetBackTransform_Shell));
  PetscFunctionReturn(PETSC_SUCCESS);
}
