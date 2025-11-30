/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Routines to set ST methods and options
*/

#include <slepc/private/stimpl.h>      /*I "slepcst.h" I*/

PetscBool         STRegisterAllCalled = PETSC_FALSE;
PetscFunctionList STList = NULL;

/*@
   STSetType - Selects the type of spectral transformation.

   Logically Collective

   Input Parameters:
+  st   - the spectral transformation context.
-  type - a known type

   Options Database Key:
.  -st_type \<type\> - sets the `ST` type

   Notes:
   Normally, it is best to call the `EPSSetFromOptions()` function and
   then set the `ST` type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different transformations.

   Level: beginner

.seealso: [](ch:st), `STGetType()`, `EPSSetType()`, `EPSSetFromOptions()`

@*/
PetscErrorCode STSetType(ST st,STType type)
{
  PetscErrorCode (*r)(ST);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)st,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);
  STCheckNotSeized(st,1);

  PetscCall(PetscFunctionListFind(STList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested ST type %s",type);

  PetscTryTypeMethod(st,destroy);
  PetscCall(PetscMemzero(st->ops,sizeof(struct _STOps)));

  st->state   = ST_STATE_INITIAL;
  st->opready = PETSC_FALSE;
  PetscCall(PetscObjectChangeTypeName((PetscObject)st,type));
  PetscCall((*r)(st));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetType - Gets the `ST` type name (as a string) from the `ST` context.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  type - name of the spectral transformation

   Level: intermediate

.seealso: [](ch:st), `STSetType()`

@*/
PetscErrorCode STGetType(ST st,STType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)st)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetFromOptions - Sets `ST` options from the options database.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Note:
   To see all options, run your program with the `-help` option.

   Level: beginner

.seealso: [](ch:st), `STSetOptionsPrefix()`
@*/
PetscErrorCode STSetFromOptions(ST st)
{
  PetscScalar    s;
  char           type[256];
  PetscBool      flg,bval;
  STMatMode      mode;
  MatStructure   mstr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscCall(STRegisterAll());
  PetscObjectOptionsBegin((PetscObject)st);
    PetscCall(PetscOptionsFList("-st_type","Spectral transformation","STSetType",STList,(char*)(((PetscObject)st)->type_name?((PetscObject)st)->type_name:STSHIFT),type,sizeof(type),&flg));
    if (flg) PetscCall(STSetType(st,type));
    else if (!((PetscObject)st)->type_name) PetscCall(STSetType(st,STSHIFT));

    PetscCall(PetscOptionsScalar("-st_shift","Value of the shift","STSetShift",st->sigma,&s,&flg));
    if (flg) PetscCall(STSetShift(st,s));

    PetscCall(PetscOptionsEnum("-st_matmode","Matrix mode for transformed matrices","STSetMatMode",STMatModes,(PetscEnum)st->matmode,(PetscEnum*)&mode,&flg));
    if (flg) PetscCall(STSetMatMode(st,mode));

    PetscCall(PetscOptionsEnum("-st_matstructure","Relation of the sparsity pattern of the matrices","STSetMatStructure",MatStructures,(PetscEnum)st->str,(PetscEnum*)&mstr,&flg));
    if (flg) PetscCall(STSetMatStructure(st,mstr));

    PetscCall(PetscOptionsBool("-st_transform","Whether transformed matrices are computed or not","STSetTransform",st->transform,&bval,&flg));
    if (flg) PetscCall(STSetTransform(st,bval));

    PetscTryTypeMethod(st,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)st,PetscOptionsObject));
  PetscOptionsEnd();

  if (st->usesksp) {
    PetscCall(STSetDefaultKSP(st));
    PetscCall(KSPSetFromOptions(st->ksp));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetMatStructure - Sets an internal `MatStructure` attribute to
   indicate which is the relation of the sparsity pattern of all `ST` matrices.

   Logically Collective

   Input Parameters:
+  st  - the spectral transformation context
-  str - either `SAME_NONZERO_PATTERN`, `DIFFERENT_NONZERO_PATTERN`,
         `SUBSET_NONZERO_PATTERN`, or `UNKNOWN_NONZERO_PATTERN`

   Options Database Key:
.  -st_matstructure \<str\> - the structure flag, one of `same`, `different`, `subset`, `unknown`

   Notes:
   If the sparsity pattern of the second matrix is equal or a subset of the
   pattern of the first matrix then it is recommended to set this attribute
   for efficiency reasons (in particular, for internal `MatAXPY()` operations).
   If not set, the default is `UNKNOWN_NONZERO_PATTERN`, in which case the patterns
   will be compared to determine if they are equal.

   This function has no effect in the case of standard eigenproblems.

   In case of polynomial eigenproblems, the flag applies to all matrices
   relative to the first one.

   Level: advanced

.seealso: [](ch:st), `STGetMatStructure()`, `STSetMatrices()`, `MatAXPY()`
@*/
PetscErrorCode STSetMatStructure(ST st,MatStructure str)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveEnum(st,str,2);
  switch (str) {
    case SAME_NONZERO_PATTERN:
    case DIFFERENT_NONZERO_PATTERN:
    case SUBSET_NONZERO_PATTERN:
    case UNKNOWN_NONZERO_PATTERN:
      st->str = str;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Invalid matrix structure flag");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetMatStructure - Gets the internal `MatStructure` attribute to
   indicate which is the relation of the sparsity pattern of the matrices.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  str - either `SAME_NONZERO_PATTERN`, `DIFFERENT_NONZERO_PATTERN`,
         `SUBSET_NONZERO_PATTERN`, or `UNKNOWN_NONZERO_PATTERN`

   Level: advanced

.seealso: [](ch:st), `STSetMatStructure()`, `STSetMatrices()`, `MatAXPY()`
@*/
PetscErrorCode STGetMatStructure(ST st,MatStructure *str)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(str,2);
  *str = st->str;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetMatMode - Sets a flag to indicate how the transformed matrices are
   being stored in the spectral transformations.

   Logically Collective

   Input Parameters:
+  st - the spectral transformation context
-  mode - the mode

   Options Database Key:
.  -st_matmode \<mode\> - set the mode, one of `copy`, `inplace`, `shell`

   Notes:
   By default (`ST_MATMODE_COPY`), a copy of matrix $A$ is made and then
   this copy is modified explicitly, e.g., $A \leftarrow A - \sigma B$.

   With `ST_MATMODE_INPLACE`, the original matrix $A$ is modified at `STSetUp()`
   and changes are reverted at the end of the computations. With respect to
   the previous one, this mode avoids a copy of matrix $A$. However, a
   drawback is that the recovered matrix might be slightly different
   from the original one (due to roundoff).

   With `ST_MATMODE_SHELL`, the solver works with an implicit shell
   matrix that represents the shifted matrix. This mode is the most efficient
   in creating the transformed matrix but it places serious limitations to the
   linear solves performed in each iteration of the eigensolver (typically,
   only iterative solvers with Jacobi preconditioning can be used).

   In the two first modes the efficiency of the computation can be controlled
   with `STSetMatStructure()`.

   Level: intermediate

.seealso: [](ch:st), `STSetMatrices()`, `STSetMatStructure()`, `STGetMatMode()`, `STMatMode`
@*/
PetscErrorCode STSetMatMode(ST st,STMatMode mode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveEnum(st,mode,2);
  if (st->matmode != mode) {
    STCheckNotSeized(st,1);
    st->matmode = mode;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetMatMode - Gets a flag that indicates how the transformed matrices
   are stored in spectral transformations.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  mode - the mode flag

   Level: intermediate

.seealso: [](ch:st), `STSetMatMode()`, `STMatMode`
@*/
PetscErrorCode STGetMatMode(ST st,STMatMode *mode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(mode,2);
  *mode = st->matmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetTransform - Sets a flag to indicate whether the transformed matrices are
   computed or not.

   Logically Collective

   Input Parameters:
+  st  - the spectral transformation context
-  flg - the boolean flag

   Options Database Key:
.  -st_transform - toggle the computation of matrices

   Notes:
   This flag is intended for the case of polynomial eigenproblems solved
   via linearization. If this flag is off (default) the spectral transformation
   is applied to the linearization (handled by the eigensolver), otherwise
   it is applied to the original problem.

   Level: developer

.seealso: [](ch:st), `STMatSolve()`, `STMatMult()`, `STSetMatStructure()`, `STGetTransform()`
@*/
PetscErrorCode STSetTransform(ST st,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveBool(st,flg,2);
  if (st->transform != flg) {
    st->transform = flg;
    st->state     = ST_STATE_INITIAL;
    st->opready   = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetTransform - Gets a flag that indicates whether the transformed
   matrices are computed or not.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  flg - the flag

   Level: developer

.seealso: [](ch:st), `STSetTransform()`
@*/
PetscErrorCode STGetTransform(ST st,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(flg,2);
  *flg = st->transform;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetStructured - Sets a flag to indicate that the application of the
   operator must be done taking into account its structure.

   Logically Collective

   Input Parameters:
+  st  - the spectral transformation context
-  flg - the boolean flag

   Note:
   This flag is intended for the case of structured eigenproblems. It is set
   internally by the solver, the user should not modify its value.

   Level: developer

.seealso: [](ch:st), `STApply()`, `STGetStructured()`
@*/
PetscErrorCode STSetStructured(ST st,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveBool(st,flg,2);
  if (st->structured != flg) {
    st->structured = flg;
    st->state      = ST_STATE_INITIAL;
    st->opready    = PETSC_FALSE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetStructured - Gets a flag that indicates if the application of the
   operator is done using its structure.

   Not Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  flg - the flag

   Level: developer

.seealso: [](ch:st), `STSetStructured()`
@*/
PetscErrorCode STGetStructured(ST st,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(flg,2);
  *flg = st->structured;
  PetscFunctionReturn(PETSC_SUCCESS);
}
