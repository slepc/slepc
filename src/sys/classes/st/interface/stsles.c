/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   ST interface routines related to the KSP object associated with it
*/

#include <slepc/private/stimpl.h>            /*I "slepcst.h" I*/

/*
   This is used to set a default type for the KSP and PC objects.
   It is called at STSetFromOptions (before KSPSetFromOptions)
   and also at STSetUp (in case STSetFromOptions was not called).
*/
PetscErrorCode STSetDefaultKSP(ST st)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
  PetscTryTypeMethod(st,setdefaultksp);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   This is done by all ST types except PRECOND.
   The default is an LU direct solver, or GMRES+Jacobi if matmode=shell.
*/
PetscErrorCode STSetDefaultKSP_Default(ST st)
{
  PC             pc;
  PCType         pctype;
  KSPType        ksptype;

  PetscFunctionBegin;
  PetscCall(KSPGetPC(st->ksp,&pc));
  PetscCall(KSPGetType(st->ksp,&ksptype));
  PetscCall(PCGetType(pc,&pctype));
  if (!pctype && !ksptype) {
    if (st->Pmat || st->Psplit) {
      PetscCall(KSPSetType(st->ksp,KSPBCGS));
      PetscCall(PCSetType(pc,PCBJACOBI));
    } else if (st->matmode == ST_MATMODE_SHELL) {
      PetscCall(KSPSetType(st->ksp,KSPGMRES));
      PetscCall(PCSetType(pc,PCJACOBI));
    } else {
      PetscCall(KSPSetType(st->ksp,KSPPREONLY));
      PetscCall(PCSetType(pc,st->asymm?PCCHOLESKY:PCLU));
    }
  }
  PetscCall(KSPSetErrorIfNotConverged(st->ksp,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatMult - Computes the matrix-vector product $y = T_k x$, where $T_k$ is
   the $k$-th (transformed) matrix of the spectral transformation.

   Neighbor-wise Collective

   Input Parameters:
+  st - the spectral transformation context
.  k  - index of matrix to use
-  x  - the vector to be multiplied

   Output Parameter:
.  y - the result

   Level: developer

.seealso: [](ch:st), `STMatMultTranspose()`, `STMatMultHermitianTranspose()`, `STGetMatrixTransformed()`
@*/
PetscErrorCode STMatMult(ST st,PetscInt k,Vec x,Vec y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  STCheckMatrices(st,1);
  PetscCheck(k>=0 && k<PetscMax(2,st->nmat),PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nmat);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(x));
  PetscCall(PetscLogEventBegin(ST_MatMult,st,x,y,0));
  if (!st->T[k]) PetscCall(VecCopy(x,y)); /* T[k]=NULL means identity matrix */
  else PetscCall(MatMult(st->T[k],x,y));
  PetscCall(PetscLogEventEnd(ST_MatMult,st,x,y,0));
  PetscCall(VecLockReadPop(x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatMultTranspose - Computes the matrix-vector product $y = T_k^T x$, where $T_k$
   is the $k$-th (transformed) matrix of the spectral transformation.

   Neighbor-wise Collective

   Input Parameters:
+  st - the spectral transformation context
.  k  - index of matrix to use
-  x  - the vector to be multiplied

   Output Parameter:
.  y - the result

   Level: developer

.seealso: [](ch:st), `STMatMult()`, `STMatMultHermitianTranspose()`, `STGetMatrixTransformed()`
@*/
PetscErrorCode STMatMultTranspose(ST st,PetscInt k,Vec x,Vec y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  STCheckMatrices(st,1);
  PetscCheck(k>=0 && k<PetscMax(2,st->nmat),PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nmat);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(x));
  PetscCall(PetscLogEventBegin(ST_MatMultTranspose,st,x,y,0));
  if (!st->T[k]) PetscCall(VecCopy(x,y)); /* T[k]=NULL means identity matrix */
  else PetscCall(MatMultTranspose(st->T[k],x,y));
  PetscCall(PetscLogEventEnd(ST_MatMultTranspose,st,x,y,0));
  PetscCall(VecLockReadPop(x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatMultHermitianTranspose - Computes the matrix-vector product $y = T_k^* x$,
   where $T_k$ is the $k$-th (transformed) matrix of the spectral transformation.

   Neighbor-wise Collective

   Input Parameters:
+  st - the spectral transformation context
.  k  - index of matrix to use
-  x  - the vector to be multiplied

   Output Parameter:
.  y - the result

   Level: developer

.seealso: [](ch:st), `STMatMult()`, `STMatMultTranspose()`, `STGetMatrixTransformed()`
@*/
PetscErrorCode STMatMultHermitianTranspose(ST st,PetscInt k,Vec x,Vec y)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,k,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  STCheckMatrices(st,1);
  PetscCheck(k>=0 && k<PetscMax(2,st->nmat),PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"k must be between 0 and %" PetscInt_FMT,st->nmat);
  PetscCheck(x!=y,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  PetscCall(VecSetErrorIfLocked(y,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(x));
  PetscCall(PetscLogEventBegin(ST_MatMultTranspose,st,x,y,0));
  if (!st->T[k]) PetscCall(VecCopy(x,y)); /* T[k]=NULL means identity matrix */
  else PetscCall(MatMultHermitianTranspose(st->T[k],x,y));
  PetscCall(PetscLogEventEnd(ST_MatMultTranspose,st,x,y,0));
  PetscCall(VecLockReadPop(x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatSolve - Solves $P x = b$, where $P$ is the preconditioner matrix of
   the spectral transformation, using a `KSP` object stored internally.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  b  - right hand side vector

   Output Parameter:
.  x - computed solution

   Level: developer

.seealso: [](ch:st), `STMatSolveTranspose()`, `STMatSolveHermitianTranspose()`, `STMatMatSolve()`
@*/
PetscErrorCode STMatSolve(ST st,Vec b,Vec x)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  STCheckMatrices(st,1);
  PetscCheck(x!=b,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and b must be different vectors");
  PetscCall(VecSetErrorIfLocked(x,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(b));
  PetscCall(PetscLogEventBegin(ST_MatSolve,st,b,x,0));
  if (!st->P) PetscCall(VecCopy(b,x)); /* P=NULL means identity matrix */
  else PetscCall(KSPSolve(st->ksp,b,x));
  PetscCall(PetscLogEventEnd(ST_MatSolve,st,b,x,0));
  PetscCall(VecLockReadPop(b));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatMatSolve - Solves $P X = B$, where $P$ is the preconditioner matrix of
   the spectral transformation, using a `KSP` object stored internally.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  B  - right hand side vectors

   Output Parameter:
.  X - computed solutions

   Level: developer

.seealso: [](ch:st), `STMatSolve()`
@*/
PetscErrorCode STMatMatSolve(ST st,Mat B,Mat X)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(B,MAT_CLASSID,2);
  PetscValidHeaderSpecific(X,MAT_CLASSID,3);
  STCheckMatrices(st,1);

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(PetscLogEventBegin(ST_MatSolve,st,B,X,0));
  if (!st->P) PetscCall(MatCopy(B,X,SAME_NONZERO_PATTERN)); /* P=NULL means identity matrix */
  else PetscCall(KSPMatSolve(st->ksp,B,X));
  PetscCall(PetscLogEventEnd(ST_MatSolve,st,B,X,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatSolveTranspose - Solves $P^T x = b$, where $P$ is the preconditioner matrix of
   the spectral transformation, using a `KSP` object stored internally.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  b  - right hand side vector

   Output Parameter:
.  x - computed solution

   Level: developer

.seealso: [](ch:st), `STMatSolve()`
@*/
PetscErrorCode STMatSolveTranspose(ST st,Vec b,Vec x)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  STCheckMatrices(st,1);
  PetscCheck(x!=b,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and b must be different vectors");
  PetscCall(VecSetErrorIfLocked(x,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(b));
  PetscCall(PetscLogEventBegin(ST_MatSolveTranspose,st,b,x,0));
  if (!st->P) PetscCall(VecCopy(b,x)); /* P=NULL means identity matrix */
  else PetscCall(KSPSolveTranspose(st->ksp,b,x));
  PetscCall(PetscLogEventEnd(ST_MatSolveTranspose,st,b,x,0));
  PetscCall(VecLockReadPop(b));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STMatSolveHermitianTranspose - Solves $P^* x = b$, where $P$ is the preconditioner matrix of
   the spectral transformation, using a `KSP` object stored internally.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  b  - right hand side vector

   Output Parameter:
.  x - computed solution

   Level: developer

.seealso: [](ch:st), `STMatSolve()`
@*/
PetscErrorCode STMatSolveHermitianTranspose(ST st,Vec b,Vec x)
{
  Vec  w;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  STCheckMatrices(st,1);
  PetscCheck(x!=b,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and b must be different vectors");
  PetscCall(VecSetErrorIfLocked(x,3));

  if (st->state!=ST_STATE_SETUP) PetscCall(STSetUp(st));
  PetscCall(VecLockReadPush(b));
  PetscCall(PetscLogEventBegin(ST_MatSolveTranspose,st,b,x,0));
  if (!st->P) PetscCall(VecCopy(b,x)); /* P=NULL means identity matrix */
  else {
    PetscCall(VecDuplicate(b,&w));
    PetscCall(VecCopy(b,w));
    PetscCall(VecConjugate(w));
    PetscCall(KSPSolveTranspose(st->ksp,w,x));
    PetscCall(VecDestroy(&w));
    PetscCall(VecConjugate(x));
  }
  PetscCall(PetscLogEventEnd(ST_MatSolveTranspose,st,b,x,0));
  PetscCall(VecLockReadPop(b));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STCheckFactorPackage(ST st)
{
  PC             pc;
  PetscMPIInt    size;
  PetscBool      flg;
  MatSolverType  stype;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)st),&size));
  if (size==1) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(KSPGetPC(st->ksp,&pc));
  PetscCall(PCFactorGetMatSolverType(pc,&stype));
  if (stype) {   /* currently selected PC is a factorization */
    PetscCall(PetscStrcmp(stype,MATSOLVERPETSC,&flg));
    PetscCheck(!flg,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"You chose to solve linear systems with a factorization, but in parallel runs you need to select an external package; see the users guide for details");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STSetKSP - Sets the `KSP` object associated with the spectral
   transformation.

   Collective

   Input Parameters:
+  st   - the spectral transformation context
-  ksp  - the linear solver context

   Level: advanced

.seealso: [](ch:st), `STGetKSP()`
@*/
PetscErrorCode STSetKSP(ST st,KSP ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(st,1,ksp,2);
  STCheckNotSeized(st,1);
  PetscCall(PetscObjectReference((PetscObject)ksp));
  PetscCall(KSPDestroy(&st->ksp));
  st->ksp = ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STGetKSP - Gets the KSP object associated with the spectral
   transformation.

   Collective

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  ksp  - the linear system context

   Level: intermediate

.seealso: [](ch:st), `STSetKSP()`
@*/
PetscErrorCode STGetKSP(ST st,KSP* ksp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(ksp,2);
  if (!st->ksp) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)st),&st->ksp));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)st->ksp,(PetscObject)st,1));
    PetscCall(KSPSetOptionsPrefix(st->ksp,((PetscObject)st)->prefix));
    PetscCall(KSPAppendOptionsPrefix(st->ksp,"st_"));
    PetscCall(PetscObjectSetOptions((PetscObject)st->ksp,((PetscObject)st)->options));
    PetscCall(KSPSetTolerances(st->ksp,SLEPC_DEFAULT_TOL,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
  }
  *ksp = st->ksp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode STCheckNullSpace_Default(ST st,BV V)
{
  PetscInt       nc,i,c;
  PetscReal      norm;
  Vec            *T,w,vi;
  Mat            A;
  PC             pc;
  MatNullSpace   nullsp;

  PetscFunctionBegin;
  PetscCall(BVGetNumConstraints(V,&nc));
  PetscCall(PetscMalloc1(nc,&T));
  if (!st->ksp) PetscCall(STGetKSP(st,&st->ksp));
  PetscCall(KSPGetPC(st->ksp,&pc));
  PetscCall(PCGetOperators(pc,&A,NULL));
  PetscCall(MatCreateVecs(A,NULL,&w));
  c = 0;
  for (i=0;i<nc;i++) {
    PetscCall(BVGetColumn(V,-nc+i,&vi));
    PetscCall(MatMult(A,vi,w));
    PetscCall(VecNorm(w,NORM_2,&norm));
    if (norm < 10.0*PETSC_SQRT_MACHINE_EPSILON) {
      PetscCall(PetscInfo(st,"Vector %" PetscInt_FMT " included in the nullspace of OP, norm=%g\n",i,(double)norm));
      PetscCall(BVCreateVec(V,T+c));
      PetscCall(VecCopy(vi,T[c]));
      PetscCall(VecNormalize(T[c],NULL));
      c++;
    } else PetscCall(PetscInfo(st,"Vector %" PetscInt_FMT " discarded as possible nullspace of OP, norm=%g\n",i,(double)norm));
    PetscCall(BVRestoreColumn(V,-nc+i,&vi));
  }
  PetscCall(VecDestroy(&w));
  if (c>0) {
    PetscCall(MatNullSpaceCreate(PetscObjectComm((PetscObject)st),PETSC_FALSE,c,T,&nullsp));
    PetscCall(MatSetNullSpace(A,nullsp));
    PetscCall(MatNullSpaceDestroy(&nullsp));
    PetscCall(VecDestroyVecs(c,&T));
  } else PetscCall(PetscFree(T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STCheckNullSpace - Tests if constraint vectors are nullspace vectors.

   Collective

   Input Parameters:
+  st - the spectral transformation context
-  V  - basis vectors to be checked

   Notes:
   Given a basis vectors object, this function tests each of its constraint
   vectors to be a nullspace vector of the coefficient matrix of the
   associated `KSP` object. All these nullspace vectors are passed to the `KSP`
   object with `MatSetNullSpace()`.

   This function allows handling singular pencils (i.e., when $A$ and $B$ share
   a common nullspace), and solving some problems in which the nullspace is important.

   Level: developer

.seealso: [](ch:st), `EPSSetDeflationSpace()`, `MatSetNullSpace()`
@*/
PetscErrorCode STCheckNullSpace(ST st,BV V)
{
  PetscInt       nc;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);
  PetscValidType(st,1);
  PetscCheckSameComm(st,1,V,2);
  PetscCheck(st->state,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONGSTATE,"Must call STSetUp() first");

  PetscCall(BVGetNumConstraints(V,&nc));
  if (nc) PetscTryTypeMethod(st,checknullspace,V);
  PetscFunctionReturn(PETSC_SUCCESS);
}
