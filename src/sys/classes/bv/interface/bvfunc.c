/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV (basis vectors) interface routines, callable by users
*/

#include <slepc/private/bvimpl.h>            /*I "slepcbv.h" I*/

PetscClassId     BV_CLASSID = 0;
PetscLogEvent    BV_Create = 0,BV_Copy = 0,BV_Mult = 0,BV_MultVec = 0,BV_MultInPlace = 0,BV_Dot = 0,BV_DotVec = 0,BV_Orthogonalize = 0,BV_OrthogonalizeVec = 0,BV_Scale = 0,BV_Norm = 0,BV_NormVec = 0,BV_Normalize = 0,BV_SetRandom = 0,BV_MatMult = 0,BV_MatMultVec = 0,BV_MatProject = 0,BV_SVDAndRank = 0;
static PetscBool BVPackageInitialized = PETSC_FALSE;
MPI_Op MPIU_TSQR = 0,MPIU_LAPY2;

const char *BVOrthogTypes[] = {"CGS","MGS","BVOrthogType","BV_ORTHOG_",NULL};
const char *BVOrthogRefineTypes[] = {"IFNEEDED","NEVER","ALWAYS","BVOrthogRefineType","BV_ORTHOG_REFINE_",NULL};
const char *BVOrthogBlockTypes[] = {"GS","CHOL","TSQR","TSQRCHOL","SVQB","BVOrthogBlockType","BV_ORTHOG_BLOCK_",NULL};
const char *BVMatMultTypes[] = {"VECS","MAT","MAT_SAVE","BVMatMultType","BV_MATMULT_",NULL};
const char *BVSVDMethods[] = {"REFINE","QR","QR_CAA","BVSVDMethod","BV_SVD_METHOD_",NULL};

/*@C
   BVFinalizePackage - This function destroys everything in the SLEPc interface
   to the `BV` package. It is called from `SlepcFinalize()`.

   Level: developer

.seealso: [](sec:bv), `SlepcFinalize()`, `BVInitializePackage()`
@*/
PetscErrorCode BVFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&BVList));
  PetscCallMPI(MPI_Op_free(&MPIU_TSQR));
  PetscCallMPI(MPI_Op_free(&MPIU_LAPY2));
  BVPackageInitialized = PETSC_FALSE;
  BVRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVInitializePackage - This function initializes everything in the `BV` package.
   It is called from `PetscDLLibraryRegister_slepc()` when using dynamic libraries, and
   on the first call to `BVCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](sec:bv), `BV`, `SlepcInitialize()`, `BVFinalizePackage()`
@*/
PetscErrorCode BVInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (BVPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  BVPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Basis Vectors",&BV_CLASSID));
  /* Register Constructors */
  PetscCall(BVRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("BVCreate",BV_CLASSID,&BV_Create));
  PetscCall(PetscLogEventRegister("BVCopy",BV_CLASSID,&BV_Copy));
  PetscCall(PetscLogEventRegister("BVMult",BV_CLASSID,&BV_Mult));
  PetscCall(PetscLogEventRegister("BVMultVec",BV_CLASSID,&BV_MultVec));
  PetscCall(PetscLogEventRegister("BVMultInPlace",BV_CLASSID,&BV_MultInPlace));
  PetscCall(PetscLogEventRegister("BVDot",BV_CLASSID,&BV_Dot));
  PetscCall(PetscLogEventRegister("BVDotVec",BV_CLASSID,&BV_DotVec));
  PetscCall(PetscLogEventRegister("BVOrthogonalize",BV_CLASSID,&BV_Orthogonalize));
  PetscCall(PetscLogEventRegister("BVOrthogonalizeV",BV_CLASSID,&BV_OrthogonalizeVec));
  PetscCall(PetscLogEventRegister("BVScale",BV_CLASSID,&BV_Scale));
  PetscCall(PetscLogEventRegister("BVNorm",BV_CLASSID,&BV_Norm));
  PetscCall(PetscLogEventRegister("BVNormVec",BV_CLASSID,&BV_NormVec));
  PetscCall(PetscLogEventRegister("BVNormalize",BV_CLASSID,&BV_Normalize));
  PetscCall(PetscLogEventRegister("BVSetRandom",BV_CLASSID,&BV_SetRandom));
  PetscCall(PetscLogEventRegister("BVMatMult",BV_CLASSID,&BV_MatMult));
  PetscCall(PetscLogEventRegister("BVMatMultVec",BV_CLASSID,&BV_MatMultVec));
  PetscCall(PetscLogEventRegister("BVMatProject",BV_CLASSID,&BV_MatProject));
  PetscCall(PetscLogEventRegister("BVSVDAndRank",BV_CLASSID,&BV_SVDAndRank));
  /* MPI reduction operation used in BVOrthogonalize */
  PetscCallMPI(MPI_Op_create(SlepcGivensPacked,PETSC_FALSE,&MPIU_TSQR));
  PetscCallMPI(MPI_Op_create(SlepcPythag,PETSC_TRUE,&MPIU_LAPY2));
  /* Process Info */
  classids[0] = BV_CLASSID;
  PetscCall(PetscInfoProcessClass("bv",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("bv",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(BV_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(BVFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVDestroy - Destroys a `BV` context that was created with `BVCreate()`.

   Collective

   Input Parameter:
.  bv - the basis vectors context

   Level: beginner

.seealso: [](sec:bv), `BVCreate()`
@*/
PetscErrorCode BVDestroy(BV *bv)
{
  PetscFunctionBegin;
  if (!*bv) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*bv,BV_CLASSID,1);
  PetscCheck(!(*bv)->lsplit,PetscObjectComm((PetscObject)*bv),PETSC_ERR_ARG_WRONGSTATE,"Must call BVRestoreSplit before destroying the BV");
  if (--((PetscObject)*bv)->refct > 0) { *bv = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscTryTypeMethod(*bv,destroy);
  PetscCall(PetscLayoutDestroy(&(*bv)->map));
  PetscCall(PetscFree((*bv)->vtype));
  PetscCall(MatDestroy(&(*bv)->matrix));
  PetscCall(VecDestroy(&(*bv)->Bx));
  PetscCall(VecDestroy(&(*bv)->buffer));
  PetscCall(BVDestroy(&(*bv)->cached));
  PetscCall(BVDestroy(&(*bv)->L));
  PetscCall(BVDestroy(&(*bv)->R));
  PetscCall(PetscFree((*bv)->work));
  PetscCall(PetscFree2((*bv)->h,(*bv)->c));
  PetscCall(VecDestroy(&(*bv)->omega));
  PetscCall(MatDestroy(&(*bv)->Acreate));
  PetscCall(MatDestroy(&(*bv)->Aget));
  PetscCall(MatDestroy(&(*bv)->Abuffer));
  PetscCall(PetscRandomDestroy(&(*bv)->rand));
  PetscCall(PetscHeaderDestroy(bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreate - Creates a basis vectors context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newbv - location to put the basis vectors context

   Note:
   A basis vectors object of type `BVTENSOR` should be created using
   the helper function `BVCreateTensor()` instead of this one.

   Level: beginner

.seealso: [](sec:bv), `BV`, `BVDuplicate()`, `BVDestroy()`, `BVCreateTensor()` , `BVCreateFromMat()`
@*/
PetscErrorCode BVCreate(MPI_Comm comm,BV *newbv)
{
  BV             bv;

  PetscFunctionBegin;
  PetscAssertPointer(newbv,2);
  PetscCall(BVInitializePackage());
  PetscCall(SlepcHeaderCreate(bv,BV_CLASSID,"BV","Basis Vectors","BV",comm,BVDestroy,BVView));

  bv->map          = NULL;
  bv->vtype        = NULL;
  bv->n            = -1;
  bv->N            = -1;
  bv->m            = 0;
  bv->l            = 0;
  bv->k            = 0;
  bv->nc           = 0;
  bv->ld           = 0;
  bv->orthog_type  = BV_ORTHOG_CGS;
  bv->orthog_ref   = BV_ORTHOG_REFINE_IFNEEDED;
  bv->orthog_eta   = 0.7071;
  bv->orthog_block = BV_ORTHOG_BLOCK_GS;
  bv->matrix       = NULL;
  bv->indef        = PETSC_FALSE;
  bv->vmm          = BV_MATMULT_MAT;
  bv->rrandom      = PETSC_FALSE;
  bv->deftol       = 10*PETSC_MACHINE_EPSILON;

  bv->buffer       = NULL;
  bv->Abuffer      = NULL;
  bv->Bx           = NULL;
  bv->xid          = 0;
  bv->xstate       = 0;
  bv->cv[0]        = NULL;
  bv->cv[1]        = NULL;
  bv->ci[0]        = -1;
  bv->ci[1]        = -1;
  bv->st[0]        = -1;
  bv->st[1]        = -1;
  bv->id[0]        = 0;
  bv->id[1]        = 0;
  bv->h            = NULL;
  bv->c            = NULL;
  bv->omega        = NULL;
  bv->defersfo     = PETSC_FALSE;
  bv->cached       = NULL;
  bv->bvstate      = 0;
  bv->L            = NULL;
  bv->R            = NULL;
  bv->lstate       = 0;
  bv->rstate       = 0;
  bv->lsplit       = 0;
  bv->issplit      = 0;
  bv->splitparent  = NULL;
  bv->rand         = NULL;
  bv->Acreate      = NULL;
  bv->Aget         = NULL;
  bv->cuda         = PETSC_FALSE;
  bv->hip          = PETSC_FALSE;
  bv->sfocalled    = PETSC_FALSE;
  bv->work         = NULL;
  bv->lwork        = 0;
  bv->data         = NULL;

  *newbv = bv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVCreateFromMat - Creates a basis vectors object from a dense `Mat` object.

   Collective

   Input Parameter:
.  A - a dense tall-skinny matrix

   Output Parameter:
.  bv - the new basis vectors context

   Notes:
   The matrix values are copied to the `BV` data storage, memory is not shared.

   The communicator of the `BV` object will be the same as `A`, and so will be
   the dimensions.

   Level: intermediate

.seealso: [](sec:bv), `BVCreate()`, `BVDestroy()`, `BVCreateMat()`
@*/
PetscErrorCode BVCreateFromMat(Mat A,BV *bv)
{
  PetscInt  n,N,k;
  VecType   vtype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_CLASSID,1);

  PetscCall(MatGetSize(A,&N,&k));
  PetscCall(MatGetLocalSize(A,&n,NULL));
  PetscCall(MatGetVecType(A,&vtype));
  PetscCall(BVCreate(PetscObjectComm((PetscObject)A),bv));
  PetscCall(BVSetSizes(*bv,n,N,k));
  PetscCall(BVSetVecType(*bv,vtype));

  (*bv)->Acreate = A;
  PetscCall(PetscObjectReference((PetscObject)A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVInsertVec - Insert a vector into the specified column.

   Logically Collective

   Input Parameters:
+  V - basis vectors
.  j - the column of `V` to be overwritten
-  w - the vector to be copied

   Level: intermediate

.seealso: [](sec:bv), `BVInsertVecs()`
@*/
PetscErrorCode BVInsertVec(BV V,PetscInt j,Vec w)
{
  PetscInt       n,N;
  Vec            v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,j,2);
  PetscValidHeaderSpecific(w,VEC_CLASSID,3);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscCheckSameComm(V,1,w,3);

  PetscCall(VecGetSize(w,&N));
  PetscCall(VecGetLocalSize(w,&n));
  PetscCheck(N==V->N && n==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ") do not match BV sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ")",N,n,V->N,V->n);
  PetscCheck(j>=-V->nc && j<V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %" PetscInt_FMT ", should be between %" PetscInt_FMT " and %" PetscInt_FMT,j,-V->nc,V->m-1);

  PetscCall(BVGetColumn(V,j,&v));
  PetscCall(VecCopy(w,v));
  PetscCall(BVRestoreColumn(V,j,&v));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVInsertVecs - Insert a set of vectors into the specified columns.

   Collective

   Input Parameters:
+  V - basis vectors
.  s - first column of `V` to be overwritten
.  W - array of vectors to be copied
-  orth - flag indicating if the vectors must be orthogonalized

   Input/Output Parameter:
.  m - number of input vectors, on output the number of linearly independent
       vectors

   Notes:
   Copies the contents of vectors `W` to `V(:,s:s+n)`. If the orthogonalization
   flag is set, then the vectors are copied one by one and then orthogonalized
   against the previous ones. If any of them is linearly dependent then it
   is discarded and the value of `m` is decreased.

   Level: intermediate

.seealso: [](sec:bv), `BVInsertVec()`, `BVOrthogonalizeColumn()`
@*/
PetscErrorCode BVInsertVecs(BV V,PetscInt s,PetscInt *m,Vec W[],PetscBool orth)
{
  PetscInt       n,N,i,ndep;
  PetscBool      lindep;
  PetscReal      norm;
  Vec            v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,s,2);
  PetscAssertPointer(m,3);
  PetscValidLogicalCollectiveInt(V,*m,3);
  if (!*m) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(*m>0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Number of vectors (given %" PetscInt_FMT ") cannot be negative",*m);
  PetscAssertPointer(W,4);
  PetscValidHeaderSpecific(*W,VEC_CLASSID,4);
  PetscValidLogicalCollectiveBool(V,orth,5);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscCheckSameComm(V,1,*W,4);

  PetscCall(VecGetSize(*W,&N));
  PetscCall(VecGetLocalSize(*W,&n));
  PetscCheck(N==V->N && n==V->n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ") do not match BV sizes (global %" PetscInt_FMT ", local %" PetscInt_FMT ")",N,n,V->N,V->n);
  PetscCheck(s>=0 && s<V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument s has wrong value %" PetscInt_FMT ", should be between 0 and %" PetscInt_FMT,s,V->m-1);
  PetscCheck(s+(*m)<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Too many vectors provided, there is only room for %" PetscInt_FMT,V->m);

  ndep = 0;
  for (i=0;i<*m;i++) {
    PetscCall(BVGetColumn(V,s+i-ndep,&v));
    PetscCall(VecCopy(W[i],v));
    PetscCall(BVRestoreColumn(V,s+i-ndep,&v));
    if (orth) {
      PetscCall(BVOrthogonalizeColumn(V,s+i-ndep,NULL,&norm,&lindep));
      if (norm==0.0 || lindep) {
        PetscCall(PetscInfo(V,"Removing linearly dependent vector %" PetscInt_FMT "\n",i));
        ndep++;
      } else PetscCall(BVScaleColumn(V,s+i-ndep,1.0/norm));
    }
  }
  *m -= ndep;
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVInsertConstraints - Insert a set of vectors as constraints.

   Collective

   Input Parameters:
+  V - basis vectors
-  C - set of vectors to be inserted as constraints

   Input/Output Parameter:
.  nc - number of input vectors, on output the number of linearly independent
       vectors

   Notes:
   The constraints are relevant only during orthogonalization. Constraint
   vectors span a subspace that is deflated in every orthogonalization
   operation, so they are intended for removing those directions from the
   orthogonal basis computed in regular `BV` columns.

   Constraints are not stored in regular `BV` columns, but in a special part of
   the storage. They can be accessed with negative indices in `BVGetColumn()`.

   This operation is DESTRUCTIVE, meaning that all data contained in the
   columns of `V` is lost. This is typically invoked just after creating the `BV`.
   Once a set of constraints has been set, it is not allowed to call this
   function again.

   The vectors are copied one by one and then orthogonalized against the
   previous ones. If any of them is linearly dependent then it is discarded
   and the value of `nc` is decreased. The behavior is similar to `BVInsertVecs()`.

   Level: advanced

.seealso: [](sec:bv), `BVInsertVecs()`, `BVOrthogonalizeColumn()`, `BVGetColumn()`, `BVGetNumConstraints()`
@*/
PetscErrorCode BVInsertConstraints(BV V,PetscInt *nc,Vec C[])
{
  PetscInt       msave;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscAssertPointer(nc,2);
  PetscValidLogicalCollectiveInt(V,*nc,2);
  if (!*nc) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCheck(*nc>0,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Number of constraints (given %" PetscInt_FMT ") cannot be negative",*nc);
  PetscAssertPointer(C,3);
  PetscValidHeaderSpecific(*C,VEC_CLASSID,3);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscCheckSameComm(V,1,*C,3);
  PetscCheck(!V->issplit,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"Operation not permitted for a BV obtained from BVGetSplit");
  PetscCheck(!V->nc,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"Constraints already present in this BV object");
  PetscCheck(V->ci[0]==-1 && V->ci[1]==-1,PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Cannot call BVInsertConstraints after BVGetColumn");

  msave = V->m;
  PetscCall(BVResize(V,*nc+V->m,PETSC_FALSE));
  PetscCall(BVInsertVecs(V,0,nc,C,PETSC_TRUE));
  V->nc = *nc;
  V->m  = msave;
  V->ci[0] = -V->nc-1;
  V->ci[1] = -V->nc-1;
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVSetOptionsPrefix - Sets the prefix used for searching for all
   `BV` options in the database.

   Logically Collective

   Input Parameters:
+  bv     - the basis vectors context
-  prefix - the prefix string to prepend to all `BV` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: [](sec:bv), `BVAppendOptionsPrefix()`, `BVGetOptionsPrefix()`
@*/
PetscErrorCode BVSetOptionsPrefix(BV bv,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)bv,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVAppendOptionsPrefix - Appends to the prefix used for searching for all
   `BV` options in the database.

   Logically Collective

   Input Parameters:
+  bv     - the basis vectors context
-  prefix - the prefix string to prepend to all `BV` option requests

   Note:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: [](sec:bv), `BVSetOptionsPrefix()`, `BVGetOptionsPrefix()`
@*/
PetscErrorCode BVAppendOptionsPrefix(BV bv,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)bv,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVGetOptionsPrefix - Gets the prefix used for searching for all
   `BV` options in the database.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  prefix - pointer to the prefix string used is returned

   Level: advanced

.seealso: [](sec:bv), `BVSetOptionsPrefix()`, `BVAppendOptionsPrefix()`
@*/
PetscErrorCode BVGetOptionsPrefix(BV bv,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)bv,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVView - Prints the `BV` data structure.

   Collective

   Input Parameters:
+  bv     - the basis vectors context
-  viewer - optional visualization context

   Notes:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Use `BVViewFromOptions()` to allow the user to select many different `PetscViewerType`
   and formats from the options database.

   Level: beginner

.seealso: [](sec:bv), `BVCreate()`
@*/
PetscErrorCode BVView(BV bv,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;
  const char        *orthname[2] = {"classical","modified"};
  const char        *refname[3] = {"if needed","never","always"};

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)bv),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)bv,viewer));
    PetscCall(PetscViewerGetFormat(viewer,&format));
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  %" PetscInt_FMT " columns of global length %" PetscInt_FMT "%s\n",bv->m,bv->N,bv->cuda?" (CUDA)":bv->hip?" (HIP)":""));
      if (bv->nc>0) PetscCall(PetscViewerASCIIPrintf(viewer,"  number of constraints: %" PetscInt_FMT "\n",bv->nc));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  vector orthogonalization method: %s Gram-Schmidt\n",orthname[bv->orthog_type]));
      switch (bv->orthog_ref) {
        case BV_ORTHOG_REFINE_IFNEEDED:
          PetscCall(PetscViewerASCIIPrintf(viewer,"  orthogonalization refinement: %s (eta: %g)\n",refname[bv->orthog_ref],(double)bv->orthog_eta));
          break;
        case BV_ORTHOG_REFINE_NEVER:
        case BV_ORTHOG_REFINE_ALWAYS:
          PetscCall(PetscViewerASCIIPrintf(viewer,"  orthogonalization refinement: %s\n",refname[bv->orthog_ref]));
          break;
      }
      PetscCall(PetscViewerASCIIPrintf(viewer,"  block orthogonalization method: %s\n",BVOrthogBlockTypes[bv->orthog_block]));
      if (bv->matrix) {
        if (bv->indef) PetscCall(PetscViewerASCIIPrintf(viewer,"  indefinite inner product\n"));
        else PetscCall(PetscViewerASCIIPrintf(viewer,"  non-standard inner product\n"));
        PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance for definite inner product: %g\n",(double)bv->deftol));
        PetscCall(PetscViewerASCIIPrintf(viewer,"  inner product matrix:\n"));
        PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO));
        PetscCall(PetscViewerASCIIPushTab(viewer));
        PetscCall(MatView(bv->matrix,viewer));
        PetscCall(PetscViewerASCIIPopTab(viewer));
        PetscCall(PetscViewerPopFormat(viewer));
      }
      switch (bv->vmm) {
        case BV_MATMULT_VECS:
          PetscCall(PetscViewerASCIIPrintf(viewer,"  doing matmult as matrix-vector products\n"));
          break;
        case BV_MATMULT_MAT:
          PetscCall(PetscViewerASCIIPrintf(viewer,"  doing matmult as a single matrix-matrix product\n"));
          break;
        case BV_MATMULT_MAT_SAVE:
          PetscCall(PetscViewerASCIIPrintf(viewer,"  mat_save is deprecated, use mat\n"));
          break;
      }
      if (bv->rrandom) PetscCall(PetscViewerASCIIPrintf(viewer,"  generating random vectors independent of the number of processes\n"));
    }
  }
  PetscTryTypeMethod(bv,view,viewer);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVViewFromOptions - View (print) a `BV` object based on values in the options database.

   Collective

   Input Parameters:
+  bv   - the basis vectors context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](sec:bv), `BVView()`, `BVCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode BVViewFromOptions(BV bv,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)bv,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   BVRegister - Adds a new storage format to the `BV` package.

   Not Collective

   Input Parameters:
+  name     - name of a new user-defined `BV`
-  function - routine to create context

   Notes:
   `BVRegister()` may be called multiple times to add several user-defined
   basis vectors.

   Level: advanced

.seealso: [](sec:bv), `BVRegisterAll()`
@*/
PetscErrorCode BVRegister(const char *name,PetscErrorCode (*function)(BV))
{
  PetscFunctionBegin;
  PetscCall(BVInitializePackage());
  PetscCall(PetscFunctionListAdd(&BVList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BVAllocateWork_Private(BV bv,PetscInt s)
{
  PetscFunctionBegin;
  if (s>bv->lwork) {
    PetscCall(PetscFree(bv->work));
    PetscCall(PetscMalloc1(s,&bv->work));
    bv->lwork = s;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_USE_DEBUG) && !defined(PETSC_CLANG_STATIC_ANALYZER)
/*
   SlepcDebugBVView - partially view a BV object, to be used from within a debugger.

     ini, end: columns to be viewed
     s: name of Matlab variable
     filename: optionally write output to a file
 */
PETSC_UNUSED PetscErrorCode SlepcDebugBVView(BV bv,PetscInt ini,PetscInt end,const char *s,const char *filename)
{
  PetscInt       N,m;
  PetscScalar    *array;

  PetscFunctionBegin;
  PetscCall(BVGetArray(bv,&array));
  PetscCall(BVGetSizes(bv,NULL,&N,&m));
  PetscCall(SlepcDebugViewMatrix(N,end-ini+1,array+ini*N,NULL,bv->ld,s,filename));
  PetscCall(BVRestoreArray(bv,&array));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
