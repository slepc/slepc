/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic DS routines
*/

#include <slepc/private/dsimpl.h>      /*I "slepcds.h" I*/

PetscFunctionList DSList = NULL;
PetscBool         DSRegisterAllCalled = PETSC_FALSE;
PetscClassId      DS_CLASSID = 0;
PetscLogEvent     DS_Solve = 0,DS_Vectors = 0,DS_Synchronize = 0,DS_Other = 0;
static PetscBool  DSPackageInitialized = PETSC_FALSE;

const char *DSStateTypes[] = {"RAW","INTERMEDIATE","CONDENSED","TRUNCATED","DSStateType","DS_STATE_",NULL};
const char *DSParallelTypes[] = {"REDUNDANT","SYNCHRONIZED","DISTRIBUTED","DSParallelType","DS_PARALLEL_",NULL};
const char *DSMatName[DS_NUM_MAT] = {"A","B","C","T","D","Q","Z","X","Y","U","V","W","E0","E1","E2","E3","E4","E5","E6","E7","E8","E9"};
DSMatType  DSMatExtra[DS_NUM_EXTRA] = {DS_MAT_E0,DS_MAT_E1,DS_MAT_E2,DS_MAT_E3,DS_MAT_E4,DS_MAT_E5,DS_MAT_E6,DS_MAT_E7,DS_MAT_E8,DS_MAT_E9};

/*@C
   DSFinalizePackage - This function destroys everything in the SLEPc interface
   to the DS package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode DSFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&DSList));
  DSPackageInitialized = PETSC_FALSE;
  DSRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  DSInitializePackage - This function initializes everything in the DS package.
  It is called from PetscDLLibraryRegister() when using dynamic libraries, and
  on the first call to DSCreate() when using static libraries.

  Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode DSInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (DSPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  DSPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Direct Solver",&DS_CLASSID));
  /* Register Constructors */
  PetscCall(DSRegisterAll());
  /* Register Events */
  PetscCall(PetscLogEventRegister("DSSolve",DS_CLASSID,&DS_Solve));
  PetscCall(PetscLogEventRegister("DSVectors",DS_CLASSID,&DS_Vectors));
  PetscCall(PetscLogEventRegister("DSSynchronize",DS_CLASSID,&DS_Synchronize));
  PetscCall(PetscLogEventRegister("DSOther",DS_CLASSID,&DS_Other));
  /* Process Info */
  classids[0] = DS_CLASSID;
  PetscCall(PetscInfoProcessClass("ds",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("ds",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(DS_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(DSFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSCreate - Creates a DS context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newds - location to put the DS context

   Level: beginner

   Note:
   DS objects are not intended for normal users but only for
   advanced user that for instance implement their own solvers.

.seealso: DSDestroy(), DS
@*/
PetscErrorCode DSCreate(MPI_Comm comm,DS *newds)
{
  DS             ds;
  PetscInt       i;

  PetscFunctionBegin;
  PetscAssertPointer(newds,2);
  PetscCall(DSInitializePackage());
  PetscCall(SlepcHeaderCreate(ds,DS_CLASSID,"DS","Direct Solver (or Dense System)","DS",comm,DSDestroy,DSView));

  ds->state         = DS_STATE_RAW;
  ds->method        = 0;
  ds->compact       = PETSC_FALSE;
  ds->refined       = PETSC_FALSE;
  ds->extrarow      = PETSC_FALSE;
  ds->ld            = 0;
  ds->l             = 0;
  ds->n             = 0;
  ds->k             = 0;
  ds->t             = 0;
  ds->bs            = 1;
  ds->sc            = NULL;
  ds->pmode         = DS_PARALLEL_REDUNDANT;

  for (i=0;i<DS_NUM_MAT;i++) ds->omat[i] = NULL;
  ds->perm          = NULL;
  ds->data          = NULL;
  ds->scset         = PETSC_FALSE;
  ds->work          = NULL;
  ds->rwork         = NULL;
  ds->iwork         = NULL;
  ds->lwork         = 0;
  ds->lrwork        = 0;
  ds->liwork        = 0;

  *newds = ds;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetOptionsPrefix - Sets the prefix used for searching for all
   DS options in the database.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  prefix - the prefix string to prepend to all DS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: DSAppendOptionsPrefix()
@*/
PetscErrorCode DSSetOptionsPrefix(DS ds,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)ds,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSAppendOptionsPrefix - Appends to the prefix used for searching for all
   DS options in the database.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  prefix - the prefix string to prepend to all DS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: DSSetOptionsPrefix()
@*/
PetscErrorCode DSAppendOptionsPrefix(DS ds,const char *prefix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)ds,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetOptionsPrefix - Gets the prefix used for searching for all
   DS options in the database.

   Not Collective

   Input Parameters:
.  ds - the direct solver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: DSSetOptionsPrefix(), DSAppendOptionsPrefix()
@*/
PetscErrorCode DSGetOptionsPrefix(DS ds,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)ds,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetType - Selects the type for the DS object.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
-  type - a known type

   Level: intermediate

.seealso: DSGetType()
@*/
PetscErrorCode DSSetType(DS ds,DSType type)
{
  PetscErrorCode (*r)(DS);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)ds,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(DSList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested DS type %s",type);

  PetscTryTypeMethod(ds,destroy);
  PetscCall(DSReset(ds));
  PetscCall(PetscMemzero(ds->ops,sizeof(struct _DSOps)));

  PetscCall(PetscObjectChangeTypeName((PetscObject)ds,type));
  PetscCall((*r)(ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetType - Gets the DS type name (as a string) from the DS context.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  type - name of the direct solver

   Level: intermediate

.seealso: DSSetType()
@*/
PetscErrorCode DSGetType(DS ds,DSType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)ds)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSDuplicate - Creates a new direct solver object with the same options as
   an existing one.

   Collective

   Input Parameter:
.  ds - direct solver context

   Output Parameter:
.  dsnew - location to put the new DS

   Notes:
   DSDuplicate() DOES NOT COPY the matrices, and the new DS does not even have
   internal arrays allocated. Use DSAllocate() to use the new DS.

   The sorting criterion options are not copied, see DSSetSlepcSC().

   Level: intermediate

.seealso: DSCreate(), DSAllocate(), DSSetSlepcSC()
@*/
PetscErrorCode DSDuplicate(DS ds,DS *dsnew)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(dsnew,2);
  PetscCall(DSCreate(PetscObjectComm((PetscObject)ds),dsnew));
  if (((PetscObject)ds)->type_name) PetscCall(DSSetType(*dsnew,((PetscObject)ds)->type_name));
  (*dsnew)->method   = ds->method;
  (*dsnew)->compact  = ds->compact;
  (*dsnew)->refined  = ds->refined;
  (*dsnew)->extrarow = ds->extrarow;
  (*dsnew)->bs       = ds->bs;
  (*dsnew)->pmode    = ds->pmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetMethod - Selects the method to be used to solve the problem.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
-  meth - an index identifying the method

   Options Database Key:
.  -ds_method <meth> - Sets the method

   Level: intermediate

.seealso: DSGetMethod()
@*/
PetscErrorCode DSSetMethod(DS ds,PetscInt meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,meth,2);
  PetscCheck(meth>=0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The method must be a non-negative integer");
  PetscCheck(meth<=DS_MAX_SOLVE,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Too large value for the method");
  ds->method = meth;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetMethod - Gets the method currently used in the DS.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  meth - identifier of the method

   Level: intermediate

.seealso: DSSetMethod()
@*/
PetscErrorCode DSGetMethod(DS ds,PetscInt *meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(meth,2);
  *meth = ds->method;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetParallel - Selects the mode of operation in parallel runs.

   Logically Collective

   Input Parameters:
+  ds    - the direct solver context
-  pmode - the parallel mode

   Options Database Key:
.  -ds_parallel <mode> - Sets the parallel mode, 'redundant', 'synchronized'
   or 'distributed'

   Notes:
   In the 'redundant' parallel mode, all processes will make the computation
   redundantly, starting from the same data, and producing the same result.
   This result may be slightly different in the different processes if using a
   multithreaded BLAS library, which may cause issues in ill-conditioned problems.

   In the 'synchronized' parallel mode, only the first MPI process performs the
   computation and then the computed quantities are broadcast to the other
   processes in the communicator. This communication is not done automatically,
   an explicit call to DSSynchronize() is required.

   The 'distributed' parallel mode can be used in some DS types only, such
   as the contour integral method of DSNEP. In this case, every MPI process
   will be in charge of part of the computation.

   Level: advanced

.seealso: DSSynchronize(), DSGetParallel()
@*/
PetscErrorCode DSSetParallel(DS ds,DSParallelType pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(ds,pmode,2);
  ds->pmode = pmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetParallel - Gets the mode of operation in parallel runs.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  pmode - the parallel mode

   Level: advanced

.seealso: DSSetParallel()
@*/
PetscErrorCode DSGetParallel(DS ds,DSParallelType *pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(pmode,2);
  *pmode = ds->pmode;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetCompact - Switch to compact storage of matrices.

   Logically Collective

   Input Parameters:
+  ds   - the direct solver context
-  comp - a boolean flag

   Notes:
   Compact storage is used in some DS types such as DSHEP when the matrix
   is tridiagonal. This flag can be used to indicate whether the user
   provides the matrix entries via the compact form (the tridiagonal DS_MAT_T)
   or the non-compact one (DS_MAT_A).

   The default is PETSC_FALSE.

   Level: advanced

.seealso: DSGetCompact()
@*/
PetscErrorCode DSSetCompact(DS ds,PetscBool comp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveBool(ds,comp,2);
  if (ds->compact != comp && ds->ld) PetscTryTypeMethod(ds,setcompact,comp);
  ds->compact = comp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetCompact - Gets the compact storage flag.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  comp - the flag

   Level: advanced

.seealso: DSSetCompact()
@*/
PetscErrorCode DSGetCompact(DS ds,PetscBool *comp)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(comp,2);
  *comp = ds->compact;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetExtraRow - Sets a flag to indicate that the matrix has one extra
   row.

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
-  ext - a boolean flag

   Notes:
   In Krylov methods it is useful that the matrix representing the direct solver
   has one extra row, i.e., has dimension (n+1) x n. If this flag is activated, all
   transformations applied to the right of the matrix also affect this additional
   row. In that case, (n+1) must be less or equal than the leading dimension.

   The default is PETSC_FALSE.

   Level: advanced

.seealso: DSSolve(), DSAllocate(), DSGetExtraRow()
@*/
PetscErrorCode DSSetExtraRow(DS ds,PetscBool ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveBool(ds,ext,2);
  PetscCheck(!ext || ds->n==0 || ds->n!=ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"Cannot set extra row after setting n=ld");
  ds->extrarow = ext;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetExtraRow - Gets the extra row flag.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  ext - the flag

   Level: advanced

.seealso: DSSetExtraRow()
@*/
PetscErrorCode DSGetExtraRow(DS ds,PetscBool *ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(ext,2);
  *ext = ds->extrarow;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetRefined - Sets a flag to indicate that refined vectors must be
   computed.

   Logically Collective

   Input Parameters:
+  ds  - the direct solver context
-  ref - a boolean flag

   Notes:
   Normally the vectors returned in DS_MAT_X are eigenvectors of the
   projected matrix. With this flag activated, DSVectors() will return
   the right singular vector of the smallest singular value of matrix
   \tilde{A}-theta*I, where \tilde{A} is the extended (n+1)xn matrix
   and theta is the Ritz value. This is used in the refined Ritz
   approximation.

   The default is PETSC_FALSE.

   Level: advanced

.seealso: DSVectors(), DSGetRefined()
@*/
PetscErrorCode DSSetRefined(DS ds,PetscBool ref)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveBool(ds,ref,2);
  ds->refined = ref;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetRefined - Gets the refined vectors flag.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  ref - the flag

   Level: advanced

.seealso: DSSetRefined()
@*/
PetscErrorCode DSGetRefined(DS ds,PetscBool *ref)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(ref,2);
  *ref = ds->refined;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetBlockSize - Sets the block size.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  bs - the block size

   Options Database Key:
.  -ds_block_size <bs> - Sets the block size

   Level: intermediate

.seealso: DSGetBlockSize()
@*/
PetscErrorCode DSSetBlockSize(DS ds,PetscInt bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,bs,2);
  PetscCheck(bs>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The block size must be at least one");
  ds->bs = bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSGetBlockSize - Gets the block size.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  bs - block size

   Level: intermediate

.seealso: DSSetBlockSize()
@*/
PetscErrorCode DSGetBlockSize(DS ds,PetscInt *bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(bs,2);
  *bs = ds->bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSSetSlepcSC - Sets the sorting criterion context.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  sc - a pointer to the sorting criterion context

   Note:
   Not available in Fortran.

   Level: developer

.seealso: DSGetSlepcSC(), DSSort()
@*/
PetscErrorCode DSSetSlepcSC(DS ds,SlepcSC sc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(sc,2);
  if (ds->sc && !ds->scset) PetscCall(PetscFree(ds->sc));
  ds->sc    = sc;
  ds->scset = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSGetSlepcSC - Gets the sorting criterion context.

   Not Collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameter:
.  sc - a pointer to the sorting criterion context

   Note:
   Not available in Fortran.

   Level: developer

.seealso: DSSetSlepcSC(), DSSort()
@*/
PetscErrorCode DSGetSlepcSC(DS ds,SlepcSC *sc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscAssertPointer(sc,2);
  if (!ds->sc) PetscCall(PetscNew(&ds->sc));
  *sc = ds->sc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSSetFromOptions - Sets DS options from the options database.

   Collective

   Input Parameters:
.  ds - the direct solver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso: DSSetOptionsPrefix()
@*/
PetscErrorCode DSSetFromOptions(DS ds)
{
  PetscInt       bs,meth;
  PetscBool      flag;
  DSParallelType pmode;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscCall(DSRegisterAll());
  /* Set default type (we do not allow changing it with -ds_type) */
  if (!((PetscObject)ds)->type_name) PetscCall(DSSetType(ds,DSNHEP));
  PetscObjectOptionsBegin((PetscObject)ds);

    PetscCall(PetscOptionsInt("-ds_block_size","Block size for the dense system solver","DSSetBlockSize",ds->bs,&bs,&flag));
    if (flag) PetscCall(DSSetBlockSize(ds,bs));

    PetscCall(PetscOptionsInt("-ds_method","Method to be used for the dense system","DSSetMethod",ds->method,&meth,&flag));
    if (flag) PetscCall(DSSetMethod(ds,meth));

    PetscCall(PetscOptionsEnum("-ds_parallel","Operation mode in parallel runs","DSSetParallel",DSParallelTypes,(PetscEnum)ds->pmode,(PetscEnum*)&pmode,&flag));
    if (flag) PetscCall(DSSetParallel(ds,pmode));

    PetscTryTypeMethod(ds,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)ds,PetscOptionsObject));
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSView - Prints the DS data structure.

   Collective

   Input Parameters:
+  ds - the direct solver context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization context with
   PetscViewerASCIIOpen() - output to a specified file.

   Level: beginner

.seealso: DSViewMat()
@*/
PetscErrorCode DSView(DS ds,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;
  PetscMPIInt       size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ds),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(ds,1,viewer,2);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)ds,viewer));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ds),&size));
    if (size>1) PetscCall(PetscViewerASCIIPrintf(viewer,"  parallel operation mode: %s\n",DSParallelTypes[ds->pmode]));
    if (format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  current state: %s\n",DSStateTypes[ds->state]));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  dimensions: ld=%" PetscInt_FMT ", n=%" PetscInt_FMT ", l=%" PetscInt_FMT ", k=%" PetscInt_FMT,ds->ld,ds->n,ds->l,ds->k));
      if (ds->state==DS_STATE_TRUNCATED) PetscCall(PetscViewerASCIIPrintf(viewer,", t=%" PetscInt_FMT "\n",ds->t));
      else PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  flags:%s%s%s\n",ds->compact?" compact":"",ds->extrarow?" extrarow":"",ds->refined?" refined":""));
    }
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(ds,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSViewFromOptions - View from options

   Collective

   Input Parameters:
+  ds   - the direct solver context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: DSView(), DSCreate()
@*/
PetscErrorCode DSViewFromOptions(DS ds,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)ds,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSAllocate - Allocates memory for internal storage or matrices in DS.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  ld - leading dimension (maximum allowed dimension for the matrices, including
        the extra row if present)

   Note:
   If the leading dimension is different from a previously set value, then
   all matrices are destroyed with DSReset().

   Level: intermediate

.seealso: DSGetLeadingDimension(), DSSetDimensions(), DSSetExtraRow(), DSReset(), DSReallocate()
@*/
PetscErrorCode DSAllocate(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,ld,2);
  PetscValidType(ds,1);
  PetscCheck(ld>0,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Leading dimension should be at least one");
  if (ld!=ds->ld) {
    PetscCall(PetscInfo(ds,"Allocating memory with leading dimension=%" PetscInt_FMT "\n",ld));
    PetscCall(DSReset(ds));
    ds->ld = ld;
    PetscUseTypeMethod(ds,allocate,ld);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSReallocate - Reallocates memory for internal storage or matrices in DS,
   keeping the previously set data.

   Logically Collective

   Input Parameters:
+  ds - the direct solver context
-  ld - new leading dimension

   Notes:
   The new leading dimension must be larger than the previous one. The relevant
   data previously set is copied over to the new data structures.

   This operation is not available in all DS types.

   Level: developer

.seealso: DSAllocate()
@*/
PetscErrorCode DSReallocate(DS ds,PetscInt ld)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,ld,2);
  PetscValidType(ds,1);
  PetscCheck(ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ORDER,"DSAllocate() must be called first");
  PetscCheck(ld>ds->ld,PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"New leading dimension %" PetscInt_FMT " must be larger than the previous one %" PetscInt_FMT,ld,ds->ld);
  PetscCall(PetscInfo(ds,"Reallocating memory with new leading dimension=%" PetscInt_FMT "\n",ld));
  PetscUseTypeMethod(ds,reallocate,ld);
  ds->ld = ld;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSReset - Resets the DS context to the initial state.

   Collective

   Input Parameter:
.  ds - the direct solver context

   Note:
   All data structures with size depending on the leading dimension
   of DSAllocate() are released.

   Level: advanced

.seealso: DSDestroy(), DSAllocate()
@*/
PetscErrorCode DSReset(DS ds)
{
  PetscInt       i;

  PetscFunctionBegin;
  if (ds) PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  if (!ds) PetscFunctionReturn(PETSC_SUCCESS);
  ds->state    = DS_STATE_RAW;
  ds->ld       = 0;
  ds->l        = 0;
  ds->n        = 0;
  ds->k        = 0;
  for (i=0;i<DS_NUM_MAT;i++) PetscCall(MatDestroy(&ds->omat[i]));
  PetscCall(PetscFree(ds->perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   DSDestroy - Destroys DS context that was created with DSCreate().

   Collective

   Input Parameter:
.  ds - the direct solver context

   Level: beginner

.seealso: DSCreate()
@*/
PetscErrorCode DSDestroy(DS *ds)
{
  PetscFunctionBegin;
  if (!*ds) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*ds,DS_CLASSID,1);
  if (--((PetscObject)*ds)->refct > 0) { *ds = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(DSReset(*ds));
  PetscTryTypeMethod(*ds,destroy);
  PetscCall(PetscFree((*ds)->work));
  PetscCall(PetscFree((*ds)->rwork));
  PetscCall(PetscFree((*ds)->iwork));
  if (!(*ds)->scset) PetscCall(PetscFree((*ds)->sc));
  PetscCall(PetscHeaderDestroy(ds));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   DSRegister - Adds a direct solver to the DS package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined DS
-  function - routine to create context

   Note:
   DSRegister() may be called multiple times to add several user-defined
   direct solvers.

   Level: advanced

.seealso: DSRegisterAll()
@*/
PetscErrorCode DSRegister(const char *name,PetscErrorCode (*function)(DS))
{
  PetscFunctionBegin;
  PetscCall(DSInitializePackage());
  PetscCall(PetscFunctionListAdd(&DSList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode DSCreate_HEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_NHEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_GHEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_GHIEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_GNHEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_NHEPTS(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_SVD(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_HSVD(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_GSVD(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_PEP(DS);
SLEPC_EXTERN PetscErrorCode DSCreate_NEP(DS);

/*@C
   DSRegisterAll - Registers all of the direct solvers in the DS package.

   Not Collective

   Level: advanced

.seealso: DSRegister()
@*/
PetscErrorCode DSRegisterAll(void)
{
  PetscFunctionBegin;
  if (DSRegisterAllCalled) PetscFunctionReturn(PETSC_SUCCESS);
  DSRegisterAllCalled = PETSC_TRUE;
  PetscCall(DSRegister(DSHEP,DSCreate_HEP));
  PetscCall(DSRegister(DSNHEP,DSCreate_NHEP));
  PetscCall(DSRegister(DSGHEP,DSCreate_GHEP));
  PetscCall(DSRegister(DSGHIEP,DSCreate_GHIEP));
  PetscCall(DSRegister(DSGNHEP,DSCreate_GNHEP));
  PetscCall(DSRegister(DSNHEPTS,DSCreate_NHEPTS));
  PetscCall(DSRegister(DSSVD,DSCreate_SVD));
  PetscCall(DSRegister(DSHSVD,DSCreate_HSVD));
  PetscCall(DSRegister(DSGSVD,DSCreate_GSVD));
  PetscCall(DSRegister(DSPEP,DSCreate_PEP));
  PetscCall(DSRegister(DSNEP,DSCreate_NEP));
  PetscFunctionReturn(PETSC_SUCCESS);
}
