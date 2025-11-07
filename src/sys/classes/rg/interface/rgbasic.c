/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic RG routines
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/

PetscFunctionList RGList = NULL;
PetscBool         RGRegisterAllCalled = PETSC_FALSE;
PetscClassId      RG_CLASSID = 0;
static PetscBool  RGPackageInitialized = PETSC_FALSE;

/*@C
   RGFinalizePackage - This function destroys everything in the SLEPc interface
   to the `RG` package. It is called from `SlepcFinalize()`.

   Level: developer

.seealso: [](sec:rg), `SlepcFinalize()`, `BVInitializePackage()`
@*/
PetscErrorCode RGFinalizePackage(void)
{
  PetscFunctionBegin;
  PetscCall(PetscFunctionListDestroy(&RGList));
  RGPackageInitialized = PETSC_FALSE;
  RGRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   RGInitializePackage - This function initializes everything in the `RG` package.
   It is called from `PetscDLLibraryRegister_slepc()` when using dynamic libraries, and
   on the first call to `RGCreate()` when using shared or static libraries.

   Note:
   This function never needs to be called by SLEPc users.

   Level: developer

.seealso: [](sec:rg), `RG`, `SlepcInitialize()`, `RGFinalizePackage()`
@*/
PetscErrorCode RGInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];

  PetscFunctionBegin;
  if (RGPackageInitialized) PetscFunctionReturn(PETSC_SUCCESS);
  RGPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  PetscCall(PetscClassIdRegister("Region",&RG_CLASSID));
  /* Register Constructors */
  PetscCall(RGRegisterAll());
  /* Process Info */
  classids[0] = RG_CLASSID;
  PetscCall(PetscInfoProcessClass("rg",1,&classids[0]));
  /* Process summary exclusions */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt));
  if (opt) {
    PetscCall(PetscStrInList("rg",logList,',',&pkg));
    if (pkg) PetscCall(PetscLogEventDeactivateClass(RG_CLASSID));
  }
  /* Register package finalizer */
  PetscCall(PetscRegisterFinalize(RGFinalizePackage));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGCreate - Creates an `RG` context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newrg - location to put the `RG` context

   Level: beginner

.seealso: [](sec:rg), `RG`, `RGDestroy()`
@*/
PetscErrorCode RGCreate(MPI_Comm comm,RG *newrg)
{
  RG             rg;

  PetscFunctionBegin;
  PetscAssertPointer(newrg,2);
  PetscCall(RGInitializePackage());
  PetscCall(SlepcHeaderCreate(rg,RG_CLASSID,"RG","Region","RG",comm,RGDestroy,RGView));
  rg->complement = PETSC_FALSE;
  rg->sfactor    = 1.0;
  rg->osfactor   = 0.0;
  rg->data       = NULL;

  *newrg = rg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGSetOptionsPrefix - Sets the prefix used for searching for all
   `RG` options in the database.

   Logically Collective

   Input Parameters:
+  rg     - the region context
-  prefix - the prefix string to prepend to all `RG` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: [](sec:rg), `RGAppendOptionsPrefix()`
@*/
PetscErrorCode RGSetOptionsPrefix(RG rg,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)rg,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGAppendOptionsPrefix - Appends to the prefix used for searching for all
   `RG` options in the database.

   Logically Collective

   Input Parameters:
+  rg     - the region context
-  prefix - the prefix string to prepend to all `RG` option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: [](sec:rg), `RGSetOptionsPrefix()`
@*/
PetscErrorCode RGAppendOptionsPrefix(RG rg,const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscCall(PetscObjectAppendOptionsPrefix((PetscObject)rg,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGGetOptionsPrefix - Gets the prefix used for searching for all
   `RG` options in the database.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  prefix - pointer to the prefix string used is returned

   Level: advanced

.seealso: [](sec:rg), `RGSetOptionsPrefix()`, `RGAppendOptionsPrefix()`
@*/
PetscErrorCode RGGetOptionsPrefix(RG rg,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscAssertPointer(prefix,2);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)rg,prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGSetType - Selects the type for the `RG` object.

   Logically Collective

   Input Parameters:
+  rg   - the region context
-  type - a known type

   Options Database Key:
.  -rg_type \<type\> - sets `RG` type

   Level: beginner

.seealso: [](sec:rg), `RGGetType()`
@*/
PetscErrorCode RGSetType(RG rg,RGType type)
{
  PetscErrorCode (*r)(RG);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscAssertPointer(type,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)rg,type,&match));
  if (match) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscFunctionListFind(RGList,type,&r));
  PetscCheck(r,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested RG type %s",type);

  PetscTryTypeMethod(rg,destroy);
  PetscCall(PetscMemzero(rg->ops,sizeof(struct _RGOps)));

  PetscCall(PetscObjectChangeTypeName((PetscObject)rg,type));
  PetscCall((*r)(rg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGGetType - Gets the `RG` type name (as a string) from the `RG` context.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  type - name of the region

   Level: beginner

.seealso: [](sec:rg), `RGSetType()`
@*/
PetscErrorCode RGGetType(RG rg,RGType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscAssertPointer(type,2);
  *type = ((PetscObject)rg)->type_name;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGSetFromOptions - Sets `RG` options from the options database.

   Collective

   Input Parameter:
.  rg - the region context

   Note:
   To see all options, run your program with the `-help` option.

   Level: beginner

.seealso: [](sec:rg), `RGSetOptionsPrefix()`
@*/
PetscErrorCode RGSetFromOptions(RG rg)
{
  char           type[256];
  PetscBool      flg;
  PetscReal      sfactor;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscCall(RGRegisterAll());
  PetscObjectOptionsBegin((PetscObject)rg);
    PetscCall(PetscOptionsFList("-rg_type","Region type","RGSetType",RGList,(char*)(((PetscObject)rg)->type_name?((PetscObject)rg)->type_name:RGINTERVAL),type,sizeof(type),&flg));
    if (flg) PetscCall(RGSetType(rg,type));
    else if (!((PetscObject)rg)->type_name) PetscCall(RGSetType(rg,RGINTERVAL));

    PetscCall(PetscOptionsBool("-rg_complement","Whether region is complemented or not","RGSetComplement",rg->complement,&rg->complement,NULL));

    PetscCall(PetscOptionsReal("-rg_scale","Scaling factor","RGSetScale",1.0,&sfactor,&flg));
    if (flg) PetscCall(RGSetScale(rg,sfactor));

    PetscTryTypeMethod(rg,setfromoptions,PetscOptionsObject);
    PetscCall(PetscObjectProcessOptionsHandlers((PetscObject)rg,PetscOptionsObject));
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGView - Prints the `RG` data structure.

   Collective

   Input Parameters:
+  rg - the region context
-  viewer - optional visualization context

   Notes:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Use `RGViewFromOptions()` to allow the user to select many different `PetscViewerType`
   and formats from the options database.

   Level: beginner

.seealso: [](sec:rg), `RGCreate()`, `RGViewFromOptions()`
@*/
PetscErrorCode RGView(RG rg,PetscViewer viewer)
{
  PetscBool      isdraw,isascii;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)rg),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(rg,1,viewer,2);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)rg,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(rg,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    if (rg->complement) PetscCall(PetscViewerASCIIPrintf(viewer,"  selected region is the complement of the specified one\n"));
    if (rg->sfactor!=1.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  scaling factor = %g\n",(double)rg->sfactor));
  } else if (isdraw) PetscTryTypeMethod(rg,view,viewer);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGViewFromOptions - View (print) an `RG` object based on values in the options database.

   Collective

   Input Parameters:
+  rg   - the region context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](sec:rg), `RGView()`, `RGCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode RGViewFromOptions(RG rg,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)rg,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGIsTrivial - Whether it is the trivial region (whole complex plane).

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  trivial - true if the region is equal to the whole complex plane, e.g.,
             an interval region with all four endpoints unbounded or an
             ellipse with infinite radius.

   Level: beginner

.seealso: [](sec:rg), `RGCheckInside()`
@*/
PetscErrorCode RGIsTrivial(RG rg,PetscBool *trivial)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscAssertPointer(trivial,2);
  *trivial = PETSC_FALSE;
  PetscTryTypeMethod(rg,istrivial,trivial);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGCheckInside - Determines if a set of given points are inside the region or not.

   Not Collective

   Input Parameters:
+  rg - the region context
.  n  - number of points to check
.  ar - array of real parts
-  ai - array of imaginary parts

   Output Parameter:
.  inside - array of results (1=inside, 0=on the contour, -1=outside)

   Notes:
   The point `a` is expressed as a couple of `PetscScalar` variables `ar`, `ai`.
   If built with complex scalars, the point is supposed to be stored in `ar`,
   otherwise `ar`, `ai` contain the real and imaginary parts, respectively.

   If a scaling factor was set, the points are scaled before checking.

   Level: intermediate

.seealso: [](sec:rg), `RGSetScale()`, `RGSetComplement()`
@*/
PetscErrorCode RGCheckInside(RG rg,PetscInt n,PetscScalar ar[],PetscScalar ai[],PetscInt inside[])
{
  PetscReal      px,py;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscAssertPointer(ar,3);
#if !defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(ai,4);
#endif
  PetscAssertPointer(inside,5);

  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    px = PetscRealPart(ar[i]);
    py = PetscImaginaryPart(ar[i]);
#else
    px = ar[i];
    py = ai[i];
#endif
    if (PetscUnlikely(rg->sfactor != 1.0)) {
      px /= rg->sfactor;
      py /= rg->sfactor;
    }
    PetscUseTypeMethod(rg,checkinside,px,py,inside+i);
    if (PetscUnlikely(rg->complement)) inside[i] = -inside[i];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGIsAxisymmetric - Determines if the region is symmetric with respect
   to the real or imaginary axis.

   Not Collective

   Input Parameters:
+  rg       - the region context
-  vertical - `PETSC_TRUE` if symmetry must be checked against the vertical axis

   Output Parameter:
.  symm - `PETSC_TRUE` if the region is axisymmetric

   Note:
   If the vertical argument is true, symmetry is checked with respect to
   the vertical axis, otherwise with respect to the horizontal axis.

   Level: intermediate

.seealso: [](sec:rg), `RGCanUseConjugates()`
@*/
PetscErrorCode RGIsAxisymmetric(RG rg,PetscBool vertical,PetscBool *symm)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscAssertPointer(symm,3);
  *symm = PETSC_FALSE;
  PetscTryTypeMethod(rg,isaxisymmetric,vertical,symm);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGCanUseConjugates - Used in contour integral methods to determine whether
   half of integration points can be avoided (use their conjugates).

   Not Collective

   Input Parameters:
+  rg       - the region context
-  realmats - `PETSC_TRUE` if the problem matrices are real

   Output Parameter:
.  useconj  - whether it is possible to use conjugates

   Notes:
   If some integration points are the conjugates of other points, then the
   associated computational cost can be saved. This depends on the problem
   matrices being real and also the region being symmetric with respect to
   the horizontal axis. The result is false if using real arithmetic or
   in the case of a flat region (height equal to zero).

   Level: developer

.seealso: [](sec:rg), `RGIsAxisymmetric()`
@*/
PetscErrorCode RGCanUseConjugates(RG rg,PetscBool realmats,PetscBool *useconj)
{
#if defined(PETSC_USE_COMPLEX)
  PetscReal      c,d;
  PetscBool      isaxisymm;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscAssertPointer(useconj,3);
  *useconj = PETSC_FALSE;
#if defined(PETSC_USE_COMPLEX)
  if (realmats) {
    PetscCall(RGIsAxisymmetric(rg,PETSC_FALSE,&isaxisymm));
    if (isaxisymm) {
      PetscCall(RGComputeBoundingBox(rg,NULL,NULL,&c,&d));
      if (c!=d) *useconj = PETSC_TRUE;
    }
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGComputeContour - Computes the coordinates of several points lying on the
   contour of the region.

   Not Collective

   Input Parameters:
+  rg - the region context
-  n  - number of points to compute

   Output Parameters:
+  cr - location to store real parts
-  ci - location to store imaginary parts

   Note:
   In real scalars, either `cr` or `ci` can be `NULL` (but not both). In complex
   scalars, the coordinates are stored in `cr`, which cannot be `NULL` (`ci` is
   not referenced).

   Level: intermediate

.seealso: [](sec:rg), `RGComputeBoundingBox()`, `RGSetScale()`
@*/
PetscErrorCode RGComputeContour(RG rg,PetscInt n,PetscScalar cr[],PetscScalar ci[])
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
#if defined(PETSC_USE_COMPLEX)
  PetscAssertPointer(cr,3);
#else
  PetscCheck(cr || ci,PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"cr and ci cannot be NULL simultaneously");
#endif
  PetscCheck(!rg->complement,PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Cannot compute contour of region with complement flag set");
  PetscUseTypeMethod(rg,computecontour,n,cr,ci);
  for (i=0;i<n;i++) {
    if (cr) cr[i] *= rg->sfactor;
    if (ci) ci[i] *= rg->sfactor;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGComputeBoundingBox - Determines the endpoints of a rectangle in the complex plane that
   contains the region.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameters:
+  a - left endpoint of the bounding box in the real axis
.  b - right endpoint of the bounding box in the real axis
.  c - bottom endpoint of the bounding box in the imaginary axis
-  d - top endpoint of the bounding box in the imaginary axis

   Note:
   The bounding box is defined as $[a,b]\times[c,d]$. In regions that are not bounded (e.g., an
   open interval) or with the complement flag set, it makes no sense to compute a bounding
   box, so the return values are infinite.

   Level: intermediate

.seealso: [](sec:rg), `RGComputeContour()`, `RGSetScale()`, `RGSetComplement()`
@*/
PetscErrorCode RGComputeBoundingBox(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);

  if (rg->complement) {  /* cannot compute bounding box */
    if (a) *a = -PETSC_MAX_REAL;
    if (b) *b =  PETSC_MAX_REAL;
    if (c) *c = -PETSC_MAX_REAL;
    if (d) *d =  PETSC_MAX_REAL;
  } else {
    PetscUseTypeMethod(rg,computebbox,a,b,c,d);
    if (a && *a!=-PETSC_MAX_REAL) *a *= rg->sfactor;
    if (b && *b!= PETSC_MAX_REAL) *b *= rg->sfactor;
    if (c && *c!=-PETSC_MAX_REAL) *c *= rg->sfactor;
    if (d && *d!= PETSC_MAX_REAL) *d *= rg->sfactor;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGComputeQuadrature - Computes the values of the parameters used in a
   quadrature rule for a contour integral around the boundary of the region.

   Not Collective

   Input Parameters:
+  rg   - the region context
.  quad - the type of quadrature
-  n    - number of quadrature points to compute

   Output Parameters:
+  z  - quadrature points
.  zn - normalized quadrature points
-  w  - quadrature weights

   Notes:
   In complex scalars, the values returned in `z` are often the same as those
   computed by `RGComputeContour()`, but this is not the case in real scalars
   where all output arguments are real.

   The computed values change for different quadrature rules.

   Level: advanced

.seealso: [](sec:rg), `RGComputeContour()`
@*/
PetscErrorCode RGComputeQuadrature(RG rg,RGQuadRule quad,PetscInt n,PetscScalar z[],PetscScalar zn[],PetscScalar w[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidType(rg,1);
  PetscAssertPointer(z,4);
  PetscAssertPointer(zn,5);
  PetscAssertPointer(w,6);

  PetscCall(RGComputeContour(rg,n,z,NULL));
  PetscUseTypeMethod(rg,computequadrature,quad,n,z,zn,w);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGSetComplement - Sets a flag to indicate that the region is the complement
   of the specified one.

   Logically Collective

   Input Parameters:
+  rg  - the region context
-  flg - the boolean flag

   Options Database Key:
.  -rg_complement - activate/deactivate complementing of the region

   Level: intermediate

.seealso: [](sec:rg), `RGGetComplement()`
@*/
PetscErrorCode RGSetComplement(RG rg,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveBool(rg,flg,2);
  rg->complement = flg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGGetComplement - Gets a flag that indicates whether the region
   is complemented or not.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  flg - the flag

   Level: intermediate

.seealso: [](sec:rg), `RGSetComplement()`
@*/
PetscErrorCode RGGetComplement(RG rg,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscAssertPointer(flg,2);
  *flg = rg->complement;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGSetScale - Sets the scaling factor to be used when checking that a
   point is inside the region and when computing the contour.

   Logically Collective

   Input Parameters:
+  rg      - the region context
-  sfactor - the scaling factor

   Options Database Key:
.  -rg_scale \<sfactor\> - sets the scaling factor

   Level: intermediate

.seealso: [](sec:rg), `RGGetScale()`, `RGCheckInside()`, `RGComputeContour()`
@*/
PetscErrorCode RGSetScale(RG rg,PetscReal sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveReal(rg,sfactor,2);
  if (sfactor == (PetscReal)PETSC_DEFAULT || sfactor == (PetscReal)PETSC_DECIDE) sfactor = 1.0;
  PetscCheck(sfactor>0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of scaling factor. Must be > 0");
  rg->sfactor = sfactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGGetScale - Gets the scaling factor.

   Not Collective

   Input Parameter:
.  rg - the region context

   Output Parameter:
.  sfactor - the scaling factor

   Level: intermediate

.seealso: [](sec:rg), `RGSetScale()`
@*/
PetscErrorCode RGGetScale(RG rg,PetscReal *sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscAssertPointer(sfactor,2);
  *sfactor = rg->sfactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGPushScale - Sets an additional scaling factor, that will multiply the
   user-defined scaling factor.

   Logically Collective

   Input Parameters:
+  rg      - the region context
-  sfactor - the scaling factor

   Notes:
   The current implementation does not allow pushing several scaling factors.

   This is intended for internal use, for instance in polynomial eigensolvers
   that use parameter scaling.

   Level: developer

.seealso: [](sec:rg), `RGPopScale()`, `RGSetScale()`
@*/
PetscErrorCode RGPushScale(RG rg,PetscReal sfactor)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscValidLogicalCollectiveReal(rg,sfactor,2);
  PetscCheck(sfactor>0.0,PetscObjectComm((PetscObject)rg),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of scaling factor. Must be > 0");
  PetscCheck(!rg->osfactor,PetscObjectComm((PetscObject)rg),PETSC_ERR_SUP,"Current implementation does not allow pushing several scaling factors");
  rg->osfactor = rg->sfactor;
  rg->sfactor *= sfactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGPopScale - Pops the scaling factor set with `RGPushScale()`.

   Logically Collective

   Input Parameter:
.  rg - the region context

   Level: developer

.seealso: [](sec:rg), `RGPushScale()`
@*/
PetscErrorCode RGPopScale(RG rg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(rg,RG_CLASSID,1);
  PetscCheck(rg->osfactor,PetscObjectComm((PetscObject)rg),PETSC_ERR_ORDER,"Must call RGPushScale first");
  rg->sfactor  = rg->osfactor;
  rg->osfactor = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   RGDestroy - Destroys `RG` context that was created with `RGCreate()`.

   Collective

   Input Parameter:
.  rg - the region context

   Level: beginner

.seealso: [](sec:rg), `RGCreate()`
@*/
PetscErrorCode RGDestroy(RG *rg)
{
  PetscFunctionBegin;
  if (!*rg) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific(*rg,RG_CLASSID,1);
  if (--((PetscObject)*rg)->refct > 0) { *rg = NULL; PetscFunctionReturn(PETSC_SUCCESS); }
  PetscTryTypeMethod(*rg,destroy);
  PetscCall(PetscHeaderDestroy(rg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   RGRegister - Adds a region to the `RG` package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined RG
-  function - routine to create context

   Note:
   `RGRegister()` may be called multiple times to add several user-defined regions.

   Level: advanced

.seealso: [](sec:rg), `RGRegisterAll()`
@*/
PetscErrorCode RGRegister(const char *name,PetscErrorCode (*function)(RG))
{
  PetscFunctionBegin;
  PetscCall(RGInitializePackage());
  PetscCall(PetscFunctionListAdd(&RGList,name,function));
  PetscFunctionReturn(PETSC_SUCCESS);
}
