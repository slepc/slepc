/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to monitors
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode PEPMonitor(PEP pep,PetscInt it,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest)
{
  PetscInt       i,n = pep->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) PetscCall((*pep->monitor[i])(pep,it,nconv,eigr,eigi,errest,nest,pep->monitorcontext[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested eigenpair.

   Logically Collective

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL),
             see PetscCtxDestroyFn for the calling sequence

   Calling sequence of monitor:
$  PetscErrorCode monitor(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx)
+  pep    - polynomial eigensolver context obtained from PEPCreate()
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by PEPMonitorSet()

   Options Database Keys:
+    -pep_monitor        - print only the first error estimate
.    -pep_monitor_all    - print error estimates at each iteration
.    -pep_monitor_conv   - print the eigenvalue approximations only when
      convergence has been reached
.    -pep_monitor draw::draw_lg - sets line graph monitor for the first unconverged
      approximate eigenvalue
.    -pep_monitor_all draw::draw_lg - sets line graph monitor for all unconverged
      approximate eigenvalues
.    -pep_monitor_conv draw::draw_lg - sets line graph monitor for convergence history
-    -pep_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to PEPMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   PEPMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: PEPMonitorFirst(), PEPMonitorAll(), PEPMonitorCancel()
@*/
PetscErrorCode PEPMonitorSet(PEP pep,PetscErrorCode (*monitor)(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscCtxDestroyFn *monitordestroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscCheck(pep->numbermonitors<MAXPEPMONITORS,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Too many PEP monitors set");
  pep->monitor[pep->numbermonitors]           = monitor;
  pep->monitorcontext[pep->numbermonitors]    = (void*)mctx;
  pep->monitordestroy[pep->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPMonitorCancel - Clears all monitors for a PEP object.

   Logically Collective

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Options Database Key:
.    -pep_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to PEPMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorCancel(PEP pep)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  for (i=0; i<pep->numbermonitors; i++) {
    if (pep->monitordestroy[i]) PetscCall((*pep->monitordestroy[i])(&pep->monitorcontext[i]));
  }
  pep->numbermonitors = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPGetMonitorContext - Gets the monitor context, as set by
   PEPMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: PEPMonitorSet(), PEPDefaultMonitor()
@*/
PetscErrorCode PEPGetMonitorContext(PEP pep,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  *(void**)ctx = pep->monitorcontext[0];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Helper function to compute eigenvalue that must be viewed in monitor
 */
static PetscErrorCode PEPMonitorGetTrueEig(PEP pep,PetscScalar *er,PetscScalar *ei)
{
  PetscBool      flg,sinv;

  PetscFunctionBegin;
  PetscCall(STBackTransform(pep->st,1,er,ei));
  if (pep->sfactor!=1.0) {
    PetscCall(STGetTransform(pep->st,&flg));
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
    if (flg && sinv) {
      *er /= pep->sfactor; *ei /= pep->sfactor;
    } else {
      *er *= pep->sfactor; *ei *= pep->sfactor;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor - activates PEPMonitorFirst()

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorAll(), PEPMonitorConverged()
@*/
PetscErrorCode PEPMonitorFirst(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)pep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)pep)->prefix));
  if (nconv<nest) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " PEP nconv=%" PetscInt_FMT " first unconverged value (error)",its,nconv));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    er = eigr[nconv]; ei = eigi[nconv];
    PetscCall(PEPMonitorGetTrueEig(pep,&er,&ei));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er)));
#else
    PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)er));
    if (ei!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei));
#endif
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[nconv]));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorAll - Print the current approximate values and
   error estimates at each iteration of the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor_all - activates PEPMonitorAll()

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorFirst(), PEPMonitorConverged()
@*/
PetscErrorCode PEPMonitorAll(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel));
  if (its==1 && ((PetscObject)pep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)pep)->prefix));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " PEP nconv=%" PetscInt_FMT " Values (Errors)",its,nconv));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
  for (i=0;i<nest;i++) {
    er = eigr[i]; ei = eigi[i];
    PetscCall(PEPMonitorGetTrueEig(pep,&er,&ei));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er)));
#else
    PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)er));
    if (ei!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei));
#endif
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)",(double)errest[i]));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
  PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor_conv - activates PEPMonitorConverged()

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorFirst(), PEPMonitorAll()
@*/
PetscErrorCode PEPMonitorConverged(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;
  SlepcConvMon   ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  ctx = (SlepcConvMon)vf->data;
  if (its==1 && ((PetscObject)pep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)pep)->prefix));
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel));
    for (i=ctx->oldnconv;i<nconv;i++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " PEP converged value (error) #%" PetscInt_FMT,its,i));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      er = eigr[i]; ei = eigi[i];
      PetscCall(PEPMonitorGetTrueEig(pep,&er,&ei));
#if defined(PETSC_USE_COMPLEX)
      PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er)));
#else
      PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)er));
      if (ei!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei));
#endif
      PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[i]));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPMonitorConvergedCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPMonitorConvergedDestroy(PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  if (!*vf) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree((*vf)->data));
  PetscCall(PetscViewerDestroy(&(*vf)->viewer));
  PetscCall(PetscFree(*vf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorFirstDrawLG - Plots the error estimate of the first unconverged
   approximation at each iteration of the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor draw::draw_lg - activates PEPMonitorFirstDrawLG()

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorFirstDrawLG(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscReal      x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(pep->tol)-2,0.0));
  }
  if (nconv<nest) {
    x = (PetscReal)its;
    if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
    else y = 0.0;
    PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
    if (its <= 20 || !(its % 5) || pep->reason) {
      PetscCall(PetscDrawLGDraw(lg));
      PetscCall(PetscDrawLGSave(lg));
    }
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorFirstDrawLGCreate - Creates the plotter for the first error estimate.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorFirstDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"First Error Estimate","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorAllDrawLG - Plots the error estimate of all unconverged
   approximations at each iteration of the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor_all draw::draw_lg - activates PEPMonitorAllDrawLG()

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorAllDrawLG(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscInt       i,n = PetscMin(pep->nev,255);
  PetscReal      *x,*y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,n));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(pep->tol)-2,0.0));
  }
  PetscCall(PetscMalloc2(n,&x,n,&y));
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  PetscCall(PetscDrawLGAddPoint(lg,x,y));
  if (its <= 20 || !(its % 5) || pep->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscFree2(x,y));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorAllDrawLGCreate - Creates the plotter for all the error estimates.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorAllDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"All Error Estimates","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorConvergedDrawLG - Plots the number of converged eigenvalues
   at each iteration of the polynomial eigensolver.

   Collective

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -pep_monitor_conv draw::draw_lg - activates PEPMonitorConvergedDrawLG()

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorConvergedDrawLG(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer      viewer = vf->viewer;
  PetscDrawLG      lg;
  PetscReal        x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,0,pep->nev));
  }
  x = (PetscReal)its;
  y = (PetscReal)pep->nconv;
  PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
  if (its <= 20 || !(its % 5) || pep->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   PEPMonitorConvergedDrawLGCreate - Creates the plotter for the convergence history.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorConvergedDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"Convergence History","Number Converged",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}
