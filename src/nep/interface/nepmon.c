/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to monitors
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode NEPMonitor(NEP nep,PetscInt it,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest)
{
  PetscInt       i,n = nep->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) PetscCall((*nep->monitor[i])(nep,it,nconv,eigr,eigi,errest,nest,nep->monitorcontext[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested eigenpair.

   Logically Collective

   Input Parameters:
+  nep     - eigensolver context obtained from NEPCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL),
             see PetscCtxDestroyFn for the calling sequence

   Calling sequence of monitor:
$  PetscErrorCode monitor(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx)
+  nep    - nonlinear eigensolver context obtained from NEPCreate()
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates for each eigenpair
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by NEPMonitorSet()

   Options Database Keys:
+    -nep_monitor        - print only the first error estimate
.    -nep_monitor_all    - print error estimates at each iteration
.    -nep_monitor_conv   - print the eigenvalue approximations only when
      convergence has been reached
.    -nep_monitor draw::draw_lg - sets line graph monitor for the first unconverged
      approximate eigenvalue
.    -nep_monitor_all draw::draw_lg - sets line graph monitor for all unconverged
      approximate eigenvalues
.    -nep_monitor_conv draw::draw_lg - sets line graph monitor for convergence history
-    -nep_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to NEPMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   NEPMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: NEPMonitorFirst(), NEPMonitorAll(), NEPMonitorCancel()
@*/
PetscErrorCode NEPMonitorSet(NEP nep,PetscErrorCode (*monitor)(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscCtxDestroyFn *monitordestroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscCheck(nep->numbermonitors<MAXNEPMONITORS,PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Too many NEP monitors set");
  nep->monitor[nep->numbermonitors]           = monitor;
  nep->monitorcontext[nep->numbermonitors]    = (void*)mctx;
  nep->monitordestroy[nep->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPMonitorCancel - Clears all monitors for a NEP object.

   Logically Collective

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Options Database Key:
.    -nep_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to NEPMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorCancel(NEP nep)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  for (i=0; i<nep->numbermonitors; i++) {
    if (nep->monitordestroy[i]) PetscCall((*nep->monitordestroy[i])(&nep->monitorcontext[i]));
  }
  nep->numbermonitors = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPGetMonitorContext - Gets the monitor context, as set by
   NEPMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: NEPMonitorSet(), NEPDefaultMonitor()
@*/
PetscErrorCode NEPGetMonitorContext(NEP nep,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  *(void**)ctx = nep->monitorcontext[0];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor - activates NEPMonitorFirst()

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorAll(), NEPMonitorConverged()
@*/
PetscErrorCode NEPMonitorFirst(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)nep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)nep)->prefix));
  if (nconv<nest) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " NEP nconv=%" PetscInt_FMT " first unconverged value (error)",its,nconv));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[nconv]),(double)PetscImaginaryPart(eigr[nconv])));
#else
    PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)eigr[nconv]));
    if (eigi[nconv]!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[nconv]));
#endif
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[nconv]));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorAll - Print the current approximate values and
   error estimates at each iteration of the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor_all - activates NEPMonitorAll()

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorFirst(), NEPMonitorConverged()
@*/
PetscErrorCode NEPMonitorAll(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel));
  if (its==1 && ((PetscObject)nep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)nep)->prefix));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " NEP nconv=%" PetscInt_FMT " Values (Errors)",its,nconv));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
  for (i=0;i<nest;i++) {
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[i]),(double)PetscImaginaryPart(eigr[i])));
#else
    PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)eigr[i]));
    if (eigi[i]!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[i]));
#endif
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)",(double)errest[i]));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
  PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor_conv - activates NEPMonitorConverged()

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorFirst(), NEPMonitorAll()
@*/
PetscErrorCode NEPMonitorConverged(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscViewer    viewer = vf->viewer;
  SlepcConvMon   ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  ctx = (SlepcConvMon)vf->data;
  if (its==1 && ((PetscObject)nep)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)nep)->prefix));
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel));
    for (i=ctx->oldnconv;i<nconv;i++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " NEP converged value (error) #%" PetscInt_FMT,its,i));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
#if defined(PETSC_USE_COMPLEX)
      PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[i]),(double)PetscImaginaryPart(eigr[i])));
#else
      PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)eigr[i]));
      if (eigi[i]!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[i]));
#endif
      PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[i]));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPMonitorConvergedCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPMonitorConvergedDestroy(PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  if (!*vf) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree((*vf)->data));
  PetscCall(PetscViewerDestroy(&(*vf)->viewer));
  PetscCall(PetscFree(*vf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorFirstDrawLG - Plots the error estimate of the first unconverged
   approximation at each iteration of the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor draw::draw_lg - activates NEPMonitorFirstDrawLG()

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorFirstDrawLG(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscReal      x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(nep->tol)-2,0.0));
  }
  if (nconv<nest) {
    x = (PetscReal)its;
    if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
    else y = 0.0;
    PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
    if (its <= 20 || !(its % 5) || nep->reason) {
      PetscCall(PetscDrawLGDraw(lg));
      PetscCall(PetscDrawLGSave(lg));
    }
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorFirstDrawLGCreate - Creates the plotter for the first error estimate.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorFirstDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"First Error Estimate","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorAllDrawLG - Plots the error estimate of all unconverged
   approximations at each iteration of the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor_all draw::draw_lg - activates NEPMonitorAllDrawLG()

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorAllDrawLG(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscInt       i,n = PetscMin(nep->nev,255);
  PetscReal      *x,*y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,n));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(nep->tol)-2,0.0));
  }
  PetscCall(PetscMalloc2(n,&x,n,&y));
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  PetscCall(PetscDrawLGAddPoint(lg,x,y));
  if (its <= 20 || !(its % 5) || nep->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscFree2(x,y));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorAllDrawLGCreate - Creates the plotter for all the error estimates.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorAllDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"All Error Estimates","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorConvergedDrawLG - Plots the number of converged eigenvalues
   at each iteration of the nonlinear eigensolver.

   Collective

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -nep_monitor_conv draw::draw_lg - activates NEPMonitorConvergedDrawLG()

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorConvergedDrawLG(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer      viewer = vf->viewer;
  PetscDrawLG      lg;
  PetscReal        x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,0,nep->nev));
  }
  x = (PetscReal)its;
  y = (PetscReal)nep->nconv;
  PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
  if (its <= 20 || !(its % 5) || nep->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   NEPMonitorConvergedDrawLGCreate - Creates the plotter for the convergence history.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorConvergedDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
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
