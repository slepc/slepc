/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SVD routines related to monitors
*/

#include <slepc/private/svdimpl.h>   /*I "slepcsvd.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode SVDMonitor(SVD svd,PetscInt it,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest)
{
  PetscInt       i,n = svd->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) PetscCall((*svd->monitor[i])(svd,it,nconv,sigma,errest,nest,svd->monitorcontext[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested singular triplet.

   Logically Collective

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL),
             see PetscCtxDestroyFn for the calling sequence

   Calling sequence of monitor:
$  PetscErrorCode monitor(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,void *mctx)
+  svd    - singular value solver context obtained from SVDCreate()
.  its    - iteration number
.  nconv  - number of converged singular triplets
.  sigma  - singular values
.  errest - relative error estimates for each singular triplet
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by SVDMonitorSet()

   Options Database Keys:
+    -svd_monitor        - print only the first error estimate
.    -svd_monitor_all    - print error estimates at each iteration
.    -svd_monitor_conv   - print the singular value approximations only when
      convergence has been reached
.    -svd_monitor_conditioning - print the condition number when available
.    -svd_monitor draw::draw_lg - sets line graph monitor for the first unconverged
      approximate singular value
.    -svd_monitor_all draw::draw_lg - sets line graph monitor for all unconverged
      approximate singular values
.    -svd_monitor_conv draw::draw_lg - sets line graph monitor for convergence history
-    -svd_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to SVDMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   SVDMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: SVDMonitorFirst(), SVDMonitorAll(), SVDMonitorConditioning(), SVDMonitorCancel()
@*/
PetscErrorCode SVDMonitorSet(SVD svd,PetscErrorCode (*monitor)(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscCtxDestroyFn *monitordestroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscCheck(svd->numbermonitors<MAXSVDMONITORS,PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Too many SVD monitors set");
  svd->monitor[svd->numbermonitors]           = monitor;
  svd->monitorcontext[svd->numbermonitors]    = (void*)mctx;
  svd->monitordestroy[svd->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDMonitorCancel - Clears all monitors for an SVD object.

   Logically Collective

   Input Parameters:
.  svd - singular value solver context obtained from SVDCreate()

   Options Database Key:
.    -svd_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to SVDMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorCancel(SVD svd)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  for (i=0; i<svd->numbermonitors; i++) {
    if (svd->monitordestroy[i]) PetscCall((*svd->monitordestroy[i])(&svd->monitorcontext[i]));
  }
  svd->numbermonitors = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDGetMonitorContext - Gets the monitor context, as set by
   SVDMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDGetMonitorContext(SVD svd,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  *(void**)ctx = svd->monitorcontext[0];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the singular value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor - activates SVDMonitorFirst()

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorAll(), SVDMonitorConditioning(), SVDMonitorConverged()
@*/
PetscErrorCode SVDMonitorFirst(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  if (its==1 && ((PetscObject)svd)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Singular value approximations and residual norms for %s solve.\n",((PetscObject)svd)->prefix));
  if (nconv<nest) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " SVD nconv=%" PetscInt_FMT " first unconverged value (error)",its,nconv));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    PetscCall(PetscViewerASCIIPrintf(viewer," %g (%10.8e)\n",(double)sigma[nconv],(double)errest[nconv]));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorAll - Print the current approximate values and
   error estimates at each iteration of the singular value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor_all - activates SVDMonitorAll()

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorFirst(), SVDMonitorConditioning(), SVDMonitorConverged()
@*/
PetscErrorCode SVDMonitorAll(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel));
  if (its==1 && ((PetscObject)svd)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Singular value approximations and residual norms for %s solve.\n",((PetscObject)svd)->prefix));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " SVD nconv=%" PetscInt_FMT " Values (Errors)",its,nconv));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
  for (i=0;i<nest;i++) PetscCall(PetscViewerASCIIPrintf(viewer," %g (%10.8e)",(double)sigma[i],(double)errest[i]));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
  PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor_conv - activates SVDMonitorConverged()

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorFirst(), SVDMonitorConditioning(), SVDMonitorAll()
@*/
PetscErrorCode SVDMonitorConverged(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscViewer    viewer = vf->viewer;
  SlepcConvMon   ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  ctx = (SlepcConvMon)vf->data;
  if (its==1 && ((PetscObject)svd)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)svd)->prefix));
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel));
    for (i=ctx->oldnconv;i<nconv;i++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " SVD converged value (error) #%" PetscInt_FMT,its,i));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      PetscCall(PetscViewerASCIIPrintf(viewer," %g (%10.8e)\n",(double)sigma[i],(double)errest[i]));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SVDMonitorConvergedCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SVDMonitorConvergedDestroy(PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  if (!*vf) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree((*vf)->data));
  PetscCall(PetscViewerDestroy(&(*vf)->viewer));
  PetscCall(PetscFree(*vf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorFirstDrawLG - Plots the error estimate of the first unconverged
   approximation at each iteration of the singular value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor draw::draw_lg - activates SVDMonitorFirstDrawLG()

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorFirstDrawLG(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscReal      x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(svd->tol)-2,0.0));
  }
  if (nconv<nest) {
    x = (PetscReal)its;
    if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
    else y = 0.0;
    PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
    if (its <= 20 || !(its % 5) || svd->reason) {
      PetscCall(PetscDrawLGDraw(lg));
      PetscCall(PetscDrawLGSave(lg));
    }
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorFirstDrawLGCreate - Creates the plotter for the first error estimate.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorFirstDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"First Error Estimate","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorAllDrawLG - Plots the error estimate of all unconverged
   approximations at each iteration of the singular value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor_all draw::draw_lg - activates SVDMonitorAllDrawLG()

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorAllDrawLG(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg;
  PetscInt       i,n = PetscMin(svd->nsv,255);
  PetscReal      *x,*y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,n));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(svd->tol)-2,0.0));
  }
  PetscCall(PetscMalloc2(n,&x,n,&y));
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  PetscCall(PetscDrawLGAddPoint(lg,x,y));
  if (its <= 20 || !(its % 5) || svd->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscFree2(x,y));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorAllDrawLGCreate - Creates the plotter for all the error estimates.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorAllDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(PetscViewerMonitorLGSetUp(viewer,NULL,"All Error Estimates","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorConvergedDrawLG - Plots the number of converged eigenvalues
   at each iteration of the singular value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor_conv draw::draw_lg - activates SVDMonitorConvergedDrawLG()

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorConvergedDrawLG(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer      viewer = vf->viewer;
  PetscDrawLG      lg;
  PetscReal        x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerDrawGetDrawLG(viewer,0,&lg));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,0,svd->nsv));
  }
  x = (PetscReal)its;
  y = (PetscReal)svd->nconv;
  PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
  if (its <= 20 || !(its % 5) || svd->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SVDMonitorConvergedDrawLGCreate - Creates the plotter for the convergence history.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorConvergedDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
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

/*@C
   SVDMonitorConditioning - Print the condition number at each iteration of the singular
   value solver.

   Collective

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - (unused) number of converged singular triplets so far
.  sigma  - (unused) singular values
.  errest - (unused) error estimates
.  nest   - (unused) number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -svd_monitor_conditioning - activates SVDMonitorConditioning()

   Note:
   Works only for solvers that use a DS of type GSVD. The printed information corresponds
   to the maximum of the condition number of the two generated bidiagonal matrices.

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorAll(), SVDMonitorFirst(), SVDMonitorConverged()
@*/
PetscErrorCode SVDMonitorConditioning(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer viewer = vf->viewer;
  PetscBool   isgsvd;
  PetscReal   cond;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  PetscCall(PetscObjectTypeCompare((PetscObject)svd->ds,DSGSVD,&isgsvd));
  if (!isgsvd) PetscFunctionReturn(PETSC_SUCCESS);
  if (its==1 && ((PetscObject)svd)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Condition number of bidiagonal matrices for %s solve.\n",((PetscObject)svd)->prefix));
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel));
  PetscCall(DSCond(svd->ds,&cond));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " SVD condition number = %g\n",its,(double)cond));
  PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}
