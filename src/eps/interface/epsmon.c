/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to monitors
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <petscdraw.h>

PetscErrorCode EPSMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],const char metric[],PetscInt l,const char *names[],int x,int y,int m,int n,PetscDrawLG *lgctx)
{
  PetscDraw      draw;
  PetscDrawAxis  axis;
  PetscDrawLG    lg;

  PetscFunctionBegin;
  PetscCall(PetscDrawCreate(comm,host,label,x,y,m,n,&draw));
  PetscCall(PetscDrawSetFromOptions(draw));
  PetscCall(PetscDrawLGCreate(draw,l,&lg));
  if (names) PetscCall(PetscDrawLGSetLegend(lg,names));
  PetscCall(PetscDrawLGSetFromOptions(lg));
  PetscCall(PetscDrawLGGetAxis(lg,&axis));
  PetscCall(PetscDrawAxisSetLabels(axis,"Convergence","Iteration",metric));
  PetscCall(PetscDrawDestroy(&draw));
  *lgctx = lg;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode EPSMonitor(EPS eps,PetscInt it,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest)
{
  PetscInt       i,n = eps->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) PetscCall((*eps->monitor[i])(eps,it,nconv,eigr,eigi,errest,nest,eps->monitorcontext[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested eigenpair.

   Logically Collective

   Input Parameters:
+  eps     - eigensolver context obtained from EPSCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL),
             see PetscCtxDestroyFn for the calling sequence

   Calling sequence of monitor:
$  PetscErrorCode monitor(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx)
+  eps    - eigensolver context obtained from EPSCreate()
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by EPSMonitorSet()

   Options Database Keys:
+    -eps_monitor        - print only the first error estimate
.    -eps_monitor_all    - print error estimates at each iteration
.    -eps_monitor_conv   - print the eigenvalue approximations only when
      convergence has been reached
.    -eps_monitor draw::draw_lg - sets line graph monitor for the first unconverged
      approximate eigenvalue
.    -eps_monitor_all draw::draw_lg - sets line graph monitor for all unconverged
      approximate eigenvalues
.    -eps_monitor_conv draw::draw_lg - sets line graph monitor for convergence history
-    -eps_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to EPSMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   EPSMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: EPSMonitorFirst(), EPSMonitorAll(), EPSMonitorCancel()
@*/
PetscErrorCode EPSMonitorSet(EPS eps,PetscErrorCode (*monitor)(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscCtxDestroyFn *monitordestroy)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscCheck(eps->numbermonitors<MAXEPSMONITORS,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Too many EPS monitors set");
  eps->monitor[eps->numbermonitors]           = monitor;
  eps->monitorcontext[eps->numbermonitors]    = (void*)mctx;
  eps->monitordestroy[eps->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSMonitorCancel - Clears all monitors for an EPS object.

   Logically Collective

   Input Parameters:
.  eps - eigensolver context obtained from EPSCreate()

   Options Database Key:
.    -eps_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to EPSMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorCancel(EPS eps)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  for (i=0; i<eps->numbermonitors; i++) {
    if (eps->monitordestroy[i]) PetscCall((*eps->monitordestroy[i])(&eps->monitorcontext[i]));
  }
  eps->numbermonitors = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSGetMonitorContext - Gets the monitor context, as set by
   EPSMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSGetMonitorContext(EPS eps,void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  *(void**)ctx = eps->monitorcontext[0];
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode EPSMonitorPrintEval(EPS eps,PetscViewer viewer,PetscScalar er,PetscScalar ei)
{
  PetscFunctionBegin;
  if (eps->problem_type==EPS_HEP || eps->problem_type==EPS_GHEP || eps->problem_type==EPS_BSE) PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)PetscRealPart(er)));
  else {
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er)));
#else
    PetscCall(PetscViewerASCIIPrintf(viewer," %g",(double)er));
    if (ei!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei));
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the eigensolver.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor - activates EPSMonitorFirst()

   Level: intermediate

.seealso: EPSMonitorSet(), EPSMonitorAll(), EPSMonitorConverged()
@*/
PetscErrorCode EPSMonitorFirst(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)eps)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)eps)->prefix));
  if (nconv<nest) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)eps)->tablevel));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " EPS nconv=%" PetscInt_FMT " first unconverged value (error)",its,nconv));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    er = eigr[nconv]; ei = eigi[nconv];
    PetscCall(STBackTransform(eps->st,1,&er,&ei));
    PetscCall(EPSMonitorPrintEval(eps,viewer,er,ei));
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[nconv]));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)eps)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorAll - Print the current approximate values and
   error estimates at each iteration of the eigensolver.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor_all - activates EPSMonitorAll()

   Level: intermediate

.seealso: EPSMonitorSet(), EPSMonitorFirst(), EPSMonitorConverged()
@*/
PetscErrorCode EPSMonitorAll(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)eps)->tablevel));
  if (its==1 && ((PetscObject)eps)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)eps)->prefix));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " EPS nconv=%" PetscInt_FMT " Values (Errors)",its,nconv));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
  for (i=0;i<nest;i++) {
    er = eigr[i]; ei = eigi[i];
    PetscCall(STBackTransform(eps->st,1,&er,&ei));
    PetscCall(EPSMonitorPrintEval(eps,viewer,er,ei));
    PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)",(double)errest[i]));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
  PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)eps)->tablevel));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor_conv - activates EPSMonitorConverged()

   Level: intermediate

.seealso: EPSMonitorSet(), EPSMonitorFirst(), EPSMonitorAll()
@*/
PetscErrorCode EPSMonitorConverged(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer = vf->viewer;
  SlepcConvMon   ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  ctx = (SlepcConvMon)vf->data;
  if (its==1 && ((PetscObject)eps)->prefix) PetscCall(PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)eps)->prefix));
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    PetscCall(PetscViewerPushFormat(viewer,vf->format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)eps)->tablevel));
    for (i=ctx->oldnconv;i<nconv;i++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%3" PetscInt_FMT " EPS converged value (error) #%" PetscInt_FMT,its,i));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
      er = eigr[i]; ei = eigi[i];
      PetscCall(STBackTransform(eps->st,1,&er,&ei));
      PetscCall(EPSMonitorPrintEval(eps,viewer,er,ei));
      PetscCall(PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[i]));
      PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    }
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)eps)->tablevel));
    PetscCall(PetscViewerPopFormat(viewer));
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSMonitorConvergedCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSMonitorConvergedDestroy(PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  if (!*vf) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree((*vf)->data));
  PetscCall(PetscViewerDestroy(&(*vf)->viewer));
  PetscCall(PetscDrawLGDestroy(&(*vf)->lg));
  PetscCall(PetscFree(*vf));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorFirstDrawLG - Plots the error estimate of the first unconverged
   approximation at each iteration of the eigensolver.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor draw::draw_lg - activates EPSMonitorFirstDrawLG()

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorFirstDrawLG(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg = vf->lg;
  PetscReal      x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(eps->tol)-2,0.0));
  }
  if (nconv<nest) {
    x = (PetscReal)its;
    if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
    else y = 0.0;
    PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
    if (its <= 20 || !(its % 5) || eps->reason) {
      PetscCall(PetscDrawLGDraw(lg));
      PetscCall(PetscDrawLGSave(lg));
    }
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorFirstDrawLGCreate - Creates the plotter for the first error estimate.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorFirstDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(EPSMonitorLGCreate(PetscObjectComm((PetscObject)viewer),NULL,"First Error Estimate","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300,&(*vf)->lg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorAllDrawLG - Plots the error estimate of all unconverged
   approximations at each iteration of the eigensolver.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor_all draw::draw_lg - activates EPSMonitorAllDrawLG()

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorAllDrawLG(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer    viewer = vf->viewer;
  PetscDrawLG    lg = vf->lg;
  PetscInt       i,n = PetscMin(eps->nev,255);
  PetscReal      *x,*y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,n));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,PetscLog10Real(eps->tol)-2,0.0));
  }
  PetscCall(PetscMalloc2(n,&x,n,&y));
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  PetscCall(PetscDrawLGAddPoint(lg,x,y));
  if (its <= 20 || !(its % 5) || eps->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscFree2(x,y));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorAllDrawLGCreate - Creates the plotter for all the error estimates.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorAllDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscCall(EPSMonitorLGCreate(PetscObjectComm((PetscObject)viewer),NULL,"All Error Estimates","Log Error Estimate",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300,&(*vf)->lg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorConvergedDrawLG - Plots the number of converged eigenvalues
   at each iteration of the eigensolver.

   Collective

   Input Parameters:
+  eps    - eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Options Database Key:
.  -eps_monitor_conv draw::draw_lg - activates EPSMonitorConvergedDrawLG()

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorConvergedDrawLG(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscViewer      viewer = vf->viewer;
  PetscDrawLG      lg = vf->lg;
  PetscReal        x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  PetscCall(PetscViewerPushFormat(viewer,vf->format));
  if (its==1) {
    PetscCall(PetscDrawLGReset(lg));
    PetscCall(PetscDrawLGSetDimension(lg,1));
    PetscCall(PetscDrawLGSetLimits(lg,1,1,0,eps->nev));
  }
  x = (PetscReal)its;
  y = (PetscReal)eps->nconv;
  PetscCall(PetscDrawLGAddPoint(lg,&x,&y));
  if (its <= 20 || !(its % 5) || eps->reason) {
    PetscCall(PetscDrawLGDraw(lg));
    PetscCall(PetscDrawLGSave(lg));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   EPSMonitorConvergedDrawLGCreate - Creates the plotter for the convergence history.

   Collective

   Input Parameters:
+  viewer - the viewer
.  format - the viewer format
-  ctx    - an optional user context

   Output Parameter:
.  vf     - the viewer and format context

   Level: intermediate

.seealso: EPSMonitorSet()
@*/
PetscErrorCode EPSMonitorConvergedDrawLGCreate(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  SlepcConvMon   mctx;

  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  PetscCall(PetscNew(&mctx));
  mctx->ctx = ctx;
  (*vf)->data = (void*)mctx;
  PetscCall(EPSMonitorLGCreate(PetscObjectComm((PetscObject)viewer),NULL,"Convergence History","Number Converged",1,NULL,PETSC_DECIDE,PETSC_DECIDE,400,300,&(*vf)->lg));
  PetscFunctionReturn(PETSC_SUCCESS);
}
