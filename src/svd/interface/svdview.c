/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   The SVD routines related to various viewers
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/
#include <petscdraw.h>

/*@
   SVDView - Prints the SVD data structure.

   Collective

   Input Parameters:
+  svd - the singular value solver context
-  viewer - optional visualization context

   Options Database Key:
.  -svd_view -  Calls SVDView() at end of SVDSolve()

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

.seealso: EPSView()
@*/
PetscErrorCode SVDView(SVD svd,PetscViewer viewer)
{
  const char     *type=NULL;
  PetscBool      isascii,isshell,isexternal;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)svd),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(svd,1,viewer,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)svd,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(svd,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    if (svd->problem_type) {
      switch (svd->problem_type) {
        case SVD_STANDARD:    type = "(standard) singular value problem"; break;
        case SVD_GENERALIZED: type = "generalized singular value problem"; break;
        case SVD_HYPERBOLIC:  type = "hyperbolic singular value problem"; break;
      }
    } else type = "not yet set";
    PetscCall(PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  transpose mode: %s\n",svd->impltrans?"implicit":"explicit"));
    if (svd->which == SVD_LARGEST) PetscCall(PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: largest\n"));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: smallest\n"));
    if (svd->stop==SVD_STOP_THRESHOLD) PetscCall(PetscViewerASCIIPrintf(viewer,"  computing singular values %s the threshold: %g%s\n",svd->which==SVD_LARGEST?"above":"below",(double)svd->thres,svd->threlative?" (relative)":""));
    if (svd->nsv) PetscCall(PetscViewerASCIIPrintf(viewer,"  number of singular values (nsv): %" PetscInt_FMT "\n",svd->nsv));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %" PetscInt_FMT "\n",svd->ncv));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %" PetscInt_FMT "\n",svd->mpd));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %" PetscInt_FMT "\n",svd->max_it));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)svd->tol));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  convergence test: "));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    switch (svd->conv) {
    case SVD_CONV_ABS:
      PetscCall(PetscViewerASCIIPrintf(viewer,"absolute\n"));break;
    case SVD_CONV_REL:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the singular value\n"));break;
    case SVD_CONV_NORM:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the matrix norms\n"));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  computed matrix norms: norm(A)=%g",(double)svd->nrma));
      if (svd->isgeneralized) PetscCall(PetscViewerASCIIPrintf(viewer,", norm(B)=%g",(double)svd->nrmb));
      PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      break;
    case SVD_CONV_MAXIT:
      PetscCall(PetscViewerASCIIPrintf(viewer,"maximum number of iterations\n"));break;
    case SVD_CONV_USER:
      PetscCall(PetscViewerASCIIPrintf(viewer,"user-defined\n"));break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    if (svd->nini) PetscCall(PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %" PetscInt_FMT "\n",PetscAbs(svd->nini)));
    if (svd->ninil) PetscCall(PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial left space: %" PetscInt_FMT "\n",PetscAbs(svd->ninil)));
  } else PetscTryTypeMethod(svd,view,viewer);
  PetscCall(PetscObjectTypeCompareAny((PetscObject)svd,&isshell,SVDCROSS,SVDCYCLIC,""));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)svd,&isexternal,SVDSCALAPACK,SVDKSVD,SVDELEMENTAL,SVDPRIMME,""));
  if (!isshell && !isexternal) {
    PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO));
    if (!svd->V) PetscCall(SVDGetBV(svd,&svd->V,NULL));
    PetscCall(BVView(svd->V,viewer));
    if (!svd->ds) PetscCall(SVDGetDS(svd,&svd->ds));
    PetscCall(DSView(svd->ds,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDViewFromOptions - View from options

   Collective

   Input Parameters:
+  svd  - the singular value solver context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: SVDView(), SVDCreate()
@*/
PetscErrorCode SVDViewFromOptions(SVD svd,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)svd,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDConvergedReasonView - Displays the reason an SVD solve converged or diverged.

   Collective

   Input Parameters:
+  svd - the singular value solver context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -svd_converged_reason - print reason for convergence, and number of iterations

   Note:
   To change the format of the output call PetscViewerPushFormat(viewer,format) before
   this call. Use PETSC_VIEWER_DEFAULT for the default, use PETSC_VIEWER_FAILED to only
   display a reason if it fails. The latter can be set in the command line with
   -svd_converged_reason ::failed

   Level: intermediate

.seealso: SVDSetTolerances(), SVDGetIterationNumber(), SVDConvergedReasonViewFromOptions()
@*/
PetscErrorCode SVDConvergedReasonView(SVD svd,PetscViewer viewer)
{
  PetscBool         isAscii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (!viewer) viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)svd));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii));
  if (isAscii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel));
    if (svd->reason > 0 && format != PETSC_VIEWER_FAILED) PetscCall(PetscViewerASCIIPrintf(viewer,"%s SVD solve converged (%" PetscInt_FMT " singular triplet%s) due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)svd)->prefix?((PetscObject)svd)->prefix:"",svd->nconv,(svd->nconv>1)?"s":"",SVDConvergedReasons[svd->reason],svd->its));
    else if (svd->reason <= 0) PetscCall(PetscViewerASCIIPrintf(viewer,"%s SVD solve did not converge due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)svd)->prefix?((PetscObject)svd)->prefix:"",SVDConvergedReasons[svd->reason],svd->its));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDConvergedReasonViewFromOptions - Processes command line options to determine if/how
   the SVD converged reason is to be viewed.

   Collective

   Input Parameter:
.  svd - the singular value solver context

   Level: developer

.seealso: SVDConvergedReasonView()
@*/
PetscErrorCode SVDConvergedReasonViewFromOptions(SVD svd)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_converged_reason",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDConvergedReasonView(svd,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDErrorView_ASCII(SVD svd,SVDErrorType etype,PetscViewer viewer)
{
  PetscReal      error,sigma;
  PetscInt       i,j,numsv;

  PetscFunctionBegin;
  numsv = svd->stop==SVD_STOP_THRESHOLD? svd->nconv: svd->nsv;
  if (svd->nconv<numsv) {
    PetscCall(PetscViewerASCIIPrintf(viewer," Problem: less than %" PetscInt_FMT " singular values converged\n\n",numsv));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  for (i=0;i<numsv;i++) {
    PetscCall(SVDComputeError(svd,i,etype,&error));
    if (error>=5.0*svd->tol) {
      PetscCall(PetscViewerASCIIPrintf(viewer," Problem: some of the first %" PetscInt_FMT " relative errors are higher than the tolerance\n\n",numsv));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  PetscCall(PetscViewerASCIIPrintf(viewer," All %s%ssingular values computed up to the required tolerance:",svd->stop==SVD_STOP_THRESHOLD?"":"requested ",svd->isgeneralized?"generalized ":""));
  for (i=0;i<=(numsv-1)/8;i++) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n     "));
    for (j=0;j<PetscMin(8,numsv-8*i);j++) {
      PetscCall(SVDGetSingularTriplet(svd,8*i+j,&sigma,NULL,NULL));
      PetscCall(PetscViewerASCIIPrintf(viewer,"%.5f",(double)sigma));
      if (8*i+j+1<numsv) PetscCall(PetscViewerASCIIPrintf(viewer,", "));
    }
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDErrorView_DETAIL(SVD svd,SVDErrorType etype,PetscViewer viewer)
{
  PetscReal      error,sigma;
  PetscInt       i;
  char           ex[30],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!svd->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  switch (etype) {
    case SVD_ERROR_ABSOLUTE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex)," absolute error"));
      break;
    case SVD_ERROR_RELATIVE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex)," relative error"));
      break;
    case SVD_ERROR_NORM:
      if (svd->isgeneralized) PetscCall(PetscSNPrintf(ex,sizeof(ex)," ||r||/||[A;B]||"));
      else PetscCall(PetscSNPrintf(ex,sizeof(ex),"  ||r||/||A||"));
      break;
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s          sigma           %s\n%s",sep,ex,sep));
  for (i=0;i<svd->nconv;i++) {
    PetscCall(SVDGetSingularTriplet(svd,i,&sigma,NULL,NULL));
    PetscCall(SVDComputeError(svd,i,etype,&error));
    PetscCall(PetscViewerASCIIPrintf(viewer,"       % 6f          %12g\n",(double)sigma,(double)error));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s",sep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDErrorView_MATLAB(SVD svd,SVDErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)svd,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name));
  for (i=0;i<svd->nconv;i++) {
    PetscCall(SVDComputeError(svd,i,etype,&error));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDErrorView - Displays the errors associated with the computed solution
   (as well as the singular values).

   Collective

   Input Parameters:
+  svd    - the singular value solver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Keys:
+  -svd_error_absolute - print absolute errors of each singular triplet
.  -svd_error_relative - print relative errors of each singular triplet
-  -svd_error_norm     - print errors relative to the matrix norms of each singular triplet

   Notes:
   By default, this function checks the error of all singular triplets and prints
   the singular values if all of them are below the requested tolerance.
   If the viewer has format=PETSC_VIEWER_ASCII_INFO_DETAIL then a table with
   singular values and corresponding errors is printed.

   Level: intermediate

.seealso: SVDSolve(), SVDValuesView(), SVDVectorsView()
@*/
PetscErrorCode SVDErrorView(SVD svd,SVDErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)svd),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,3);
  PetscCheckSameComm(svd,1,viewer,3);
  SVDCheckSolved(svd,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (!isascii) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetFormat(viewer,&format));
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      PetscCall(SVDErrorView_ASCII(svd,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      PetscCall(SVDErrorView_DETAIL(svd,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      PetscCall(SVDErrorView_MATLAB(svd,etype,viewer));
      break;
    default:
      PetscCall(PetscInfo(svd,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective

   Input Parameter:
.  svd - the singular value solver context

   Level: developer

.seealso: SVDErrorView()
@*/
PetscErrorCode SVDErrorViewFromOptions(SVD svd)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_error_absolute",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDErrorView(svd,SVD_ERROR_ABSOLUTE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_error_relative",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDErrorView(svd,SVD_ERROR_RELATIVE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_error_norm",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDErrorView(svd,SVD_ERROR_NORM,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDValuesView_DRAW(SVD svd,PetscViewer viewer)
{
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im=0.0;
  PetscInt       i;

  PetscFunctionBegin;
  if (!svd->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
  PetscCall(PetscDrawSetTitle(draw,"Computed singular values"));
  PetscCall(PetscDrawSPCreate(draw,1,&drawsp));
  for (i=0;i<svd->nconv;i++) {
    re = svd->sigma[svd->perm[i]];
    PetscCall(PetscDrawSPAddPoint(drawsp,&re,&im));
  }
  PetscCall(PetscDrawSPDraw(drawsp,PETSC_TRUE));
  PetscCall(PetscDrawSPSave(drawsp));
  PetscCall(PetscDrawSPDestroy(&drawsp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDValuesView_BINARY(SVD svd,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscReal      *sv;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(svd->nconv,&sv));
  for (i=0;i<svd->nconv;i++) {
    k = svd->perm[i];
    sv[i] = svd->sigma[k];
  }
  PetscCall(PetscViewerBinaryWrite(viewer,sv,svd->nconv,PETSC_REAL));
  PetscCall(PetscFree(sv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_HDF5)
static PetscErrorCode SVDValuesView_HDF5(SVD svd,PetscViewer viewer)
{
  PetscInt       i,k,n,N;
  PetscMPIInt    rank;
  Vec            v;
  char           vname[30];
  const char     *ename;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)svd),&rank));
  N = svd->nconv;
  n = rank? 0: N;
  /* create a vector containing the singular values */
  PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)svd),n,N,&v));
  PetscCall(PetscObjectGetName((PetscObject)svd,&ename));
  PetscCall(PetscSNPrintf(vname,sizeof(vname),"sigma_%s",ename));
  PetscCall(PetscObjectSetName((PetscObject)v,vname));
  if (!rank) {
    for (i=0;i<svd->nconv;i++) {
      k = svd->perm[i];
      PetscCall(VecSetValue(v,i,svd->sigma[k],INSERT_VALUES));
    }
  }
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));
  PetscCall(VecView(v,viewer));
  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode SVDValuesView_ASCII(SVD svd,PetscViewer viewer)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(PetscViewerASCIIPrintf(viewer,"Singular values = \n"));
  for (i=0;i<svd->nconv;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"   %.5f\n",(double)svd->sigma[svd->perm[i]]));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDValuesView_MATLAB(SVD svd,PetscViewer viewer)
{
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)svd,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Sigma_%s = [\n",name));
  for (i=0;i<svd->nconv;i++) PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)svd->sigma[svd->perm[i]]));
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDValuesView - Displays the computed singular values in a viewer.

   Collective

   Input Parameters:
+  svd    - the singular value solver context
-  viewer - the viewer

   Options Database Key:
.  -svd_view_values - print computed singular values

   Level: intermediate

.seealso: SVDSolve(), SVDVectorsView(), SVDErrorView()
@*/
PetscErrorCode SVDValuesView(SVD svd,PetscViewer viewer)
{
  PetscBool         isascii,isdraw,isbinary;
  PetscViewerFormat format;
#if defined(PETSC_HAVE_HDF5)
  PetscBool         ishdf5;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)svd),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(svd,1,viewer,2);
  SVDCheckSolved(svd,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary));
#if defined(PETSC_HAVE_HDF5)
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERHDF5,&ishdf5));
#endif
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isdraw) PetscCall(SVDValuesView_DRAW(svd,viewer));
  else if (isbinary) PetscCall(SVDValuesView_BINARY(svd,viewer));
#if defined(PETSC_HAVE_HDF5)
  else if (ishdf5) PetscCall(SVDValuesView_HDF5(svd,viewer));
#endif
  else if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        PetscCall(SVDValuesView_ASCII(svd,viewer));
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        PetscCall(SVDValuesView_MATLAB(svd,viewer));
        break;
      default:
        PetscCall(PetscInfo(svd,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDValuesViewFromOptions - Processes command line options to determine if/how
   the computed singular values are to be viewed.

   Collective

   Input Parameter:
.  svd - the singular value solver context

   Level: developer

.seealso: SVDValuesView()
@*/
PetscErrorCode SVDValuesViewFromOptions(SVD svd)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_view_values",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDValuesView(svd,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDVectorsView - Outputs computed singular vectors to a viewer.

   Collective

   Input Parameters:
+  svd    - the singular value solver context
-  viewer - the viewer

   Options Database Key:
.  -svd_view_vectors - output singular vectors

   Note:
   Right and left singular vectors are interleaved, that is, the vectors are
   output in the following order V0, U0, V1, U1, V2, U2, ...

   Level: intermediate

.seealso: SVDSolve(), SVDValuesView(), SVDErrorView()
@*/
PetscErrorCode SVDVectorsView(SVD svd,PetscViewer viewer)
{
  PetscInt       i,k;
  Vec            x;
  char           vname[30];
  const char     *ename;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)svd),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(svd,1,viewer,2);
  SVDCheckSolved(svd,1);
  if (svd->nconv) {
    PetscCall(PetscObjectGetName((PetscObject)svd,&ename));
    PetscCall(SVDComputeVectors(svd));
    for (i=0;i<svd->nconv;i++) {
      k = svd->perm[i];
      PetscCall(PetscSNPrintf(vname,sizeof(vname),"V%" PetscInt_FMT "_%s",i,ename));
      PetscCall(BVGetColumn(svd->V,k,&x));
      PetscCall(PetscObjectSetName((PetscObject)x,vname));
      PetscCall(VecView(x,viewer));
      PetscCall(BVRestoreColumn(svd->V,k,&x));
      PetscCall(PetscSNPrintf(vname,sizeof(vname),"U%" PetscInt_FMT "_%s",i,ename));
      PetscCall(BVGetColumn(svd->U,k,&x));
      PetscCall(PetscObjectSetName((PetscObject)x,vname));
      PetscCall(VecView(x,viewer));
      PetscCall(BVRestoreColumn(svd->U,k,&x));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDVectorsViewFromOptions - Processes command line options to determine if/how
   the computed singular vectors are to be viewed.

   Collective

   Input Parameter:
.  svd - the singular value solver context

   Level: developer

.seealso: SVDVectorsView()
@*/
PetscErrorCode SVDVectorsViewFromOptions(SVD svd)
{
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,"-svd_view_vectors",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(SVDVectorsView(svd,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}
