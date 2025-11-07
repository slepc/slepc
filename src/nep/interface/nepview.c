/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to various viewers
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/
#include <slepc/private/bvimpl.h>
#include <petscdraw.h>

/*@
   NEPView - Prints the `NEP` data structure.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  viewer - optional visualization context

   Options Database Key:
.  -nep_view - calls `NEPView()` at end of `NEPSolve()`

   Notes:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Level: beginner

.seealso: [](ch:nep), `NEPCreate()`, `NEPViewFromOptions()`, `FNView()`
@*/
PetscErrorCode NEPView(NEP nep,PetscViewer viewer)
{
  const char     *type=NULL;
  char           str[50];
  PetscInt       i;
  PetscBool      isascii,istrivial;
  PetscViewer    sviewer;
  MPI_Comm       child;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)nep,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(nep,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    if (nep->problem_type) {
      switch (nep->problem_type) {
        case NEP_GENERAL:  type = "general nonlinear eigenvalue problem"; break;
        case NEP_RATIONAL: type = "rational eigenvalue problem"; break;
      }
    } else type = "not yet set";
    PetscCall(PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type));
    if (nep->fui) {
      switch (nep->fui) {
      case NEP_USER_INTERFACE_CALLBACK:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  nonlinear operator from user callbacks\n"));
        break;
      case NEP_USER_INTERFACE_SPLIT:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  nonlinear operator in split form\n"));
        PetscCall(PetscViewerASCIIPrintf(viewer,"    number of terms: %" PetscInt_FMT "\n",nep->nt));
        PetscCall(PetscViewerASCIIPrintf(viewer,"    nonzero pattern of the matrices: %s\n",MatStructures[nep->mstr]));
        break;
      }
    } else PetscCall(PetscViewerASCIIPrintf(viewer,"  nonlinear operator not specified yet\n"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: "));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),nep->target,PETSC_FALSE));
    if (!nep->which) PetscCall(PetscViewerASCIIPrintf(viewer,"not yet set\n"));
    else switch (nep->which) {
      case NEP_WHICH_USER:
        PetscCall(PetscViewerASCIIPrintf(viewer,"user defined\n"));
        break;
      case NEP_TARGET_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (in magnitude)\n",str));
        break;
      case NEP_TARGET_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the real axis)\n",str));
        break;
      case NEP_TARGET_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the imaginary axis)\n",str));
        break;
      case NEP_LARGEST_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest eigenvalues in magnitude\n"));
        break;
      case NEP_SMALLEST_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest eigenvalues in magnitude\n"));
        break;
      case NEP_LARGEST_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest real parts\n"));
        break;
      case NEP_SMALLEST_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest real parts\n"));
        break;
      case NEP_LARGEST_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest imaginary parts\n"));
        break;
      case NEP_SMALLEST_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest imaginary parts\n"));
        break;
      case NEP_ALL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"all eigenvalues in the region\n"));
        break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    if (nep->twosided) PetscCall(PetscViewerASCIIPrintf(viewer,"  using two-sided variant (for left eigenvectors)\n"));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of eigenvalues (nev): %" PetscInt_FMT "\n",nep->nev));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %" PetscInt_FMT "\n",nep->ncv));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %" PetscInt_FMT "\n",nep->mpd));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %" PetscInt_FMT "\n",nep->max_it));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)nep->tol));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  convergence test: "));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    switch (nep->conv) {
    case NEP_CONV_ABS:
      PetscCall(PetscViewerASCIIPrintf(viewer,"absolute\n"));break;
    case NEP_CONV_REL:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue\n"));break;
    case NEP_CONV_NORM:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the matrix norms\n"));
      if (nep->nrma) {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  computed matrix norms: %g",(double)nep->nrma[0]));
        for (i=1;i<nep->nt;i++) PetscCall(PetscViewerASCIIPrintf(viewer,", %g",(double)nep->nrma[i]));
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
      break;
    case NEP_CONV_USER:
      PetscCall(PetscViewerASCIIPrintf(viewer,"user-defined\n"));break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    if (nep->refine) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  iterative refinement: %s, with %s scheme\n",NEPRefineTypes[nep->refine],NEPRefineSchemes[nep->scheme]));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  refinement stopping criterion: tol=%g, its=%" PetscInt_FMT "\n",(double)nep->rtol,nep->rits));
      if (nep->npart>1) PetscCall(PetscViewerASCIIPrintf(viewer,"  splitting communicator in %" PetscInt_FMT " partitions for refinement\n",nep->npart));
    }
    if (nep->nini) PetscCall(PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %" PetscInt_FMT "\n",PetscAbs(nep->nini)));
  } else PetscTryTypeMethod(nep,view,viewer);
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO));
  if (!nep->V) PetscCall(NEPGetBV(nep,&nep->V));
  PetscCall(BVView(nep->V,viewer));
  if (!nep->rg) PetscCall(NEPGetRG(nep,&nep->rg));
  PetscCall(RGIsTrivial(nep->rg,&istrivial));
  if (!istrivial) PetscCall(RGView(nep->rg,viewer));
  if (nep->useds) {
    if (!nep->ds) PetscCall(NEPGetDS(nep,&nep->ds));
    PetscCall(DSView(nep->ds,viewer));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  if (nep->refine!=NEP_REFINE_NONE) {
    PetscCall(PetscViewerASCIIPushTab(viewer));
    if (nep->npart>1) {
      if (nep->refinesubc->color==0) {
        PetscCall(PetscSubcommGetChild(nep->refinesubc,&child));
        PetscCall(PetscViewerASCIIGetStdout(child,&sviewer));
        PetscCall(KSPView(nep->refineksp,sviewer));
      }
    } else PetscCall(KSPView(nep->refineksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPViewFromOptions - View (print) a `NEP` object based on values in the options database.

   Collective

   Input Parameters:
+  nep  - the nonlinear eigensolver context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](ch:nep), `NEPView()`, `NEPCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode NEPViewFromOptions(NEP nep,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)nep,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPConvergedReasonView - Displays the reason a `NEP` solve converged or diverged.

   Collective

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  viewer - the viewer to display the reason

   Options Database Key:
.  -nep_converged_reason - print reason for convergence/divergence, and number of iterations

   Notes:
   Use `NEPConvergedReasonViewFromOptions()` to display the reason based on values
   in the options database.

   To change the format of the output call `PetscViewerPushFormat()` before this
   call. Use `PETSC_VIEWER_DEFAULT` for the default, or `PETSC_VIEWER_FAILED` to only
   display a reason if it fails. The latter can be set in the command line with
   `-nep_converged_reason ::failed`.

   Level: intermediate

.seealso: [](ch:nep), `NEPSetConvergenceTest()`, `NEPSetTolerances()`, `NEPGetIterationNumber()`, `NEPConvergedReasonViewFromOptions()`
@*/
PetscErrorCode NEPConvergedReasonView(NEP nep,PetscViewer viewer)
{
  PetscBool         isAscii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (!viewer) viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)nep));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii));
  if (isAscii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel));
    if (nep->reason > 0 && format != PETSC_VIEWER_FAILED) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Nonlinear eigensolve converged (%" PetscInt_FMT " eigenpair%s) due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)nep)->prefix?((PetscObject)nep)->prefix:"",nep->nconv,(nep->nconv>1)?"s":"",NEPConvergedReasons[nep->reason],nep->its));
    else if (nep->reason <= 0) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Nonlinear eigensolve did not converge due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)nep)->prefix?((PetscObject)nep)->prefix:"",NEPConvergedReasons[nep->reason],nep->its));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPConvergedReasonViewFromOptions - Processes command line options to determine if/how
   the `NEP` converged reason is to be viewed.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: intermediate

.seealso: [](ch:nep), `NEPConvergedReasonView()`
@*/
PetscErrorCode NEPConvergedReasonViewFromOptions(NEP nep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_converged_reason",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPConvergedReasonView(nep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPErrorView_ASCII(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i,j,k,nvals;

  PetscFunctionBegin;
  nvals = (nep->which==NEP_ALL)? nep->nconv: nep->nev;
  if (nep->which!=NEP_ALL && nep->nconv<nvals) {
    PetscCall(PetscViewerASCIIPrintf(viewer," Problem: less than %" PetscInt_FMT " eigenvalues converged\n\n",nep->nev));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  if (nep->which==NEP_ALL && !nvals) {
    PetscCall(PetscViewerASCIIPrintf(viewer," No eigenvalues have been found\n\n"));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  for (i=0;i<nvals;i++) {
    PetscCall(NEPComputeError(nep,i,etype,&error));
    if (error>=5.0*nep->tol) {
      PetscCall(PetscViewerASCIIPrintf(viewer," Problem: some of the first %" PetscInt_FMT " relative errors are higher than the tolerance\n\n",nvals));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  if (nep->which==NEP_ALL) PetscCall(PetscViewerASCIIPrintf(viewer," Found %" PetscInt_FMT " eigenvalues, all of them computed up to the required tolerance:",nvals));
  else PetscCall(PetscViewerASCIIPrintf(viewer," All requested eigenvalues computed up to the required tolerance:"));
  for (i=0;i<=(nvals-1)/8;i++) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n     "));
    for (j=0;j<PetscMin(8,nvals-8*i);j++) {
      k = nep->perm[8*i+j];
      PetscCall(SlepcPrintEigenvalueASCII(viewer,nep->eigr[k],nep->eigi[k]));
      if (8*i+j+1<nvals) PetscCall(PetscViewerASCIIPrintf(viewer,", "));
    }
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPErrorView_DETAIL(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscInt       i;
  char           ex[30],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!nep->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  switch (etype) {
    case NEP_ERROR_ABSOLUTE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex),"    ||T(k)x||"));
      break;
    case NEP_ERROR_RELATIVE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex)," ||T(k)x||/||kx||"));
      break;
    case NEP_ERROR_BACKWARD:
      PetscCall(PetscSNPrintf(ex,sizeof(ex),"    eta(x,k)"));
      break;
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s            k             %s\n%s",sep,ex,sep));
  for (i=0;i<nep->nconv;i++) {
    PetscCall(NEPGetEigenpair(nep,i,&kr,&ki,NULL,NULL));
    PetscCall(NEPComputeError(nep,i,etype,&error));
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif
    if (im!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"  % 9f%+9fi      %12g\n",(double)re,(double)im,(double)error));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"    % 12f           %12g\n",(double)re,(double)error));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s",sep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPErrorView_MATLAB(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)nep,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name));
  for (i=0;i<nep->nconv;i++) {
    PetscCall(NEPComputeError(nep,i,etype,&error));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPErrorView - Displays the errors associated with the computed solution
   (as well as the eigenvalues).

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Keys:
+  -nep_error_absolute - print absolute errors of each eigenpair
.  -nep_error_relative - print relative errors of each eigenpair
-  -nep_error_backward - print backward errors of each eigenpair

   Notes:
   By default, this function checks the error of all eigenpairs and prints
   the eigenvalues if all of them are below the requested tolerance.
   If the viewer has format `PETSC_VIEWER_ASCII_INFO_DETAIL` then a table with
   eigenvalues and corresponding errors is printed.

   All the command-line options listed above admit an optional argument
   specifying the viewer type and options. For instance, use
   `-nep_error_relative :myerr.m:ascii_matlab` to save the errors in a file
   that can be executed in Matlab.

   Level: intermediate

.seealso: [](ch:nep), `NEPSolve()`, `NEPValuesView()`, `NEPVectorsView()`
@*/
PetscErrorCode NEPErrorView(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,3);
  PetscCheckSameComm(nep,1,viewer,3);
  NEPCheckSolved(nep,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (!isascii) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetFormat(viewer,&format));
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      PetscCall(NEPErrorView_ASCII(nep,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      PetscCall(NEPErrorView_DETAIL(nep,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      PetscCall(NEPErrorView_MATLAB(nep,etype,viewer));
      break;
    default:
      PetscCall(PetscInfo(nep,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer

.seealso: [](ch:nep), `NEPErrorView()`
@*/
PetscErrorCode NEPErrorViewFromOptions(NEP nep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_absolute",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPErrorView(nep,NEP_ERROR_ABSOLUTE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_relative",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPErrorView(nep,NEP_ERROR_RELATIVE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_backward",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPErrorView(nep,NEP_ERROR_BACKWARD,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPValuesView_DRAW(NEP nep,PetscViewer viewer)
{
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (!nep->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
  PetscCall(PetscDrawSetTitle(draw,"Computed Eigenvalues"));
  PetscCall(PetscDrawSPCreate(draw,1,&drawsp));
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(nep->eigr[k]);
    im = PetscImaginaryPart(nep->eigr[k]);
#else
    re = nep->eigr[k];
    im = nep->eigi[k];
#endif
    PetscCall(PetscDrawSPAddPoint(drawsp,&re,&im));
  }
  PetscCall(PetscDrawSPDraw(drawsp,PETSC_TRUE));
  PetscCall(PetscDrawSPSave(drawsp));
  PetscCall(PetscDrawSPDestroy(&drawsp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPValuesView_BINARY(NEP nep,PetscViewer viewer)
{
#if defined(PETSC_HAVE_COMPLEX)
  PetscInt       i,k;
  PetscComplex   *ev;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_COMPLEX)
  PetscCall(PetscMalloc1(nep->nconv,&ev));
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    ev[i] = nep->eigr[k];
#else
    ev[i] = PetscCMPLX(nep->eigr[k],nep->eigi[k]);
#endif
  }
  PetscCall(PetscViewerBinaryWrite(viewer,ev,nep->nconv,PETSC_COMPLEX));
  PetscCall(PetscFree(ev));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_HDF5)
static PetscErrorCode NEPValuesView_HDF5(NEP nep,PetscViewer viewer)
{
  PetscInt       i,k,n,N;
  PetscMPIInt    rank;
  Vec            v;
  char           vname[30];
  const char     *ename;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)nep),&rank));
  N = nep->nconv;
  n = rank? 0: N;
  /* create a vector containing the eigenvalues */
  PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)nep),n,N,&v));
  PetscCall(PetscObjectGetName((PetscObject)nep,&ename));
  PetscCall(PetscSNPrintf(vname,sizeof(vname),"eigr_%s",ename));
  PetscCall(PetscObjectSetName((PetscObject)v,vname));
  if (!rank) {
    for (i=0;i<nep->nconv;i++) {
      k = nep->perm[i];
      PetscCall(VecSetValue(v,i,nep->eigr[k],INSERT_VALUES));
    }
  }
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));
  PetscCall(VecView(v,viewer));
#if !defined(PETSC_USE_COMPLEX)
  /* in real scalars write the imaginary part as a separate vector */
  PetscCall(PetscSNPrintf(vname,sizeof(vname),"eigi_%s",ename));
  PetscCall(PetscObjectSetName((PetscObject)v,vname));
  if (!rank) {
    for (i=0;i<nep->nconv;i++) {
      k = nep->perm[i];
      PetscCall(VecSetValue(v,i,nep->eigi[k],INSERT_VALUES));
    }
  }
  PetscCall(VecAssemblyBegin(v));
  PetscCall(VecAssemblyEnd(v));
  PetscCall(VecView(v,viewer));
#endif
  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

static PetscErrorCode NEPValuesView_ASCII(NEP nep,PetscViewer viewer)
{
  PetscInt       i,k;

  PetscFunctionBegin;
  PetscCall(PetscViewerASCIIPrintf(viewer,"Eigenvalues = \n"));
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
    PetscCall(PetscViewerASCIIPrintf(viewer,"   "));
    PetscCall(SlepcPrintEigenvalueASCII(viewer,nep->eigr[k],nep->eigi[k]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPValuesView_MATLAB(NEP nep,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscReal      re,im;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)nep,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Lambda_%s = [\n",name));
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(nep->eigr[k]);
    im = PetscImaginaryPart(nep->eigr[k]);
#else
    re = nep->eigr[k];
    im = nep->eigi[k];
#endif
    if (im!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e%+18.16ei\n",(double)re,(double)im));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)re));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPValuesView - Displays the computed eigenvalues in a viewer.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -nep_view_values - print computed eigenvalues

   Note:
   The command-line option listed above admits an optional argument
   specifying the viewer type and options. For instance, use
   `-nep_view_values :evals.m:ascii_matlab` to save the values in a file
   that can be executed in Matlab.

   Level: intermediate

.seealso: [](ch:nep), `NEPSolve()`, `NEPVectorsView()`, `NEPErrorView()`
@*/
PetscErrorCode NEPValuesView(NEP nep,PetscViewer viewer)
{
  PetscBool         isascii,isdraw,isbinary;
  PetscViewerFormat format;
#if defined(PETSC_HAVE_HDF5)
  PetscBool         ishdf5;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);
  NEPCheckSolved(nep,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary));
#if defined(PETSC_HAVE_HDF5)
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERHDF5,&ishdf5));
#endif
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isdraw) PetscCall(NEPValuesView_DRAW(nep,viewer));
  else if (isbinary) PetscCall(NEPValuesView_BINARY(nep,viewer));
#if defined(PETSC_HAVE_HDF5)
  else if (ishdf5) PetscCall(NEPValuesView_HDF5(nep,viewer));
#endif
  else if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        PetscCall(NEPValuesView_ASCII(nep,viewer));
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        PetscCall(NEPValuesView_MATLAB(nep,viewer));
        break;
      default:
        PetscCall(PetscInfo(nep,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPValuesViewFromOptions - Processes command line options to determine if/how
   the computed eigenvalues are to be viewed.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer

.seealso: [](ch:nep), `NEPValuesView()`
@*/
PetscErrorCode NEPValuesViewFromOptions(NEP nep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_view_values",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPValuesView(nep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPVectorsView - Outputs computed eigenvectors to a viewer.

   Collective

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -nep_view_vectors - output eigenvectors

   Notes:
   If PETSc was configured with real scalars, complex conjugate eigenvectors
   will be viewed as two separate real vectors, one containing the real part
   and another one containing the imaginary part.

   If left eigenvectors were computed with a two-sided eigensolver, the right
   and left eigenvectors are interleaved, that is, the vectors are output in
   the following order\: `X0, Y0, X1, Y1, X2, Y2, ...`

   The command-line option listed above admits an optional argument
   specifying the viewer type and options. For instance, use
   `-nep_view_vectors binary:evecs.bin` to save the vectors in a binary file.

   Level: intermediate

.seealso: [](ch:nep), `NEPSolve()`, `NEPValuesView()`, `NEPErrorView()`
@*/
PetscErrorCode NEPVectorsView(NEP nep,PetscViewer viewer)
{
  PetscInt       i,k;
  Vec            xr,xi=NULL;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);
  NEPCheckSolved(nep,1);
  if (nep->nconv) {
    PetscCall(NEPComputeVectors(nep));
    PetscCall(BVCreateVec(nep->V,&xr));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(BVCreateVec(nep->V,&xi));
#endif
    for (i=0;i<nep->nconv;i++) {
      k = nep->perm[i];
      PetscCall(BV_GetEigenvector(nep->V,k,nep->eigi[k],xr,xi));
      PetscCall(SlepcViewEigenvector(viewer,xr,xi,"X",i,(PetscObject)nep));
      if (nep->twosided) {
        PetscCall(BV_GetEigenvector(nep->W,k,nep->eigi[k],xr,xi));
        PetscCall(SlepcViewEigenvector(viewer,xr,xi,"Y",i,(PetscObject)nep));
      }
    }
    PetscCall(VecDestroy(&xr));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(VecDestroy(&xi));
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   NEPVectorsViewFromOptions - Processes command line options to determine if/how
   the computed eigenvectors are to be viewed.

   Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer

.seealso: [](ch:nep), `NEPVectorsView()`
@*/
PetscErrorCode NEPVectorsViewFromOptions(NEP nep)
{
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_view_vectors",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(NEPVectorsView(nep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}
