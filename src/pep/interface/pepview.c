/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   The PEP routines related to various viewers
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/
#include <slepc/private/bvimpl.h>
#include <petscdraw.h>

/*@
   PEPView - Prints the `PEP` data structure.

   Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  viewer - optional visualization context

   Options Database Key:
.  -pep_view - calls `PEPView()` at end of `PEPSolve()`

   Notes:
   The available visualization contexts include
+     `PETSC_VIEWER_STDOUT_SELF` - standard output (default)
-     `PETSC_VIEWER_STDOUT_WORLD` - synchronized standard output where only the
         first process opens the file; all other processes send their data to the
         first one to print

   The user can open an alternative visualization context with `PetscViewerASCIIOpen()`
   to output to a specified file.

   Level: beginner

.seealso: [](ch:pep), `PEPCreate()`, `PEPViewFromOptions()`, `STView()`
@*/
PetscErrorCode PEPView(PEP pep,PetscViewer viewer)
{
  const char     *type=NULL;
  char           str[50];
  PetscBool      isascii,islinear,istrivial;
  PetscInt       i;
  PetscViewer    sviewer;
  MPI_Comm       child;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);

  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscObjectPrintClassNamePrefixType((PetscObject)pep,viewer));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscTryTypeMethod(pep,view,viewer);
    PetscCall(PetscViewerASCIIPopTab(viewer));
    if (pep->problem_type) {
      switch (pep->problem_type) {
        case PEP_GENERAL:    type = "general polynomial eigenvalue problem"; break;
        case PEP_HERMITIAN:  type = SLEPC_STRING_HERMITIAN " polynomial eigenvalue problem"; break;
        case PEP_HYPERBOLIC: type = "hyperbolic polynomial eigenvalue problem"; break;
        case PEP_GYROSCOPIC: type = "gyroscopic polynomial eigenvalue problem"; break;
      }
    } else type = "not yet set";
    PetscCall(PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  polynomial represented in %s basis\n",PEPBasisTypes[pep->basis]));
    switch (pep->scale) {
      case PEP_SCALE_NONE:
        break;
      case PEP_SCALE_SCALAR:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  parameter scaling enabled, with scaling factor=%g\n",(double)pep->sfactor));
        break;
      case PEP_SCALE_DIAGONAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  diagonal balancing enabled, with its=%" PetscInt_FMT " and lambda=%g\n",pep->sits,(double)pep->slambda));
        break;
      case PEP_SCALE_BOTH:
        PetscCall(PetscViewerASCIIPrintf(viewer,"  parameter scaling & diagonal balancing enabled, with scaling factor=%g, its=%" PetscInt_FMT " and lambda=%g\n",(double)pep->sfactor,pep->sits,(double)pep->slambda));
        break;
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: "));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    PetscCall(SlepcSNPrintfScalar(str,sizeof(str),pep->target,PETSC_FALSE));
    if (!pep->which) PetscCall(PetscViewerASCIIPrintf(viewer,"not yet set\n"));
    else switch (pep->which) {
      case PEP_WHICH_USER:
        PetscCall(PetscViewerASCIIPrintf(viewer,"user defined\n"));
        break;
      case PEP_TARGET_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (in magnitude)\n",str));
        break;
      case PEP_TARGET_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the real axis)\n",str));
        break;
      case PEP_TARGET_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the imaginary axis)\n",str));
        break;
      case PEP_LARGEST_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest eigenvalues in magnitude\n"));
        break;
      case PEP_SMALLEST_MAGNITUDE:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest eigenvalues in magnitude\n"));
        break;
      case PEP_LARGEST_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest real parts\n"));
        break;
      case PEP_SMALLEST_REAL:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest real parts\n"));
        break;
      case PEP_LARGEST_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"largest imaginary parts\n"));
        break;
      case PEP_SMALLEST_IMAGINARY:
        PetscCall(PetscViewerASCIIPrintf(viewer,"smallest imaginary parts\n"));
        break;
      case PEP_ALL:
        if (pep->inta || pep->intb) PetscCall(PetscViewerASCIIPrintf(viewer,"all eigenvalues in interval [%g,%g]\n",(double)pep->inta,(double)pep->intb));
        else PetscCall(PetscViewerASCIIPrintf(viewer,"all eigenvalues in the region\n"));
        break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of eigenvalues (nev): %" PetscInt_FMT "\n",pep->nev));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %" PetscInt_FMT "\n",pep->ncv));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %" PetscInt_FMT "\n",pep->mpd));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %" PetscInt_FMT "\n",pep->max_it));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)pep->tol));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  convergence test: "));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    switch (pep->conv) {
    case PEP_CONV_ABS:
      PetscCall(PetscViewerASCIIPrintf(viewer,"absolute\n"));break;
    case PEP_CONV_REL:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue\n"));break;
    case PEP_CONV_NORM:
      PetscCall(PetscViewerASCIIPrintf(viewer,"relative to the matrix norms\n"));
      if (pep->nrma) {
        PetscCall(PetscViewerASCIIPrintf(viewer,"  computed matrix norms: %g",(double)pep->nrma[0]));
        for (i=1;i<pep->nmat;i++) PetscCall(PetscViewerASCIIPrintf(viewer,", %g",(double)pep->nrma[i]));
        PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
      }
      break;
    case PEP_CONV_USER:
      PetscCall(PetscViewerASCIIPrintf(viewer,"user-defined\n"));break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  extraction type: %s\n",PEPExtractTypes[pep->extract]));
    if (pep->refine) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"  iterative refinement: %s, with %s scheme\n",PEPRefineTypes[pep->refine],PEPRefineSchemes[pep->scheme]));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  refinement stopping criterion: tol=%g, its=%" PetscInt_FMT "\n",(double)pep->rtol,pep->rits));
      if (pep->npart>1) PetscCall(PetscViewerASCIIPrintf(viewer,"  splitting communicator in %" PetscInt_FMT " partitions for refinement\n",pep->npart));
    }
    if (pep->nini) PetscCall(PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %" PetscInt_FMT "\n",PetscAbs(pep->nini)));
  } else PetscTryTypeMethod(pep,view,viewer);
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO));
  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  PetscCall(BVView(pep->V,viewer));
  if (!pep->rg) PetscCall(PEPGetRG(pep,&pep->rg));
  PetscCall(RGIsTrivial(pep->rg,&istrivial));
  if (!istrivial) PetscCall(RGView(pep->rg,viewer));
  PetscCall(PetscObjectTypeCompare((PetscObject)pep,PEPLINEAR,&islinear));
  if (!islinear) {
    if (!pep->ds) PetscCall(PEPGetDS(pep,&pep->ds));
    PetscCall(DSView(pep->ds,viewer));
  }
  PetscCall(PetscViewerPopFormat(viewer));
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(STView(pep->st,viewer));
  if (pep->refine!=PEP_REFINE_NONE) {
    PetscCall(PetscViewerASCIIPushTab(viewer));
    if (pep->npart>1) {
      if (pep->refinesubc->color==0) {
        PetscCall(PetscSubcommGetChild(pep->refinesubc,&child));
        PetscCall(PetscViewerASCIIGetStdout(child,&sviewer));
        PetscCall(KSPView(pep->refineksp,sviewer));
      }
    } else PetscCall(KSPView(pep->refineksp,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPViewFromOptions - View (print) a `PEP` object based on values in the options database.

   Collective

   Input Parameters:
+  pep  - the polynomial eigensolver context
.  obj  - optional object that provides the options prefix used to query the options database
-  name - command line option

   Level: intermediate

.seealso: [](ch:pep), `PEPView()`, `PEPCreate()`, `PetscObjectViewFromOptions()`
@*/
PetscErrorCode PEPViewFromOptions(PEP pep,PetscObject obj,const char name[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscCall(PetscObjectViewFromOptions((PetscObject)pep,obj,name));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPConvergedReasonView - Displays the reason a `PEP` solve converged or diverged.

   Collective

   Input Parameters:
+  pep - the polynomial eigensolver context
-  viewer - the viewer to display the reason

   Options Database Key:
.  -pep_converged_reason - print reason for convergence/divergence, and number of iterations

   Notes:
   Use `PEPConvergedReasonViewFromOptions()` to display the reason based on values
   in the options database.

   To change the format of the output call `PetscViewerPushFormat()` before this
   call. Use `PETSC_VIEWER_DEFAULT` for the default, or `PETSC_VIEWER_FAILED` to only
   display a reason if it fails. The latter can be set in the command line with
   `-pep_converged_reason ::failed`.

   Level: intermediate

.seealso: [](ch:pep), `PEPSetConvergenceTest()`, `PEPSetTolerances()`, `PEPGetIterationNumber()`, `PEPConvergedReasonViewFromOptions()`
@*/
PetscErrorCode PEPConvergedReasonView(PEP pep,PetscViewer viewer)
{
  PetscBool         isAscii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (!viewer) viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)pep));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii));
  if (isAscii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    PetscCall(PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel));
    if (pep->reason > 0 && format != PETSC_VIEWER_FAILED) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Polynomial eigensolve converged (%" PetscInt_FMT " eigenpair%s) due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)pep)->prefix?((PetscObject)pep)->prefix:"",pep->nconv,(pep->nconv>1)?"s":"",PEPConvergedReasons[pep->reason],pep->its));
    else if (pep->reason <= 0) PetscCall(PetscViewerASCIIPrintf(viewer,"%s Polynomial eigensolve did not converge due to %s; iterations %" PetscInt_FMT "\n",((PetscObject)pep)->prefix?((PetscObject)pep)->prefix:"",PEPConvergedReasons[pep->reason],pep->its));
    PetscCall(PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPConvergedReasonViewFromOptions - Processes command line options to determine if/how
   the `PEP` converged reason is to be viewed.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Level: intermediate

.seealso: [](ch:pep), `PEPConvergedReasonView()`
@*/
PetscErrorCode PEPConvergedReasonViewFromOptions(PEP pep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_converged_reason",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPConvergedReasonView(pep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPErrorView_ASCII(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i,j,k,nvals;

  PetscFunctionBegin;
  nvals = (pep->which==PEP_ALL)? pep->nconv: pep->nev;
  if (pep->which!=PEP_ALL && pep->nconv<pep->nev) {
    PetscCall(PetscViewerASCIIPrintf(viewer," Problem: less than %" PetscInt_FMT " eigenvalues converged\n\n",pep->nev));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  if (pep->which==PEP_ALL && !nvals) {
    PetscCall(PetscViewerASCIIPrintf(viewer," No eigenvalues have been found\n\n"));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  for (i=0;i<nvals;i++) {
    PetscCall(PEPComputeError(pep,i,etype,&error));
    if (error>=5.0*pep->tol) {
      PetscCall(PetscViewerASCIIPrintf(viewer," Problem: some of the first %" PetscInt_FMT " relative errors are higher than the tolerance\n\n",nvals));
      PetscFunctionReturn(PETSC_SUCCESS);
    }
  }
  if (pep->which==PEP_ALL) PetscCall(PetscViewerASCIIPrintf(viewer," Found %" PetscInt_FMT " eigenvalues, all of them computed up to the required tolerance:",nvals));
  else PetscCall(PetscViewerASCIIPrintf(viewer," All requested eigenvalues computed up to the required tolerance:"));
  for (i=0;i<=(nvals-1)/8;i++) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n     "));
    for (j=0;j<PetscMin(8,nvals-8*i);j++) {
      k = pep->perm[8*i+j];
      PetscCall(SlepcPrintEigenvalueASCII(viewer,pep->eigr[k],pep->eigi[k]));
      if (8*i+j+1<nvals) PetscCall(PetscViewerASCIIPrintf(viewer,", "));
    }
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPErrorView_DETAIL(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscInt       i;
  char           ex[30],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!pep->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  switch (etype) {
    case PEP_ERROR_ABSOLUTE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex),"   ||P(k)x||"));
      break;
    case PEP_ERROR_RELATIVE:
      PetscCall(PetscSNPrintf(ex,sizeof(ex),"||P(k)x||/||kx||"));
      break;
    case PEP_ERROR_BACKWARD:
      PetscCall(PetscSNPrintf(ex,sizeof(ex),"    eta(x,k)"));
      break;
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s            k             %s\n%s",sep,ex,sep));
  for (i=0;i<pep->nconv;i++) {
    PetscCall(PEPGetEigenpair(pep,i,&kr,&ki,NULL,NULL));
    PetscCall(PEPComputeError(pep,i,etype,&error));
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

static PetscErrorCode PEPErrorView_MATLAB(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)pep,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name));
  for (i=0;i<pep->nconv;i++) {
    PetscCall(PEPComputeError(pep,i,etype,&error));
    PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPErrorView - Displays the errors associated with the computed solution
   (as well as the eigenvalues).

   Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Keys:
+  -pep_error_absolute - print absolute errors of each eigenpair
.  -pep_error_relative - print relative errors of each eigenpair
-  -pep_error_backward - print backward errors of each eigenpair

   Notes:
   By default, this function checks the error of all eigenpairs and prints
   the eigenvalues if all of them are below the requested tolerance.
   If the viewer has format `PETSC_VIEWER_ASCII_INFO_DETAIL` then a table with
   eigenvalues and corresponding errors is printed.

   All the command-line options listed above admit an optional argument
   specifying the viewer type and options. For instance, use
   `-pep_error_relative :myerr.m:ascii_matlab` to save the errors in a file
   that can be executed in Matlab.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPValuesView()`, `PEPVectorsView()`
@*/
PetscErrorCode PEPErrorView(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,3);
  PetscCheckSameComm(pep,1,viewer,3);
  PEPCheckSolved(pep,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (!isascii) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(PetscViewerGetFormat(viewer,&format));
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      PetscCall(PEPErrorView_ASCII(pep,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      PetscCall(PEPErrorView_DETAIL(pep,etype,viewer));
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      PetscCall(PEPErrorView_MATLAB(pep,etype,viewer));
      break;
    default:
      PetscCall(PetscInfo(pep,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Level: developer

.seealso: [](ch:pep), `PEPErrorView()`
@*/
PetscErrorCode PEPErrorViewFromOptions(PEP pep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_absolute",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPErrorView(pep,PEP_ERROR_ABSOLUTE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_relative",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_backward",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPValuesView_DRAW(PEP pep,PetscViewer viewer)
{
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (!pep->nconv) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscViewerDrawGetDraw(viewer,0,&draw));
  PetscCall(PetscDrawSetTitle(draw,"Computed Eigenvalues"));
  PetscCall(PetscDrawSPCreate(draw,1,&drawsp));
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(pep->eigr[k]);
    im = PetscImaginaryPart(pep->eigr[k]);
#else
    re = pep->eigr[k];
    im = pep->eigi[k];
#endif
    PetscCall(PetscDrawSPAddPoint(drawsp,&re,&im));
  }
  PetscCall(PetscDrawSPDraw(drawsp,PETSC_TRUE));
  PetscCall(PetscDrawSPSave(drawsp));
  PetscCall(PetscDrawSPDestroy(&drawsp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPValuesView_BINARY(PEP pep,PetscViewer viewer)
{
#if defined(PETSC_HAVE_COMPLEX)
  PetscInt       i,k;
  PetscComplex   *ev;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_COMPLEX)
  PetscCall(PetscMalloc1(pep->nconv,&ev));
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    ev[i] = pep->eigr[k];
#else
    ev[i] = PetscCMPLX(pep->eigr[k],pep->eigi[k]);
#endif
  }
  PetscCall(PetscViewerBinaryWrite(viewer,ev,pep->nconv,PETSC_COMPLEX));
  PetscCall(PetscFree(ev));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_HDF5)
static PetscErrorCode PEPValuesView_HDF5(PEP pep,PetscViewer viewer)
{
  PetscInt       i,k,n,N;
  PetscMPIInt    rank;
  Vec            v;
  char           vname[30];
  const char     *ename;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rank));
  N = pep->nconv;
  n = rank? 0: N;
  /* create a vector containing the eigenvalues */
  PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)pep),n,N,&v));
  PetscCall(PetscObjectGetName((PetscObject)pep,&ename));
  PetscCall(PetscSNPrintf(vname,sizeof(vname),"eigr_%s",ename));
  PetscCall(PetscObjectSetName((PetscObject)v,vname));
  if (!rank) {
    for (i=0;i<pep->nconv;i++) {
      k = pep->perm[i];
      PetscCall(VecSetValue(v,i,pep->eigr[k],INSERT_VALUES));
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
    for (i=0;i<pep->nconv;i++) {
      k = pep->perm[i];
      PetscCall(VecSetValue(v,i,pep->eigi[k],INSERT_VALUES));
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

static PetscErrorCode PEPValuesView_ASCII(PEP pep,PetscViewer viewer)
{
  PetscInt       i,k;

  PetscFunctionBegin;
  PetscCall(PetscViewerASCIIPrintf(viewer,"Eigenvalues = \n"));
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
    PetscCall(PetscViewerASCIIPrintf(viewer,"   "));
    PetscCall(SlepcPrintEigenvalueASCII(viewer,pep->eigr[k],pep->eigi[k]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPValuesView_MATLAB(PEP pep,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscReal      re,im;
  const char     *name;

  PetscFunctionBegin;
  PetscCall(PetscObjectGetName((PetscObject)pep,&name));
  PetscCall(PetscViewerASCIIPrintf(viewer,"Lambda_%s = [\n",name));
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(pep->eigr[k]);
    im = PetscImaginaryPart(pep->eigr[k]);
#else
    re = pep->eigr[k];
    im = pep->eigi[k];
#endif
    if (im!=0.0) PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e%+18.16ei\n",(double)re,(double)im));
    else PetscCall(PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)re));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPValuesView - Displays the computed eigenvalues in a viewer.

   Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -pep_view_values - print computed eigenvalues

   Note:
   The command-line option listed above admits an optional argument
   specifying the viewer type and options. For instance, use
   `-pep_view_values :evals.m:ascii_matlab` to save the values in a file
   that can be executed in Matlab.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPVectorsView()`, `PEPErrorView()`
@*/
PetscErrorCode PEPValuesView(PEP pep,PetscViewer viewer)
{
  PetscBool         isascii,isdraw,isbinary;
  PetscViewerFormat format;
#if defined(PETSC_HAVE_HDF5)
  PetscBool         ishdf5;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);
  PEPCheckSolved(pep,1);
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw));
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary));
#if defined(PETSC_HAVE_HDF5)
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERHDF5,&ishdf5));
#endif
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isdraw) PetscCall(PEPValuesView_DRAW(pep,viewer));
  else if (isbinary) PetscCall(PEPValuesView_BINARY(pep,viewer));
#if defined(PETSC_HAVE_HDF5)
  else if (ishdf5) PetscCall(PEPValuesView_HDF5(pep,viewer));
#endif
  else if (isascii) {
    PetscCall(PetscViewerGetFormat(viewer,&format));
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        PetscCall(PEPValuesView_ASCII(pep,viewer));
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        PetscCall(PEPValuesView_MATLAB(pep,viewer));
        break;
      default:
        PetscCall(PetscInfo(pep,"Unsupported viewer format %s\n",PetscViewerFormats[format]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPValuesViewFromOptions - Processes command line options to determine if/how
   the computed eigenvalues are to be viewed.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Level: developer

.seealso: [](ch:pep), `PEPValuesView()`
@*/
PetscErrorCode PEPValuesViewFromOptions(PEP pep)
{
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_view_values",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPValuesView(pep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPVectorsView - Outputs computed eigenvectors to a viewer.

   Collective

   Input Parameters:
+  pep    - the polynomial eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -pep_view_vectors - output eigenvectors

   Notes:
   If PETSc was configured with real scalars, complex conjugate eigenvectors
   will be viewed as two separate real vectors, one containing the real part
   and another one containing the imaginary part.

   The command-line option listed above admits an optional argument
   specifying the viewer type and options. For instance, use
   `-pep_view_vectors binary:evecs.bin` to save the vectors in a binary file.

   Level: intermediate

.seealso: [](ch:pep), `PEPSolve()`, `PEPValuesView()`, `PEPErrorView()`
@*/
PetscErrorCode PEPVectorsView(PEP pep,PetscViewer viewer)
{
  PetscInt       i,k;
  Vec            xr,xi=NULL;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer));
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);
  PEPCheckSolved(pep,1);
  if (pep->nconv) {
    PetscCall(PEPComputeVectors(pep));
    PetscCall(BVCreateVec(pep->V,&xr));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(BVCreateVec(pep->V,&xi));
#endif
    for (i=0;i<pep->nconv;i++) {
      k = pep->perm[i];
      PetscCall(BV_GetEigenvector(pep->V,k,pep->eigi[k],xr,xi));
      PetscCall(SlepcViewEigenvector(viewer,xr,xi,"X",i,(PetscObject)pep));
    }
    PetscCall(VecDestroy(&xr));
#if !defined(PETSC_USE_COMPLEX)
    PetscCall(VecDestroy(&xi));
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPVectorsViewFromOptions - Processes command line options to determine if/how
   the computed eigenvectors are to be viewed.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Level: developer

.seealso: [](ch:pep), `PEPVectorsView()`
@*/
PetscErrorCode PEPVectorsViewFromOptions(PEP pep)
{
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(PETSC_SUCCESS);
  incall = PETSC_TRUE;
  PetscCall(PetscOptionsCreateViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_view_vectors",&viewer,&format,&flg));
  if (flg) {
    PetscCall(PetscViewerPushFormat(viewer,format));
    PetscCall(PEPVectorsView(pep,viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}
