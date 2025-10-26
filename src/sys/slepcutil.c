/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>            /*I "slepcsys.h" I*/

/*
   Internal functions used to register monitors.
 */
PetscErrorCode SlepcMonitorMakeKey_Internal(const char name[],PetscViewerType vtype,PetscViewerFormat format,char key[])
{
  PetscFunctionBegin;
  PetscCall(PetscStrncpy(key,name,PETSC_MAX_PATH_LEN));
  PetscCall(PetscStrlcat(key,":",PETSC_MAX_PATH_LEN));
  PetscCall(PetscStrlcat(key,vtype,PETSC_MAX_PATH_LEN));
  PetscCall(PetscStrlcat(key,":",PETSC_MAX_PATH_LEN));
  PetscCall(PetscStrlcat(key,PetscViewerFormats[format],PETSC_MAX_PATH_LEN));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PetscViewerAndFormatCreate_Internal(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **vf)
{
  PetscFunctionBegin;
  PetscCall(PetscViewerAndFormatCreate(viewer,format,vf));
  (*vf)->data = ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Given n vectors in V, this function gets references of them into W.
   If m<0 then some previous non-processed vectors remain in W and must be freed.
 */
PetscErrorCode SlepcBasisReference_Private(PetscInt n,Vec *V,PetscInt *m,Vec **W)
{
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<n;i++) PetscCall(PetscObjectReference((PetscObject)V[i]));
  PetscCall(SlepcBasisDestroy_Private(m,W));
  if (n>0) {
    PetscCall(PetscMalloc1(n,W));
    for (i=0;i<n;i++) (*W)[i] = V[i];
    *m = -n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Destroys a set of vectors.
   A negative value of m indicates that W contains vectors to be destroyed.
 */
PetscErrorCode SlepcBasisDestroy_Private(PetscInt *m,Vec **W)
{
  PetscInt       i;

  PetscFunctionBegin;
  if (*m<0) {
    for (i=0;i<-(*m);i++) PetscCall(VecDestroy(&(*W)[i]));
    PetscCall(PetscFree(*W));
  }
  *m = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   SlepcSNPrintfScalar - Prints a `PetscScalar` variable to a string of
   given length.

   Not Collective

   Input Parameters:
+  str - the string to print to
.  len - the length of str
.  val - scalar value to be printed
-  exp - to be used within an expression, print leading sign and parentheses
         in case of nonzero imaginary part

   Level: developer

.seealso: `PetscSNPrintf()`
@*/
PetscErrorCode SlepcSNPrintfScalar(char str[],size_t len,PetscScalar val,PetscBool exp)
{
#if defined(PETSC_USE_COMPLEX)
  PetscReal      re,im;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  if (exp) PetscCall(PetscSNPrintf(str,len,"%+g",(double)val));
  else PetscCall(PetscSNPrintf(str,len,"%g",(double)val));
#else
  re = PetscRealPart(val);
  im = PetscImaginaryPart(val);
  if (im!=0.0) {
    if (exp) PetscCall(PetscSNPrintf(str,len,"+(%g%+gi)",(double)re,(double)im));
    else PetscCall(PetscSNPrintf(str,len,"%g%+gi",(double)re,(double)im));
  } else {
    if (exp) PetscCall(PetscSNPrintf(str,len,"%+g",(double)re));
    else PetscCall(PetscSNPrintf(str,len,"%g",(double)re));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SlepcHasExternalPackage - Determine whether SLEPc has been configured with the
   given package.

   Not Collective

   Input Parameter:
.  pkg - external package name

   Output Parameter:
.  has - `PETSC_TRUE` if SLEPc is configured with the given package, else `PETSC_FALSE`

   Level: intermediate

   Notes:
   This is basically an alternative for `SLEPC_HAVE_XXX` whenever a preprocessor macro
   is not available/desirable, e.g., in Python.

   The external package name pkg is, e.g., `arpack` or `primme`.
   It should correspond to the name listed in  `./configure --help`.

   The lookup is case insensitive, i.e. looking for `ARPACK` or `arpack` is the same.

.seealso: [](sec:writing-prog), `EPSType`, `SVDType`
@*/
PetscErrorCode SlepcHasExternalPackage(const char pkg[], PetscBool *has)
{
  char           pkgstr[128],*loc;
  size_t         cnt;

  PetscFunctionBegin;
  PetscCall(PetscSNPrintfCount(pkgstr,sizeof(pkgstr),":%s:",&cnt,pkg));
  PetscCheck(cnt<sizeof(pkgstr),PETSC_COMM_SELF,PETSC_ERR_SUP,"Package name is too long: \"%s\"",pkg);
  PetscCall(PetscStrtolower(pkgstr));
#if defined(SLEPC_HAVE_PACKAGES)
  PetscCall(PetscStrstr(SLEPC_HAVE_PACKAGES,pkgstr,&loc));
#else
#error "SLEPC_HAVE_PACKAGES macro undefined. Please reconfigure"
#endif
  *has = loc? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   SlepcDebugViewMatrix - prints an array as a matrix, to be used from within a debugger.
   Output can be pasted to Matlab.

     nrows, ncols: size of printed matrix
     Xr, Xi: array to be printed (Xi not referenced in complex scalars)
     ldx: leading dimension
     s: name of Matlab variable
     filename: optionally write output to a file
 */
#if defined(PETSC_USE_DEBUG)
PetscErrorCode SlepcDebugViewMatrix(PetscInt nrows,PetscInt ncols,PetscScalar *Xr,PetscScalar *Xi,PetscInt ldx,const char *s,const char *filename)
{
  PetscInt       i,j;
  PetscViewer    viewer;

  PetscFunctionBegin;
  if (filename) PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer));
  else PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerASCIIPrintf(viewer,"%s = [\n",s));
  for (i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
#if defined(PETSC_USE_COMPLEX)
      PetscCall(PetscViewerASCIIPrintf(viewer,"%.18g+%.18gi ",(double)PetscRealPart(Xr[i+j*ldx]),(double)PetscImaginaryPart(Xr[i+j*ldx])));
#else
      if (Xi) PetscCall(PetscViewerASCIIPrintf(viewer,"%.18g+%.18gi ",(double)Xr[i+j*ldx],(double)Xi[i+j*ldx]));
      else PetscCall(PetscViewerASCIIPrintf(viewer,"%.18g ",(double)Xr[i+j*ldx]));
#endif
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(PetscViewerASCIIPrintf(viewer,"];\n"));
  if (filename) PetscCall(PetscViewerDestroy(&viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/*
   SlepcDebugSetMatlabStdout - sets Matlab format in stdout, to be used from within a debugger.
 */
#if defined(PETSC_USE_DEBUG) && !defined(PETSC_CLANG_STATIC_ANALYZER)
PETSC_UNUSED PetscErrorCode SlepcDebugSetMatlabStdout(void)
{
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer));
  PetscCall(PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
