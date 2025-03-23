/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcnep.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define nepmonitorset_                    NEPMONITORSET
#define nepmonitorall_                    NEPMONITORALL
#define nepmonitorfirst_                  NEPMONITORFIRST
#define nepmonitorconverged_              NEPMONITORCONVERGED
#define nepmonitorconvergedcreate_        NEPMONITORCONVERGEDCREATE
#define nepmonitorconvergeddestroy_       NEPMONITORCONVERGEDDESTROY
#define nepconvergedabsolute_             NEPCONVERGEDABSOLUTE
#define nepconvergedrelative_             NEPCONVERGEDRELATIVE
#define nepsetconvergencetestfunction_    NEPSETCONVERGENCETESTFUNCTION
#define nepsetstoppingtestfunction_       NEPSETSTOPPINGTESTFUNCTION
#define nepseteigenvaluecomparison_       NEPSETEIGENVALUECOMPARISON
#define nepsetfunction_                   NEPSETFUNCTION
#define nepgetfunction_                   NEPGETFUNCTION
#define nepsetjacobian_                   NEPSETJACOBIAN
#define nepgetjacobian_                   NEPGETJACOBIAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define nepmonitorset_                    nepmonitorset
#define nepmonitorall_                    nepmonitorall
#define nepmonitorfirst_                  nepmonitorfirst
#define nepmonitorconverged_              nepmonitorconverged
#define nepmonitorconvergedcreate_        nepmonitorconvergedcreate
#define nepmonitorconvergeddestroy_       nepmonitorconvergeddestroy
#define nepconvergedabsolute_             nepconvergedabsolute
#define nepconvergedrelative_             nepconvergedrelative
#define nepsetconvergencetestfunction_    nepsetconvergencetestfunction
#define nepsetstoppingtestfunction_       nepsetstoppingtestfunction
#define nepseteigenvaluecomparison_       nepseteigenvaluecomparison
#define nepsetfunction_                   nepsetfunction
#define nepgetfunction_                   nepgetfunction
#define nepsetjacobian_                   nepsetjacobian
#define nepgetjacobian_                   nepgetjacobian
#endif

/*
   These cannot be called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void nepmonitorall_(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void nepmonitorfirst_(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);
SLEPC_EXTERN void nepmonitorconverged_(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,PetscViewerAndFormat*,PetscErrorCode*);

SLEPC_EXTERN void nepmonitorconvergedcreate_(PetscViewer *vin,PetscViewerFormat *format,void *ctx,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(vin,v);
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = NEPMonitorConvergedCreate(v,*format,ctx,vf);
}

SLEPC_EXTERN void nepmonitorconvergeddestroy_(PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = NEPMonitorConvergedDestroy(vf);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
  PetscFortranCallbackId comparison;
  PetscFortranCallbackId function;
  PetscFortranCallbackId jacobian;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  PetscFortranCallbackId function_pgiptr;
  PetscFortranCallbackId jacobian_pgiptr;
#endif
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(NEP nep,PetscInt i,PetscInt nc,PetscScalar *er,PetscScalar *ei,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.monitor,(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&nep,&i,&nc,er,ei,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void **ctx)
{
  NEP nep = (NEP)*ctx;
  PetscObjectUseFortranCallback(nep,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.convergence,(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&nep,&eigr,&eigi,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void **ctx)
{
  NEP nep = (NEP)*ctx;
  PetscObjectUseFortranCallback(nep,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(nep,_cb.stopping,(NEP*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,NEPConvergedReason*,void*,PetscErrorCode*),(&nep,&its,&max_it,&nconv,&nev,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void **ctx)
{
  NEP nep = (NEP)*ctx;
  PetscObjectUseFortranCallback(nep,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode oureigenvaluecomparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  NEP eps = (NEP)ctx;
  PetscObjectUseFortranCallback(eps,_cb.comparison,(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*,PetscErrorCode*),(&ar,&ai,&br,&bi,r,_ctx,&ierr));
}

static PetscErrorCode ournepfunction(NEP nep,PetscScalar lambda,Mat T,Mat P,void *ctx)
{
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  void* ptr;
  PetscCall(PetscObjectGetFortranCallback((PetscObject)nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.function_pgiptr,NULL,&ptr));
#endif
  PetscObjectUseFortranCallback(nep,_cb.function,(NEP*,PetscScalar*,Mat*,Mat*,void*,PetscErrorCode* PETSC_F90_2PTR_PROTO_NOVAR),(&nep,&lambda,&T,&P,_ctx,&ierr PETSC_F90_2PTR_PARAM(ptr)));
}

static PetscErrorCode ournepjacobian(NEP nep,PetscScalar lambda,Mat J,void *ctx)
{
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  void* ptr;
  PetscCall(PetscObjectGetFortranCallback((PetscObject)nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.jacobian_pgiptr,NULL,&ptr));
#endif
  PetscObjectUseFortranCallback(nep,_cb.jacobian,(NEP*,PetscScalar*,Mat*,void*,PetscErrorCode* PETSC_F90_2PTR_PROTO_NOVAR),(&nep,&lambda,&J,_ctx,&ierr PETSC_F90_2PTR_PARAM(ptr)));
}

SLEPC_EXTERN void nepmonitorset_(NEP *nep,void (*monitor)(NEP*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorall_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorconverged_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorConverged,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)NEPMonitorConvergedDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)nepmonitorfirst_) {
    *ierr = NEPMonitorSet(*nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscCtxDestroyFn*)PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = NEPMonitorSet(*nep,ourmonitor,*nep,ourdestroy);
  }
}

SLEPC_EXTERN void nepconvergedabsolute_(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*);
SLEPC_EXTERN void nepconvergedrelative_(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*);

SLEPC_EXTERN void nepsetconvergencetestfunction_(NEP *nep,void (*func)(NEP*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)nepconvergedabsolute_) {
    *ierr = NEPSetConvergenceTest(*nep,NEP_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)nepconvergedrelative_) {
    *ierr = NEPSetConvergenceTest(*nep,NEP_CONV_REL);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = NEPSetConvergenceTestFunction(*nep,ourconvergence,*nep,ourconvdestroy);
  }
}

SLEPC_EXTERN void nepstoppingbasic_(NEP*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,NEPConvergedReason*,void*,PetscErrorCode*);

SLEPC_EXTERN void nepsetstoppingtestfunction_(NEP *nep,void (*func)(NEP*,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)nepstoppingbasic_) {
    *ierr = NEPSetStoppingTest(*nep,NEP_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = NEPSetStoppingTestFunction(*nep,ourstopping,*nep,ourstopdestroy);
  }
}

SLEPC_EXTERN void nepseteigenvaluecomparison_(NEP *nep,void (*func)(PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,void*),void* ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.comparison,(PetscVoidFunction)func,ctx); if (*ierr) return;
  *ierr = NEPSetEigenvalueComparison(*nep,oureigenvaluecomparison,*nep);
}

SLEPC_EXTERN void nepsetfunction_(NEP *nep,Mat *A,Mat *B,void (*func)(NEP*,PetscScalar*,Mat*,Mat*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr PETSC_F90_2PTR_PROTO(ptr))
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.function,(PetscVoidFunction)func,ctx);if (*ierr) return;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.function_pgiptr,NULL,ptr);if (*ierr) return;
#endif
  *ierr = NEPSetFunction(*nep,*A,*B,ournepfunction,NULL);
}

/* func is currently ignored from Fortran */
SLEPC_EXTERN void nepgetfunction_(NEP *nep,Mat *A,Mat *B,void *func,void **ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(ctx);
  CHKFORTRANNULLOBJECT(A);
  CHKFORTRANNULLOBJECT(B);
  *ierr = NEPGetFunction(*nep,A,B,NULL,NULL); if (*ierr) return;
  *ierr = PetscObjectGetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.function,NULL,ctx);
}

SLEPC_EXTERN void nepsetjacobian_(NEP *nep,Mat *J,void (*func)(NEP*,PetscScalar*,Mat*,void*,PetscErrorCode*),void *ctx,PetscErrorCode *ierr PETSC_F90_2PTR_PROTO(ptr))
{
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.jacobian,(PetscVoidFunction)func,ctx);if (*ierr) return;
#if defined(PETSC_HAVE_F90_2PTR_ARG)
  *ierr = PetscObjectSetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.jacobian_pgiptr,NULL,ptr);if (*ierr) return;
#endif
  *ierr = NEPSetJacobian(*nep,*J,ournepjacobian,NULL);
}

/* func is currently ignored from Fortran */
SLEPC_EXTERN void nepgetjacobian_(NEP *nep,Mat *J,void *func,void **ctx,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(ctx);
  CHKFORTRANNULLOBJECT(J);
  *ierr = NEPGetJacobian(*nep,J,NULL,NULL); if (*ierr) return;
  *ierr = PetscObjectGetFortranCallback((PetscObject)*nep,PETSC_FORTRAN_CALLBACK_CLASS,_cb.jacobian,NULL,ctx);
}
