/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcsvd.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define svddestroy_                       SVDDESTROY
#define svdmonitorset_                    SVDMONITORSET
#define svdmonitorall_                    SVDMONITORALL
#define svdmonitorfirst_                  SVDMONITORFIRST
#define svdmonitorconditioning_           SVDMONITORCONDITIONING
#define svdmonitorconverged_              SVDMONITORCONVERGED
#define svdmonitorconvergedcreate_        SVDMONITORCONVERGEDCREATE
#define svdmonitorconvergeddestroy_       SVDMONITORCONVERGEDDESTROY
#define svdconvergedabsolute_             SVDCONVERGEDABSOLUTE
#define svdconvergedrelative_             SVDCONVERGEDRELATIVE
#define svdconvergednorm_                 SVDCONVERGEDNORM
#define svdconvergedmaxit_                SVDCONVERGEDMAXIT
#define svdsetconvergencetestfunction_    SVDSETCONVERGENCETESTFUNCTION
#define svdsetstoppingtestfunction_       SVDSETSTOPPINGTESTFUNCTION
#define svdgetdimensions000_              SVDGETDIMENSIONS000
#define svdgetdimensions100_              SVDGETDIMENSIONS100
#define svdgetdimensions010_              SVDGETDIMENSIONS010
#define svdgetdimensions001_              SVDGETDIMENSIONS001
#define svdgetdimensions110_              SVDGETDIMENSIONS110
#define svdgetdimensions011_              SVDGETDIMENSIONS011
#define svdgetdimensions101_              SVDGETDIMENSIONS101
#define svdgetsingulartriplet0_           SVDGETSINGULARTRIPLET0
#define svdgetsingulartriplet1_           SVDGETSINGULARTRIPLET1
#define svdgettolerances00_               SVDGETTOLERANCES00
#define svdgettolerances10_               SVDGETTOLERANCES10
#define svdgettolerances01_               SVDGETTOLERANCES01
#define svdsetinitialspaces00_            SVDSETINITIALSPACES00
#define svdsetinitialspaces01_            SVDSETINITIALSPACES01
#define svdsetinitialspaces10_            SVDSETINITIALSPACES10
#define svdsetinitialspaces11_            SVDSETINITIALSPACES11
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define svddestroy_                       svddestroy
#define svdmonitorset_                    svdmonitorset
#define svdmonitorall_                    svdmonitorall
#define svdmonitorfirst_                  svdmonitorfirst
#define svdmonitorconditioning_           svdmonitorconditioning
#define svdmonitorconverged_              svdmonitorconverged
#define svdmonitorconvergedcreate_        svdmonitorconvergedcreate
#define svdmonitorconvergeddestroy_       svdmonitorconvergeddestroy
#define svdconvergedabsolute_             svdconvergedabsolute
#define svdconvergedrelative_             svdconvergedrelative
#define svdconvergednorm_                 svdconvergednorm
#define svdconvergedmaxit_                svdconvergedmaxit
#define svdsetconvergencetestfunction_    svdsetconvergencetestfunction
#define svdsetstoppingtestfunction_       svdsetstoppingtestfunction
#define svdgetdimensions000_              svdgetdimensions000
#define svdgetdimensions100_              svdgetdimensions100
#define svdgetdimensions010_              svdgetdimensions010
#define svdgetdimensions001_              svdgetdimensions001
#define svdgetdimensions110_              svdgetdimensions110
#define svdgetdimensions011_              svdgetdimensions011
#define svdgetdimensions101_              svdgetdimensions101
#define svdgetsingulartriplet0_           svdgetsingulartriplet0
#define svdgetsingulartriplet1_           svdgetsingulartriplet1
#define svdgettolerances00_               svdgettolerances00
#define svdgettolerances10_               svdgettolerances10
#define svdgettolerances01_               svdgettolerances01
#define svdsetinitialspaces00_            svdsetinitialspaces00
#define svdsetinitialspaces01_            svdsetinitialspaces01
#define svdsetinitialspaces10_            svdsetinitialspaces10
#define svdsetinitialspaces11_            svdsetinitialspaces11
#endif

/*
   These are not usually called from Fortran but allow Fortran users
   to transparently set these monitors from .F code
*/
SLEPC_EXTERN void svdmonitorall_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorAll(*svd,*it,*nconv,sigma,errest,*nest,*vf);
}

SLEPC_EXTERN void svdmonitorfirst_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorFirst(*svd,*it,*nconv,sigma,errest,*nest,*vf);
}

SLEPC_EXTERN void svdmonitorconditioning_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorConditioning(*svd,*it,*nconv,sigma,errest,*nest,*vf);
}

SLEPC_EXTERN void svdmonitorconverged_(SVD *svd,PetscInt *it,PetscInt *nconv,PetscReal *sigma,PetscReal *errest,PetscInt *nest,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorConverged(*svd,*it,*nconv,sigma,errest,*nest,*vf);
}

SLEPC_EXTERN void svdmonitorconvergedcreate_(PetscViewer *vin,PetscViewerFormat *format,void *ctx,PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  PetscViewer v;
  PetscPatchDefaultViewers_Fortran(vin,v);
  CHKFORTRANNULLOBJECT(ctx);
  *ierr = SVDMonitorConvergedCreate(v,*format,ctx,vf);
}

SLEPC_EXTERN void svdmonitorconvergeddestroy_(PetscViewerAndFormat **vf,PetscErrorCode *ierr)
{
  *ierr = SVDMonitorConvergedDestroy(vf);
}

static struct {
  PetscFortranCallbackId monitor;
  PetscFortranCallbackId monitordestroy;
  PetscFortranCallbackId convergence;
  PetscFortranCallbackId convdestroy;
  PetscFortranCallbackId stopping;
  PetscFortranCallbackId stopdestroy;
} _cb;

/* These are not extern C because they are passed into non-extern C user level functions */
static PetscErrorCode ourmonitor(SVD svd,PetscInt i,PetscInt nc,PetscReal *sigma,PetscReal *d,PetscInt l,void* ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.monitor,(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,void*,PetscErrorCode*),(&svd,&i,&nc,sigma,d,&l,_ctx,&ierr));
}

static PetscErrorCode ourdestroy(void** ctx)
{
  SVD svd = (SVD)*ctx;
  PetscObjectUseFortranCallback(svd,_cb.monitordestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourconvergence(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.convergence,(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*),(&svd,&sigma,&res,errest,_ctx,&ierr));
}

static PetscErrorCode ourconvdestroy(void *ctx)
{
  SVD svd = (SVD)ctx;
  PetscObjectUseFortranCallback(svd,_cb.convdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

static PetscErrorCode ourstopping(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)
{
  PetscObjectUseFortranCallback(svd,_cb.stopping,(SVD*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,SVDConvergedReason*,void*,PetscErrorCode*),(&svd,&its,&max_it,&nconv,&nsv,reason,_ctx,&ierr));
}

static PetscErrorCode ourstopdestroy(void *ctx)
{
  SVD svd = (SVD)ctx;
  PetscObjectUseFortranCallback(svd,_cb.stopdestroy,(void*,PetscErrorCode*),(_ctx,&ierr));
}

SLEPC_EXTERN void svddestroy_(SVD *svd,PetscErrorCode *ierr)
{
  PETSC_FORTRAN_OBJECT_F_DESTROYED_TO_C_NULL(svd);
  *ierr = SVDDestroy(svd); if (*ierr) return;
  PETSC_FORTRAN_OBJECT_C_NULL_TO_F_DESTROYED(svd);
}

SLEPC_EXTERN void svdmonitorset_(SVD *svd,void (*monitor)(SVD*,PetscInt*,PetscInt*,PetscReal*,PetscReal*,PetscInt*,void*,PetscErrorCode*),void *mctx,void (*monitordestroy)(void *,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(mctx);
  CHKFORTRANNULLFUNCTION(monitordestroy);
  if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorall_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorAll,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorconverged_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorConverged,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))SVDMonitorConvergedDestroy);
  } else if ((PetscVoidFunction)monitor == (PetscVoidFunction)svdmonitorfirst_) {
    *ierr = SVDMonitorSet(*svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorFirst,*(PetscViewerAndFormat**)mctx,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitor,(PetscVoidFunction)monitor,mctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.monitordestroy,(PetscVoidFunction)monitordestroy,mctx); if (*ierr) return;
    *ierr = SVDMonitorSet(*svd,ourmonitor,*svd,ourdestroy);
  }
}

SLEPC_EXTERN void svdconvergedabsolute_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedAbsolute(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdconvergedrelative_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedRelative(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdconvergednorm_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedNorm(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdconvergedmaxit_(SVD *svd,PetscReal *sigma,PetscReal *res,PetscReal *errest,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDConvergedMaxIt(*svd,*sigma,*res,errest,ctx);
}

SLEPC_EXTERN void svdsetconvergencetestfunction_(SVD *svd,void (*func)(SVD*,PetscReal*,PetscReal*,PetscReal*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedabsolute_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_ABS);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedrelative_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_REL);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergednorm_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_NORM);
  } else if ((PetscVoidFunction)func == (PetscVoidFunction)svdconvergedmaxit_) {
    *ierr = SVDSetConvergenceTest(*svd,SVD_CONV_MAXIT);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convergence,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.convdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetConvergenceTestFunction(*svd,ourconvergence,*svd,ourconvdestroy);
  }
}

SLEPC_EXTERN void svdstoppingbasic_(SVD *svd,PetscInt *its,PetscInt *max_it,PetscInt *nconv,PetscInt *nsv,SVDConvergedReason *reason,void *ctx,PetscErrorCode *ierr)
{
  *ierr = SVDStoppingBasic(*svd,*its,*max_it,*nconv,*nsv,reason,ctx);
}

SLEPC_EXTERN void svdsetstoppingtestfunction_(SVD *svd,void (*func)(SVD*,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*,PetscErrorCode*),void* ctx,void (*destroy)(void*,PetscErrorCode*),PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(ctx);
  CHKFORTRANNULLFUNCTION(destroy);
  if ((PetscVoidFunction)func == (PetscVoidFunction)svdstoppingbasic_) {
    *ierr = SVDSetStoppingTest(*svd,SVD_STOP_BASIC);
  } else {
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopping,(PetscVoidFunction)func,ctx); if (*ierr) return;
    *ierr = PetscObjectSetFortranCallback((PetscObject)*svd,PETSC_FORTRAN_CALLBACK_CLASS,&_cb.stopdestroy,(PetscVoidFunction)destroy,ctx); if (*ierr) return;
    *ierr = SVDSetStoppingTestFunction(*svd,ourstopping,*svd,ourstopdestroy);
  }
}

SLEPC_EXTERN void svdgetdimensions_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  CHKFORTRANNULLINTEGER(nev);
  CHKFORTRANNULLINTEGER(ncv);
  CHKFORTRANNULLINTEGER(mpd);
  *ierr = SVDGetDimensions(*svd,nev,ncv,mpd);
}

SLEPC_EXTERN void svdgetdimensions000_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions100_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions010_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions001_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions110_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions011_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetdimensions101_(SVD *svd,PetscInt *nev,PetscInt *ncv,PetscInt *mpd,PetscErrorCode *ierr)
{
  svdgetdimensions_(svd,nev,ncv,mpd,ierr);
}

SLEPC_EXTERN void svdgetsingulartriplet_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(sigma);
  *ierr = SVDGetSingularTriplet(*svd,*i,sigma,*u,*v);
}

SLEPC_EXTERN void svdgetsingulartriplet0_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v, PetscErrorCode *ierr)
{
  svdgetsingulartriplet_(svd,i,sigma,u,v,ierr);
}

SLEPC_EXTERN void svdgetsingulartriplet1_(SVD *svd,PetscInt *i,PetscReal *sigma,Vec *u,Vec *v, PetscErrorCode *ierr)
{
  svdgetsingulartriplet_(svd,i,sigma,u,v,ierr);
}

SLEPC_EXTERN void svdgettolerances_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(tol);
  CHKFORTRANNULLINTEGER(maxits);
  *ierr = SVDGetTolerances(*svd,tol,maxits);
}

SLEPC_EXTERN void svdgettolerances00_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdgettolerances10_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdgettolerances01_(SVD *svd,PetscReal *tol,PetscInt *maxits,PetscErrorCode *ierr)
{
  svdgettolerances_(svd,tol,maxits,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(isr);
  CHKFORTRANNULLOBJECT(isl);
  *ierr = SVDSetInitialSpaces(*svd,*nr,isr,*nl,isl);
}

SLEPC_EXTERN void svdsetinitialspaces00_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces01_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces10_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}

SLEPC_EXTERN void svdsetinitialspaces11_(SVD *svd,PetscInt *nr,Vec *isr,PetscInt *nl,Vec *isl,PetscErrorCode *ierr)
{
  svdsetinitialspaces_(svd,nr,isr,nl,isl,ierr);
}
