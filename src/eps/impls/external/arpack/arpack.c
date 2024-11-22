/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the ARPACK package
*/

#include <slepc/private/epsimpl.h>
#include "arpack.h"

static PetscErrorCode EPSSetUp_ARPACK(EPS eps)
{
  PetscInt       ncv;
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;

  PetscFunctionBegin;
  EPSCheckDefinite(eps);
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev+2,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The value of ncv must be at least nev+2");
  } else eps->ncv = PetscMin(PetscMax(20,2*eps->nev+1),eps->n); /* set default value of ncv */
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(300,(PetscInt)(2*eps->n/eps->ncv));
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
  PetscCheck(eps->which!=EPS_WHICH_USER,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support user-defined ordering of eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING | EPS_FEATURE_TWOSIDED);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION);

  ncv = eps->ncv;
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(ar->rwork));
  PetscCall(PetscMalloc1(ncv,&ar->rwork));
  ar->lworkl = 3*ncv*ncv+5*ncv;
  PetscCall(PetscFree(ar->workev));
  PetscCall(PetscMalloc1(3*ncv,&ar->workev));
#else
  if (eps->ishermitian) {
    ar->lworkl = ncv*(ncv+8);
  } else {
    ar->lworkl = 3*ncv*ncv+6*ncv;
    PetscCall(PetscFree(ar->workev));
    PetscCall(PetscMalloc1(3*ncv,&ar->workev));
  }
#endif
  PetscCall(PetscFree(ar->workl));
  PetscCall(PetscMalloc1(ar->lworkl,&ar->workl));
  PetscCall(PetscFree(ar->select));
  PetscCall(PetscMalloc1(ncv,&ar->select));
  PetscCall(PetscFree(ar->workd));
  PetscCall(PetscMalloc1(3*eps->nloc,&ar->workd));

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(EPS_SetInnerProduct(eps));
  PetscCall(EPSSetWorkVecs(eps,2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_ARPACK(EPS eps)
{
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;
  char           bmat[1],howmny[] = "A";
  const char     *which;
  PetscInt       n,ld,iparam[11],ipntr[14],ido,info,nev,ncv,rvec;
#if !defined(PETSC_HAVE_MPIUNI) && !defined(PETSC_HAVE_MSMPI)
  MPI_Fint       fcomm;
#endif
  PetscScalar    sigmar,*pV,*resid;
  Vec            x,y,w = eps->work[0];
  Mat            A;
  PetscBool      isSinv,isShift;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    sigmai = 0.0;
#endif

  PetscFunctionBegin;
  nev = eps->nev;
  ncv = eps->ncv;
#if !defined(PETSC_HAVE_MPIUNI) && !defined(PETSC_HAVE_MSMPI)
  fcomm = MPI_Comm_c2f(PetscObjectComm((PetscObject)eps));
#endif
  n = eps->nloc;
  PetscCall(EPSGetStartVector(eps,0,NULL));
  PetscCall(BVSetActiveColumns(eps->V,0,0));  /* just for deflation space */
  PetscCall(BVCopyVec(eps->V,0,eps->work[1]));
  PetscCall(BVGetLeadingDimension(eps->V,&ld));
  PetscCall(BVGetArray(eps->V,&pV));
  PetscCall(VecGetArray(eps->work[1],&resid));

  ido  = 0;            /* first call to reverse communication interface */
  info = 1;            /* indicates an initial vector is provided */
  iparam[0] = 1;       /* use exact shifts */
  iparam[2] = eps->max_it;  /* max Arnoldi iterations */
  iparam[3] = 1;       /* blocksize */
  iparam[4] = 0;       /* number of converged Ritz values */

  /*
     Computational modes ([]=not supported):
            symmetric    non-symmetric    complex
        1     1  'I'        1  'I'         1  'I'
        2     3  'I'        3  'I'         3  'I'
        3     2  'G'        2  'G'         2  'G'
        4     3  'G'        3  'G'         3  'G'
        5   [ 4  'G' ]    [ 3  'G' ]
        6   [ 5  'G' ]    [ 4  'G' ]
   */
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&isSinv));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isShift));
  PetscCall(STGetShift(eps->st,&sigmar));
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatCreateVecsEmpty(A,&x,&y));

  if (isSinv) {
    /* shift-and-invert mode */
    iparam[6] = 3;
    if (eps->ispositive) bmat[0] = 'G';
    else bmat[0] = 'I';
  } else if (isShift && eps->ispositive) {
    /* generalized shift mode with B positive definite */
    iparam[6] = 2;
    bmat[0] = 'G';
  } else {
    /* regular mode */
    PetscCheck(!eps->ishermitian || !eps->isgeneralized,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Spectral transformation not supported by ARPACK hermitian solver");
    iparam[6] = 1;
    bmat[0] = 'I';
  }

#if !defined(PETSC_USE_COMPLEX)
  if (eps->ishermitian) {
    switch (eps->which) {
      case EPS_TARGET_MAGNITUDE:
      case EPS_LARGEST_MAGNITUDE:  which = "LM"; break;
      case EPS_SMALLEST_MAGNITUDE: which = "SM"; break;
      case EPS_TARGET_REAL:
      case EPS_LARGEST_REAL:       which = "LA"; break;
      case EPS_SMALLEST_REAL:      which = "SA"; break;
      default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Wrong value of eps->which");
    }
  } else {
#endif
    switch (eps->which) {
      case EPS_TARGET_MAGNITUDE:
      case EPS_LARGEST_MAGNITUDE:  which = "LM"; break;
      case EPS_SMALLEST_MAGNITUDE: which = "SM"; break;
      case EPS_TARGET_REAL:
      case EPS_LARGEST_REAL:       which = "LR"; break;
      case EPS_SMALLEST_REAL:      which = "SR"; break;
      case EPS_TARGET_IMAGINARY:
      case EPS_LARGEST_IMAGINARY:  which = "LI"; break;
      case EPS_SMALLEST_IMAGINARY: which = "SI"; break;
      default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Wrong value of eps->which");
    }
#if !defined(PETSC_USE_COMPLEX)
  }
#endif

  do {

#if !defined(PETSC_USE_COMPLEX)
    if (eps->ishermitian) {
      PetscStackCallExternalVoid("ARPACKsaupd",ARPACKsaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    } else {
      PetscStackCallExternalVoid("ARPACKnaupd",ARPACKnaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    }
#else
    PetscStackCallExternalVoid("ARPACKnaupd",ARPACKnaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,ar->rwork,&info));
#endif

    if (ido == -1 || ido == 1 || ido == 2) {
      if (ido == 1 && iparam[6] == 3 && bmat[0] == 'G') PetscCall(VecPlaceArray(x,&ar->workd[ipntr[2]-1])); /* special case for shift-and-invert with B semi-positive definite*/
      else PetscCall(VecPlaceArray(x,&ar->workd[ipntr[0]-1]));
      PetscCall(VecPlaceArray(y,&ar->workd[ipntr[1]-1]));

      if (ido == -1) {
        /* Y = OP * X for the initialization phase to
           force the starting vector into the range of OP */
        PetscCall(STApply(eps->st,x,y));
      } else if (ido == 2) {
        /* Y = B * X */
        PetscCall(BVApplyMatrix(eps->V,x,y));
      } else { /* ido == 1 */
        if (iparam[6] == 3 && bmat[0] == 'G') {
          /* Y = OP * X for shift-and-invert with B semi-positive definite */
          PetscCall(STMatSolve(eps->st,x,y));
        } else if (iparam[6] == 2) {
          /* X=A*X Y=B^-1*X for shift with B positive definite */
          PetscCall(MatMult(A,x,y));
          if (sigmar != 0.0) {
            PetscCall(BVApplyMatrix(eps->V,x,w));
            PetscCall(VecAXPY(y,sigmar,w));
          }
          PetscCall(VecCopy(y,x));
          PetscCall(STMatSolve(eps->st,x,y));
        } else {
          /* Y = OP * X */
          PetscCall(STApply(eps->st,x,y));
        }
        PetscCall(BVOrthogonalizeVec(eps->V,y,NULL,NULL,NULL));
      }

      PetscCall(VecResetArray(x));
      PetscCall(VecResetArray(y));
    } else PetscCheck(ido==99,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Internal error in ARPACK reverse communication interface (ido=%" PetscInt_FMT ")",ido);

  } while (ido != 99);

  eps->nconv = iparam[4];
  eps->its = iparam[2];

  PetscCheck(info!=3,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"No shift could be applied in xxAUPD. Try increasing the size of NCV relative to NEV");
  PetscCheck(info==0 || info==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error reported by ARPACK subroutine xxAUPD (%" PetscInt_FMT ")",info);

  rvec = PETSC_TRUE;

  if (eps->nconv > 0) {
#if !defined(PETSC_USE_COMPLEX)
    if (eps->ishermitian) {
      PetscStackCallExternalVoid("ARPACKseupd",ARPACKseupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,pV,&ld,&sigmar,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    } else {
      PetscStackCallExternalVoid("ARPACKneupd",ARPACKneupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,eps->eigi,pV,&ld,&sigmar,&sigmai,ar->workev,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    }
#else
    PetscStackCallExternalVoid("ARPACKneupd",ARPACKneupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,pV,&ld,&sigmar,ar->workev,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&ld,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,ar->rwork,&info));
#endif
    PetscCheck(info==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error reported by ARPACK subroutine xxEUPD (%" PetscInt_FMT ")",info);
  }

  PetscCall(BVRestoreArray(eps->V,&pV));
  PetscCall(VecRestoreArray(eps->work[1],&resid));
  if (eps->nconv >= eps->nev) eps->reason = EPS_CONVERGED_TOL;
  else eps->reason = EPS_DIVERGED_ITS;

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSBackTransform_ARPACK(EPS eps)
{
  PetscBool      isSinv;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&isSinv));
  if (!isSinv) PetscCall(EPSBackTransform_Default(eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_ARPACK(EPS eps)
{
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscFree(ar->workev));
  PetscCall(PetscFree(ar->workl));
  PetscCall(PetscFree(ar->select));
  PetscCall(PetscFree(ar->workd));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree(ar->rwork));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_ARPACK(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_ARPACK(EPS eps)
{
  EPS_ARPACK     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->ops->solve          = EPSSolve_ARPACK;
  eps->ops->setup          = EPSSetUp_ARPACK;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->destroy        = EPSDestroy_ARPACK;
  eps->ops->reset          = EPSReset_ARPACK;
  eps->ops->backtransform  = EPSBackTransform_ARPACK;
  PetscFunctionReturn(PETSC_SUCCESS);
}
