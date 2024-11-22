/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to eigensolvers in EVSL.
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include <evsl.h>

#define PetscCallEVSL(func, ...) do {                                                   \
    PetscStackPushExternal(PetscStringize(func));                                                      \
    int evsl_ierr_ = func(__VA_ARGS__);                                              \
    PetscStackPop;                                                                             \
    PetscCheck(!evsl_ierr_,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error calling %s: error code %d",PetscStringize(func(__VA_ARGS__)),evsl_ierr_); \
  } while (0)

typedef struct {
  PetscBool         initialized;
  Mat               A;           /* problem matrix */
  Vec               x,y;         /* auxiliary vectors */
  PetscReal         *sli;        /* slice bounds */
  PetscInt          nev;         /* approximate number of wanted eigenvalues in each slice */
  PetscLayout       map;         /* used to distribute slices among MPI processes */
  PetscBool         estimrange;  /* the filter range was not set by the user */
  /* user parameters */
  PetscInt          nslices;     /* number of slices */
  PetscReal         lmin,lmax;   /* numerical range (min and max eigenvalue) */
  EPSEVSLDOSMethod  dos;         /* DOS method, either KPM or Lanczos */
  PetscInt          nvec;        /* number of sample vectors used for DOS */
  PetscInt          deg;         /* polynomial degree used for DOS (KPM only) */
  PetscInt          steps;       /* number of Lanczos steps used for DOS (Lanczos only) */
  PetscInt          npoints;     /* number of sample points used for DOS (Lanczos only) */
  PetscInt          max_deg;     /* maximum degree allowed for the polynomial */
  PetscReal         thresh;      /* threshold for accepting polynomial */
  EPSEVSLDamping    damping;     /* type of damping (for polynomial and for DOS-KPM) */
} EPS_EVSL;

static void AMatvec_EVSL(double *xa,double *ya,void *data)
{
  EPS_EVSL       *ctx = (EPS_EVSL*)data;
  Vec            x = ctx->x,y = ctx->y;
  Mat            A = ctx->A;

  PetscFunctionBegin;
  PetscCallAbort(PetscObjectComm((PetscObject)A),VecPlaceArray(x,(PetscScalar*)xa));
  PetscCallAbort(PetscObjectComm((PetscObject)A),VecPlaceArray(y,(PetscScalar*)ya));
  PetscCallAbort(PetscObjectComm((PetscObject)A),MatMult(A,x,y));
  PetscCallAbort(PetscObjectComm((PetscObject)A),VecResetArray(x));
  PetscCallAbort(PetscObjectComm((PetscObject)A),VecResetArray(y));
  PetscFunctionReturnVoid();
}

static PetscErrorCode EPSSetUp_EVSL(EPS eps)
{
  EPS_EVSL       *ctx = (EPS_EVSL*)eps->data;
  PetscMPIInt    size,rank;
  PetscBool      isshift;
  PetscScalar    *vinit;
  PetscReal      *mu,ecount,xintv[4],*xdos,*ydos;
  Vec            v0;
  Mat            A;
  PetscRandom    rnd;

  PetscFunctionBegin;
  EPSCheckStandard(eps);
  EPSCheckHermitian(eps);
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  PetscCheck(isshift,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support spectral transformations");

  if (ctx->initialized) EVSLFinish();
  EVSLStart();
  ctx->initialized=PETSC_TRUE;

  /* get number of slices per process */
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&rank));
  if (!ctx->nslices) ctx->nslices = size;
  PetscCall(PetscLayoutDestroy(&ctx->map));
  PetscCall(PetscLayoutCreateFromSizes(PetscObjectComm((PetscObject)eps),PETSC_DECIDE,ctx->nslices,1,&ctx->map));

  /* get matrix and prepare auxiliary vectors */
  PetscCall(MatDestroy(&ctx->A));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (size==1) {
    PetscCall(PetscObjectReference((PetscObject)A));
    ctx->A = A;
  } else PetscCall(MatCreateRedundantMatrix(A,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&ctx->A));
  SetAMatvec(eps->n,&AMatvec_EVSL,(void*)ctx);
  if (!ctx->x) PetscCall(MatCreateVecsEmpty(ctx->A,&ctx->x,&ctx->y));
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_STOPPING);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE);

  if (!eps->which) eps->which=EPS_ALL;
  PetscCheck(eps->which==EPS_ALL && eps->inta!=eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver requires setting an interval with EPSSetInterval()");

  /* estimate numerical range */
  if (ctx->estimrange || ctx->lmin == PETSC_MIN_REAL || ctx->lmax == PETSC_MAX_REAL) {
    PetscCall(MatCreateVecs(ctx->A,&v0,NULL));
    if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
    PetscCall(BVGetRandomContext(eps->V,&rnd));
    PetscCall(VecSetRandom(v0,rnd));
    PetscCall(VecGetArray(v0,&vinit));
    PetscCallEVSL(LanTrbounds,50,200,eps->tol,vinit,1,&ctx->lmin,&ctx->lmax,NULL);
    PetscCall(VecRestoreArray(v0,&vinit));
    PetscCall(VecDestroy(&v0));
    ctx->estimrange = PETSC_TRUE;   /* estimate if called again with another matrix */
  }
  PetscCheck(ctx->lmin<=eps->inta && ctx->lmax>=eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The requested interval [%g,%g] must be contained in the numerical range [%g,%g]",(double)eps->inta,(double)eps->intb,(double)ctx->lmin,(double)ctx->lmax);
  xintv[0] = eps->inta;
  xintv[1] = eps->intb;
  xintv[2] = ctx->lmin;
  xintv[3] = ctx->lmax;

  /* estimate number of eigenvalues in the interval */
  switch (ctx->dos) {
    case EPS_EVSL_DOS_KPM:
      PetscCall(PetscMalloc1(ctx->deg+1,&mu));
      if (!rank) PetscCallEVSL(kpmdos,ctx->deg,(int)ctx->damping,ctx->nvec,xintv,mu,&ecount);
      PetscCallMPI(MPI_Bcast(mu,ctx->deg+1,MPIU_REAL,0,PetscObjectComm((PetscObject)eps)));
      break;
    case EPS_EVSL_DOS_LANCZOS:
      PetscCall(PetscMalloc2(ctx->npoints,&xdos,ctx->npoints,&ydos));
      if (!rank) PetscCallEVSL(LanDos,ctx->nvec,PetscMin(ctx->steps,eps->n/2),ctx->npoints,xdos,ydos,&ecount,xintv);
      PetscCallMPI(MPI_Bcast(xdos,ctx->npoints,MPIU_REAL,0,PetscObjectComm((PetscObject)eps)));
      PetscCallMPI(MPI_Bcast(ydos,ctx->npoints,MPIU_REAL,0,PetscObjectComm((PetscObject)eps)));
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid DOS method");
  }
  PetscCallMPI(MPI_Bcast(&ecount,1,MPIU_REAL,0,PetscObjectComm((PetscObject)eps)));

  PetscCall(PetscInfo(eps,"Estimated eigenvalue count in the interval: %g\n",ecount));
  eps->ncv = (PetscInt)PetscCeilReal(1.5*ecount);

  /* slice the spectrum */
  PetscCall(PetscFree(ctx->sli));
  PetscCall(PetscMalloc1(ctx->nslices+1,&ctx->sli));
  if (ctx->dos == EPS_EVSL_DOS_KPM) {
    PetscCallEVSL(spslicer,ctx->sli,mu,ctx->deg,xintv,ctx->nslices,10*(PetscInt)ecount);
    PetscCall(PetscFree(mu));
  } else if (ctx->dos == EPS_EVSL_DOS_LANCZOS) {
    spslicer2(xdos,ydos,ctx->nslices,ctx->npoints,ctx->sli);
    PetscCall(PetscFree2(xdos,ydos));
  }

  /* approximate number of eigenvalues wanted in each slice */
  ctx->nev = (PetscInt)(1.0 + ecount/(PetscReal)ctx->nslices) + 2;

  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  PetscCall(EPSAllocateSolution(eps,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_EVSL(EPS eps)
{
  EPS_EVSL       *ctx = (EPS_EVSL*)eps->data;
  PetscInt       i,j,k=0,sl,mlan,nevout,*ind,nevmax,rstart,rend,*nevloc,*disp,N;
  PetscReal      *res,xintv[4],*errest;
  PetscScalar    *lam,*X,*Y,*vinit,*eigr;
  PetscMPIInt    size,rank;
  PetscRandom    rnd;
  Vec            v,w,v0,x;
  VecScatter     vs;
  IS             is;
  polparams      pol;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&rank));
  PetscCall(PetscLayoutGetRange(ctx->map,&rstart,&rend));
  nevmax = (rend-rstart)*ctx->nev;
  PetscCall(MatCreateVecs(ctx->A,&v0,NULL));
  PetscCall(BVGetRandomContext(eps->V,&rnd));
  PetscCall(VecSetRandom(v0,rnd));
  PetscCall(VecGetArray(v0,&vinit));
  PetscCall(PetscMalloc5(size,&nevloc,size+1,&disp,nevmax,&eigr,nevmax,&errest,nevmax*eps->n,&X));
  mlan = PetscMin(PetscMax(5*ctx->nev,300),eps->n);
  for (sl=rstart; sl<rend; sl++) {
    xintv[0] = ctx->sli[sl];
    xintv[1] = ctx->sli[sl+1];
    xintv[2] = ctx->lmin;
    xintv[3] = ctx->lmax;
    PetscCall(PetscInfo(ctx->A,"Subinterval %" PetscInt_FMT ": [%.4e, %.4e]\n",sl+1,xintv[0],xintv[1]));
    set_pol_def(&pol);
    pol.max_deg    = ctx->max_deg;
    pol.damping    = (int)ctx->damping;
    pol.thresh_int = ctx->thresh;
    find_pol(xintv,&pol);
    PetscCall(PetscInfo(ctx->A,"Polynomial [type = %" PetscInt_FMT "], deg %" PetscInt_FMT ", bar %e gam %e\n",pol.type,pol.deg,pol.bar,pol.gam));
    PetscCallEVSL(ChebLanNr,xintv,mlan,eps->tol,vinit,&pol,&nevout,&lam,&Y,&res,NULL);
    PetscCheck(k+nevout<=nevmax,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Too low estimation of eigenvalue count, try modifying the sampling parameters");
    free_pol(&pol);
    PetscCall(PetscInfo(ctx->A,"Computed %" PetscInt_FMT " eigenvalues\n",nevout));
    PetscCall(PetscMalloc1(nevout,&ind));
    sort_double(nevout,lam,ind);
    for (i=0;i<nevout;i++) {
      eigr[i+k]   = lam[i];
      errest[i+k] = res[ind[i]];
      PetscCall(PetscArraycpy(X+(i+k)*eps->n,Y+ind[i]*eps->n,eps->n));
    }
    k += nevout;
    if (lam) evsl_Free(lam);
    if (Y)   evsl_Free_device(Y);
    if (res) evsl_Free(res);
    PetscCall(PetscFree(ind));
  }
  PetscCall(VecRestoreArray(v0,&vinit));
  PetscCall(VecDestroy(&v0));

  /* gather eigenvalues computed by each MPI process */
  PetscCallMPI(MPI_Allgather(&k,1,MPIU_INT,nevloc,1,MPIU_INT,PetscObjectComm((PetscObject)eps)));
  eps->nev = nevloc[0];
  disp[0]  = 0;
  for (i=1;i<size;i++) {
    eps->nev += nevloc[i];
    disp[i]   = disp[i-1]+nevloc[i-1];
  }
  disp[size] = disp[size-1]+nevloc[size-1];
  PetscCheck(eps->nev<=eps->ncv,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Too low estimation of eigenvalue count, try modifying the sampling parameters");
  PetscCallMPI(MPI_Allgatherv(eigr,k,MPIU_SCALAR,eps->eigr,nevloc,disp,MPIU_SCALAR,PetscObjectComm((PetscObject)eps)));
  PetscCallMPI(MPI_Allgatherv(errest,k,MPIU_REAL,eps->errest,nevloc,disp,MPIU_REAL,PetscObjectComm((PetscObject)eps)));
  eps->nconv  = eps->nev;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;

  /* scatter computed eigenvectors and store them in eps->V */
  PetscCall(BVCreateVec(eps->V,&w));
  for (i=0;i<size;i++) {
    N = (rank==i)? eps->n: 0;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,N,&x));
    PetscCall(VecSetFromOptions(x));
    PetscCall(ISCreateStride(PETSC_COMM_SELF,N,0,1,&is));
    PetscCall(VecScatterCreate(x,is,w,is,&vs));
    PetscCall(ISDestroy(&is));
    for (j=disp[i];j<disp[i+1];j++) {
      PetscCall(BVGetColumn(eps->V,j,&v));
      if (rank==i) PetscCall(VecPlaceArray(x,X+(j-disp[i])*eps->n));
      PetscCall(VecScatterBegin(vs,x,v,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(vs,x,v,INSERT_VALUES,SCATTER_FORWARD));
      if (rank==i) PetscCall(VecResetArray(x));
      PetscCall(BVRestoreColumn(eps->V,j,&v));
    }
    PetscCall(VecScatterDestroy(&vs));
    PetscCall(VecDestroy(&x));
  }
  PetscCall(VecDestroy(&w));
  PetscCall(PetscFree5(nevloc,disp,eigr,errest,X));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLSetSlices_EVSL(EPS eps,PetscInt nslices)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (nslices == PETSC_DECIDE || nslices == PETSC_DEFAULT) nslices = 0;
  else PetscCheck(nslices>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Number of slices must be 1 at least");
  if (ctx->nslices != nslices) {
    ctx->nslices = nslices;
    eps->state   = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLSetSlices - Set the number of slices in which the interval must be
   subdivided.

   Logically Collective

   Input Parameters:
+  eps     - the eigensolver context
-  nslices - the number of slices

   Options Database Key:
.  -eps_evsl_slices <n> - set the number of slices to n

   Notes:
   By default, one slice per MPI process is used. Depending on the number of
   eigenvalues, using more slices may be beneficial, but very narrow subintervals
   imply higher polynomial degree.

   Level: intermediate

.seealso: EPSEVSLGetSlices()
@*/
PetscErrorCode EPSEVSLSetSlices(EPS eps,PetscInt nslices)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nslices,2);
  PetscTryMethod(eps,"EPSEVSLSetSlices_C",(EPS,PetscInt),(eps,nslices));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLGetSlices_EVSL(EPS eps,PetscInt *nslices)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  *nslices = ctx->nslices;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLGetSlices - Gets the number of slices in which the interval must be
   subdivided.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  nslices - the number of slices

   Level: intermediate

.seealso: EPSEVSLSetSlices()
@*/
PetscErrorCode EPSEVSLGetSlices(EPS eps,PetscInt *nslices)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(nslices,2);
  PetscUseMethod(eps,"EPSEVSLGetSlices_C",(EPS,PetscInt*),(eps,nslices));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLSetRange_EVSL(EPS eps,PetscReal lmin,PetscReal lmax)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  PetscCheck(lmin<lmax,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be lmin<lmax");
  if (ctx->lmin != lmin || ctx->lmax != lmax) {
    ctx->lmin  = lmin;
    ctx->lmax  = lmax;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLSetRange - Defines the numerical range (or field of values) of the problem,
   that is, the interval containing all eigenvalues.

   Logically Collective

   Input Parameters:
+  eps  - the eigensolver context
.  lmin - left end of the interval
-  lmax - right end of the interval

   Options Database Key:
.  -eps_evsl_range <a,b> - set [a,b] as the numerical range

   Notes:
   The filter will be most effective if the numerical range is tight, that is, lmin
   and lmax are good approximations to the leftmost and rightmost eigenvalues,
   respectively. If not set by the user, an approximation is computed internally.

   The wanted computational interval specified via EPSSetInterval() must be
   contained in the numerical range.

   Level: intermediate

.seealso: EPSEVSLGetRange(), EPSSetInterval()
@*/
PetscErrorCode EPSEVSLSetRange(EPS eps,PetscReal lmin,PetscReal lmax)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,lmin,2);
  PetscValidLogicalCollectiveReal(eps,lmax,3);
  PetscTryMethod(eps,"EPSEVSLSetRange_C",(EPS,PetscReal,PetscReal),(eps,lmin,lmax));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLGetRange_EVSL(EPS eps,PetscReal *lmin,PetscReal *lmax)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (lmin) *lmin = ctx->lmin;
  if (lmax) *lmax = ctx->lmax;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLGetRange - Gets the interval containing all eigenvalues.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  lmin - left end of the interval
-  lmax - right end of the interval

   Level: intermediate

.seealso: EPSEVSLSetRange()
@*/
PetscErrorCode EPSEVSLGetRange(EPS eps,PetscReal *lmin,PetscReal *lmax)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSEVSLGetRange_C",(EPS,PetscReal*,PetscReal*),(eps,lmin,lmax));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLSetDOSParameters_EVSL(EPS eps,EPSEVSLDOSMethod dos,PetscInt nvec,PetscInt deg,PetscInt steps,PetscInt npoints)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  ctx->dos = dos;
  if (nvec == PETSC_DETERMINE) ctx->nvec = 80;
  else if (nvec != PETSC_CURRENT) {
    PetscCheck(nvec>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The nvec argument must be > 0");
    ctx->nvec = nvec;
  }
  switch (dos) {
    case EPS_EVSL_DOS_KPM:
      if (deg == PETSC_DETERMINE) ctx->deg = 300;
      else if (deg != PETSC_CURRENT) {
        PetscCheck(deg>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The deg argument must be > 0");
        ctx->deg = deg;
      }
      break;
    case EPS_EVSL_DOS_LANCZOS:
      if (steps == PETSC_DETERMINE) ctx->steps = 40;
      else if (steps != PETSC_CURRENT) {
        PetscCheck(steps>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The steps argument must be > 0");
        ctx->steps = steps;
      }
      if (npoints == PETSC_DETERMINE) ctx->npoints = 200;
      else if (npoints != PETSC_CURRENT) {
        PetscCheck(npoints>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The npoints argument must be > 0");
        ctx->npoints = npoints;
      }
      break;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLSetDOSParameters - Defines the parameters used for computing the
   density of states (DOS) in the EVSL solver.

   Logically Collective

   Input Parameters:
+  eps     - the eigensolver context
.  dos     - DOS method, either KPM or Lanczos
.  nvec    - number of sample vectors
.  deg     - polynomial degree (KPM only)
.  steps   - number of Lanczos steps (Lanczos only)
-  npoints - number of sample points (Lanczos only)

   Options Database Keys:
+  -eps_evsl_dos_method <dos> - set the DOS method, either kpm or lanczos
.  -eps_evsl_dos_nvec <n> - set the number of sample vectors
.  -eps_evsl_dos_degree <n> - set the polynomial degree
.  -eps_evsl_dos_steps <n> - set the number of Lanczos steps
-  -eps_evsl_dos_npoints <n> - set the number of sample points

   Notes:
   The density of states (or spectral density) can be approximated with two
   methods, kernel polynomial method (KPM) or Lanczos. Some parameters for
   these methods can be set by the user with this function, with some of
   them being relevant for one of the methods only.

   For the integer argumens, you can use PETSC_CURRENT to keep the current
   value, and PETSC_DETERMINE to set them to a reasonable default.

   Level: intermediate

.seealso: EPSEVSLGetDOSParameters()
@*/
PetscErrorCode EPSEVSLSetDOSParameters(EPS eps,EPSEVSLDOSMethod dos,PetscInt nvec,PetscInt deg,PetscInt steps,PetscInt npoints)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,dos,2);
  PetscValidLogicalCollectiveInt(eps,nvec,3);
  PetscValidLogicalCollectiveInt(eps,deg,4);
  PetscValidLogicalCollectiveInt(eps,steps,5);
  PetscValidLogicalCollectiveInt(eps,npoints,6);
  PetscTryMethod(eps,"EPSEVSLSetDOSParameters_C",(EPS,EPSEVSLDOSMethod,PetscInt,PetscInt,PetscInt,PetscInt),(eps,dos,nvec,deg,steps,npoints));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLGetDOSParameters_EVSL(EPS eps,EPSEVSLDOSMethod *dos,PetscInt *nvec,PetscInt *deg,PetscInt *steps,PetscInt *npoints)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (dos)     *dos     = ctx->dos;
  if (nvec)    *nvec    = ctx->nvec;
  if (deg)     *deg     = ctx->deg;
  if (steps)   *steps   = ctx->steps;
  if (npoints) *npoints = ctx->npoints;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLGetDOSParameters - Gets the parameters used for computing the
   density of states (DOS) in the EVSL solver.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  dos     - DOS method, either KPM or Lanczos
.  nvec    - number of sample vectors
.  deg     - polynomial degree (KPM only)
.  steps   - number of Lanczos steps (Lanczos only)
-  npoints - number of sample points (Lanczos only)

   Level: intermediate

.seealso: EPSEVSLSetDOSParameters()
@*/
PetscErrorCode EPSEVSLGetDOSParameters(EPS eps,EPSEVSLDOSMethod *dos,PetscInt *nvec,PetscInt *deg,PetscInt *steps,PetscInt *npoints)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSEVSLGetDOSParameters_C",(EPS,EPSEVSLDOSMethod*,PetscInt*,PetscInt*,PetscInt*,PetscInt*),(eps,dos,nvec,deg,steps,npoints));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLSetPolParameters_EVSL(EPS eps,PetscInt max_deg,PetscReal thresh)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (max_deg == PETSC_DETERMINE) ctx->max_deg = 10000;
  else if (max_deg != PETSC_CURRENT) {
    PetscCheck(max_deg>2,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The max_deg argument must be > 2");
    ctx->max_deg = max_deg;
  }
  if (thresh == (PetscReal)PETSC_DETERMINE) ctx->thresh = 0.8;
  else if (thresh != (PetscReal)PETSC_CURRENT) {
    PetscCheck(thresh>0.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The thresh argument must be > 0.0");
    ctx->thresh = thresh;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLSetPolParameters - Defines the parameters used for building the
   building the polynomial in the EVSL solver.

   Logically Collective

   Input Parameters:
+  eps     - the eigensolver context
.  max_deg - maximum degree allowed for the polynomial
-  thresh  - threshold for accepting polynomial

   Options Database Keys:
+  -eps_evsl_pol_max_deg <d> - set maximum polynomial degree
-  -eps_evsl_pol_thresh <t> - set the threshold

   Note:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Level: intermediate

.seealso: EPSEVSLGetPolParameters()
@*/
PetscErrorCode EPSEVSLSetPolParameters(EPS eps,PetscInt max_deg,PetscReal thresh)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,max_deg,2);
  PetscValidLogicalCollectiveReal(eps,thresh,3);
  PetscTryMethod(eps,"EPSEVSLSetPolParameters_C",(EPS,PetscInt,PetscReal),(eps,max_deg,thresh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLGetPolParameters_EVSL(EPS eps,PetscInt *max_deg,PetscReal *thresh)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (max_deg) *max_deg = ctx->max_deg;
  if (thresh)  *thresh  = ctx->thresh;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLGetPolParameters - Gets the parameters used for building the
   polynomial in the EVSL solver.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  max_deg - the maximum degree of the polynomial
-  thresh  - the threshold

   Level: intermediate

.seealso: EPSEVSLSetPolParameters()
@*/
PetscErrorCode EPSEVSLGetPolParameters(EPS eps,PetscInt *max_deg,PetscReal *thresh)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSEVSLGetPolParameters_C",(EPS,PetscInt*,PetscReal*),(eps,max_deg,thresh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLSetDamping_EVSL(EPS eps,EPSEVSLDamping damping)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (ctx->damping != damping) {
    ctx->damping = damping;
    eps->state   = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLSetDamping - Set the type of damping to be used in EVSL.

   Logically Collective

   Input Parameters:
+  eps     - the eigensolver context
-  damping - the type of damping

   Options Database Key:
.  -eps_evsl_damping <n> - set the type of damping

   Notes:
   Damping is applied when building the polynomial to be used when solving the
   eigenproblem, and also during estimation of DOS with the KPM method.

   Level: intermediate

.seealso: EPSEVSLGetDamping(), EPSEVSLSetDOSParameters()
@*/
PetscErrorCode EPSEVSLSetDamping(EPS eps,EPSEVSLDamping damping)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,damping,2);
  PetscTryMethod(eps,"EPSEVSLSetDamping_C",(EPS,EPSEVSLDamping),(eps,damping));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSEVSLGetDamping_EVSL(EPS eps,EPSEVSLDamping *damping)
{
  EPS_EVSL *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  *damping = ctx->damping;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSEVSLGetDamping - Gets the type of damping.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  damping - the type of damping

   Level: intermediate

.seealso: EPSEVSLSetDamping()
@*/
PetscErrorCode EPSEVSLGetDamping(EPS eps,EPSEVSLDamping *damping)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(damping,2);
  PetscUseMethod(eps,"EPSEVSLGetDamping_C",(EPS,EPSEVSLDamping*),(eps,damping));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_EVSL(EPS eps,PetscViewer viewer)
{
  PetscBool      isascii;
  EPS_EVSL       *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  numerical range = [%g,%g]\n",(double)ctx->lmin,(double)ctx->lmax));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  number of slices = %" PetscInt_FMT "\n",ctx->nslices));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  type of damping = %s\n",EPSEVSLDampings[ctx->damping]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  computing DOS with %s: nvec=%" PetscInt_FMT ", ",EPSEVSLDOSMethods[ctx->dos],ctx->nvec));
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_FALSE));
    switch (ctx->dos) {
      case EPS_EVSL_DOS_KPM:
        PetscCall(PetscViewerASCIIPrintf(viewer,"degree=%" PetscInt_FMT "\n",ctx->deg));
        break;
      case EPS_EVSL_DOS_LANCZOS:
        PetscCall(PetscViewerASCIIPrintf(viewer,"steps=%" PetscInt_FMT ", npoints=%" PetscInt_FMT "\n",ctx->steps,ctx->npoints));
        break;
    }
    PetscCall(PetscViewerASCIIUseTabs(viewer,PETSC_TRUE));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  polynomial parameters: max degree = %" PetscInt_FMT ", threshold = %g\n",ctx->max_deg,(double)ctx->thresh));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_EVSL(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscReal        array[2]={0,0},th;
  PetscInt         k,i1,i2,i3,i4;
  PetscBool        flg,flg1;
  EPSEVSLDOSMethod dos;
  EPSEVSLDamping   damping;
  EPS_EVSL         *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS EVSL Options");

    k = 2;
    PetscCall(PetscOptionsRealArray("-eps_evsl_range","Interval containing all eigenvalues (two real values separated with a comma without spaces)","EPSEVSLSetRange",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_SIZ,"Must pass two values in -eps_evsl_range (comma-separated without spaces)");
      PetscCall(EPSEVSLSetRange(eps,array[0],array[1]));
    }

    PetscCall(PetscOptionsInt("-eps_evsl_slices","Number of slices","EPSEVSLSetSlices",ctx->nslices,&i1,&flg));
    if (flg) PetscCall(EPSEVSLSetSlices(eps,i1));

    PetscCall(PetscOptionsEnum("-eps_evsl_damping","Type of damping","EPSEVSLSetDamping",EPSEVSLDampings,(PetscEnum)ctx->damping,(PetscEnum*)&damping,&flg));
    if (flg) PetscCall(EPSEVSLSetDamping(eps,damping));

    PetscCall(EPSEVSLGetDOSParameters(eps,&dos,&i1,&i2,&i3,&i4));
    PetscCall(PetscOptionsEnum("-eps_evsl_dos_method","Method to compute the DOS","EPSEVSLSetDOSParameters",EPSEVSLDOSMethods,(PetscEnum)ctx->dos,(PetscEnum*)&dos,&flg));
    PetscCall(PetscOptionsInt("-eps_evsl_dos_nvec","Number of sample vectors for DOS","EPSEVSLSetDOSParameters",i1,&i1,&flg1));
    if (flg1) flg = PETSC_TRUE;
    PetscCall(PetscOptionsInt("-eps_evsl_dos_degree","Polynomial degree used for DOS","EPSEVSLSetDOSParameters",i2,&i2,&flg1));
    if (flg1) flg = PETSC_TRUE;
    PetscCall(PetscOptionsInt("-eps_evsl_dos_steps","Number of Lanczos steps in DOS","EPSEVSLSetDOSParameters",i3,&i3,&flg1));
    if (flg1) flg = PETSC_TRUE;
    PetscCall(PetscOptionsInt("-eps_evsl_dos_npoints","Number of sample points used for DOS","EPSEVSLSetDOSParameters",i4,&i4,&flg1));
    if (flg || flg1) PetscCall(EPSEVSLSetDOSParameters(eps,dos,i1,i2,i3,i4));

    PetscCall(EPSEVSLGetPolParameters(eps,&i1,&th));
    PetscCall(PetscOptionsInt("-eps_evsl_pol_max_deg","Maximum degree allowed for the polynomial","EPSEVSLSetPolParameters",i1,&i1,&flg));
    PetscCall(PetscOptionsReal("-eps_evsl_pol_threshold","Threshold for accepting polynomial","EPSEVSLSetPolParameters",th,&th,&flg1));
    if (flg || flg1) PetscCall(EPSEVSLSetPolParameters(eps,i1,th));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_EVSL(EPS eps)
{
  EPS_EVSL       *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  if (ctx->initialized) EVSLFinish();
  PetscCall(PetscLayoutDestroy(&ctx->map));
  PetscCall(PetscFree(ctx->sli));
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetRange_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetRange_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetSlices_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetSlices_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetDOSParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetDOSParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetPolParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetPolParameters_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetDamping_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetDamping_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_EVSL(EPS eps)
{
  EPS_EVSL       *ctx = (EPS_EVSL*)eps->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->A));
  PetscCall(VecDestroy(&ctx->x));
  PetscCall(VecDestroy(&ctx->y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_EVSL(EPS eps)
{
  EPS_EVSL       *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  ctx->nslices = 0;
  ctx->lmin    = PETSC_MIN_REAL;
  ctx->lmax    = PETSC_MAX_REAL;
  ctx->dos     = EPS_EVSL_DOS_KPM;
  ctx->nvec    = 80;
  ctx->deg     = 300;
  ctx->steps   = 40;
  ctx->npoints = 200;
  ctx->max_deg = 10000;
  ctx->thresh  = 0.8;
  ctx->damping = EPS_EVSL_DAMPING_SIGMA;

  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_EVSL;
  eps->ops->setup          = EPSSetUp_EVSL;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->setfromoptions = EPSSetFromOptions_EVSL;
  eps->ops->destroy        = EPSDestroy_EVSL;
  eps->ops->reset          = EPSReset_EVSL;
  eps->ops->view           = EPSView_EVSL;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_NoFactor;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetRange_C",EPSEVSLSetRange_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetRange_C",EPSEVSLGetRange_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetSlices_C",EPSEVSLSetSlices_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetSlices_C",EPSEVSLGetSlices_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetDOSParameters_C",EPSEVSLSetDOSParameters_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetDOSParameters_C",EPSEVSLGetDOSParameters_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetPolParameters_C",EPSEVSLSetPolParameters_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetPolParameters_C",EPSEVSLGetPolParameters_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLSetDamping_C",EPSEVSLSetDamping_EVSL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSEVSLGetDamping_C",EPSEVSLGetDamping_EVSL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
