/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "krylovschur"

   Method: Krylov-Schur with spectrum slicing for symmetric eigenproblems

   References:

       [1] R.G. Grimes et al., "A shifted block Lanczos algorithm for
           solving sparse symmetric generalized eigenproblems", SIAM J.
           Matrix Anal. Appl. 15(1):228-272, 1994.

       [2] C. Campos and J.E. Roman, "Spectrum slicing strategies based
           on restarted Lanczos methods", Numer. Algor. 60(2):279-295,
           2012.
*/

#include <slepc/private/epsimpl.h>
#include "krylovschur.h"

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-slice,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Strategies for spectrum slicing based on restarted {Lanczos} methods\",\n"
  "   journal = \"Numer. Algorithms\",\n"
  "   volume = \"60\",\n"
  "   number = \"2\",\n"
  "   pages = \"279--295\",\n"
  "   year = \"2012,\"\n"
  "   doi = \"https://doi.org/10.1007/s11075-012-9564-z\"\n"
  "}\n";

#define SLICE_PTOL PETSC_SQRT_MACHINE_EPSILON

static PetscErrorCode EPSSliceResetSR(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr=ctx->sr;
  EPS_shift       s;

  PetscFunctionBegin;
  if (sr) {
    if (ctx->npart>1) {
      PetscCall(BVDestroy(&sr->V));
      PetscCall(PetscFree4(sr->eigr,sr->eigi,sr->errest,sr->perm));
    }
    /* Reviewing list of shifts to free memory */
    s = sr->s0;
    if (s) {
      while (s->neighb[1]) {
        s = s->neighb[1];
        PetscCall(PetscFree(s->neighb[0]));
      }
      PetscCall(PetscFree(s));
    }
    PetscCall(PetscFree(sr));
  }
  ctx->sr = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSReset_KrylovSchur_Slice(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (!ctx->global) PetscFunctionReturn(PETSC_SUCCESS);
  /* Reset auxiliary EPS */
  PetscCall(EPSSliceResetSR(ctx->eps));
  PetscCall(EPSReset(ctx->eps));
  PetscCall(EPSSliceResetSR(eps));
  PetscCall(PetscFree(ctx->inertias));
  PetscCall(PetscFree(ctx->shifts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSDestroy_KrylovSchur_Slice(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (!ctx->global) PetscFunctionReturn(PETSC_SUCCESS);
  /* Destroy auxiliary EPS */
  PetscCall(EPSReset_KrylovSchur_Slice(eps));
  PetscCall(EPSDestroy(&ctx->eps));
  if (ctx->npart>1) {
    PetscCall(PetscSubcommDestroy(&ctx->subc));
    if (ctx->commset) {
      PetscCallMPI(MPI_Comm_free(&ctx->commrank));
      ctx->commset = PETSC_FALSE;
    }
    PetscCall(ISDestroy(&ctx->isrow));
    PetscCall(ISDestroy(&ctx->iscol));
    PetscCall(MatDestroyMatrices(1,&ctx->submata));
    PetscCall(MatDestroyMatrices(1,&ctx->submatb));
  }
  PetscCall(PetscFree(ctx->subintervals));
  PetscCall(PetscFree(ctx->nconv_loc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSSliceAllocateSolution - Allocate memory storage for common variables such
  as eigenvalues and eigenvectors. The argument extra is used for methods
  that require a working basis slightly larger than ncv.
*/
static PetscErrorCode EPSSliceAllocateSolution(EPS eps,PetscInt extra)
{
  EPS_KRYLOVSCHUR    *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscReal          eta;
  PetscInt           k;
  BVType             type;
  BVOrthogType       orthog_type;
  BVOrthogRefineType orthog_ref;
  BVOrthogBlockType  ob_type;
  Mat                matrix;
  Vec                t;
  EPS_SR             sr = ctx->sr;

  PetscFunctionBegin;
  /* allocate space for eigenvalues and friends */
  k = PetscMax(1,sr->numEigs);
  PetscCall(PetscFree4(sr->eigr,sr->eigi,sr->errest,sr->perm));
  PetscCall(PetscMalloc4(k,&sr->eigr,k,&sr->eigi,k,&sr->errest,k,&sr->perm));

  /* allocate sr->V and transfer options from eps->V */
  PetscCall(BVDestroy(&sr->V));
  PetscCall(BVCreate(PetscObjectComm((PetscObject)eps),&sr->V));
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  if (!((PetscObject)eps->V)->type_name) PetscCall(BVSetType(sr->V,BVMAT));
  else {
    PetscCall(BVGetType(eps->V,&type));
    PetscCall(BVSetType(sr->V,type));
  }
  PetscCall(STMatCreateVecsEmpty(eps->st,&t,NULL));
  PetscCall(BVSetSizesFromVec(sr->V,t,k));
  PetscCall(VecDestroy(&t));
  PetscCall(EPS_SetInnerProduct(eps));
  PetscCall(BVGetMatrix(eps->V,&matrix,NULL));
  PetscCall(BVSetMatrix(sr->V,matrix,PETSC_FALSE));
  PetscCall(BVGetOrthogonalization(eps->V,&orthog_type,&orthog_ref,&eta,&ob_type));
  PetscCall(BVSetOrthogonalization(sr->V,orthog_type,orthog_ref,eta,ob_type));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSliceGetEPS(EPS eps)
{
  EPS_KRYLOVSCHUR    *ctx=(EPS_KRYLOVSCHUR*)eps->data,*ctx_local;
  BV                 V;
  BVType             type;
  PetscReal          eta;
  BVOrthogType       orthog_type;
  BVOrthogRefineType orthog_ref;
  BVOrthogBlockType  ob_type;
  PetscInt           i;
  PetscReal          h,a,b;
  PetscRandom        rand;
  EPS_SR             sr=ctx->sr;

  PetscFunctionBegin;
  if (!ctx->eps) PetscCall(EPSKrylovSchurGetChildEPS(eps,&ctx->eps));

  /* Determine subintervals */
  if (ctx->npart==1) {
    a = eps->inta; b = eps->intb;
  } else {
    if (!ctx->subintset) { /* uniform distribution if no set by user */
      PetscCheck(sr->hasEnd,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Global interval must be bounded for splitting it in uniform subintervals");
      h = (eps->intb-eps->inta)/ctx->npart;
      a = eps->inta+ctx->subc->color*h;
      b = (ctx->subc->color==ctx->npart-1)?eps->intb:eps->inta+(ctx->subc->color+1)*h;
      PetscCall(PetscFree(ctx->subintervals));
      PetscCall(PetscMalloc1(ctx->npart+1,&ctx->subintervals));
      for (i=0;i<ctx->npart;i++) ctx->subintervals[i] = eps->inta+h*i;
      ctx->subintervals[ctx->npart] = eps->intb;
    } else {
      a = ctx->subintervals[ctx->subc->color];
      b = ctx->subintervals[ctx->subc->color+1];
    }
  }
  PetscCall(EPSSetInterval(ctx->eps,a,b));
  PetscCall(EPSSetConvergenceTest(ctx->eps,eps->conv));
  PetscCall(EPSSetDimensions(ctx->eps,ctx->nev,ctx->ncv,ctx->mpd));
  PetscCall(EPSKrylovSchurSetLocking(ctx->eps,ctx->lock));

  ctx_local = (EPS_KRYLOVSCHUR*)ctx->eps->data;
  ctx_local->detect = ctx->detect;

  /* transfer options from eps->V */
  PetscCall(EPSGetBV(ctx->eps,&V));
  PetscCall(BVGetRandomContext(V,&rand));  /* make sure the random context is available when duplicating */
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  if (!((PetscObject)eps->V)->type_name) PetscCall(BVSetType(V,BVMAT));
  else {
    PetscCall(BVGetType(eps->V,&type));
    PetscCall(BVSetType(V,type));
  }
  PetscCall(BVGetOrthogonalization(eps->V,&orthog_type,&orthog_ref,&eta,&ob_type));
  PetscCall(BVSetOrthogonalization(V,orthog_type,orthog_ref,eta,ob_type));

  ctx->eps->which = eps->which;
  ctx->eps->max_it = eps->max_it;
  ctx->eps->tol = eps->tol;
  ctx->eps->purify = eps->purify;
  if (eps->tol==(PetscReal)PETSC_DETERMINE) eps->tol = SLEPC_DEFAULT_TOL;
  PetscCall(EPSSetProblemType(ctx->eps,eps->problem_type));
  PetscCall(EPSSetUp(ctx->eps));
  ctx->eps->nconv = 0;
  ctx->eps->its   = 0;
  for (i=0;i<ctx->eps->ncv;i++) {
    ctx->eps->eigr[i]   = 0.0;
    ctx->eps->eigi[i]   = 0.0;
    ctx->eps->errest[i] = 0.0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSliceGetInertia(EPS eps,PetscReal shift,PetscInt *inertia,PetscInt *zeros)
{
  KSP            ksp,kspr;
  PC             pc;
  Mat            F;
  PetscReal      nzshift=shift;
  PetscBool      flg;

  PetscFunctionBegin;
  if (shift >= PETSC_MAX_REAL) { /* Right-open interval */
    if (inertia) *inertia = eps->n;
  } else if (shift <= PETSC_MIN_REAL) {
    if (inertia) *inertia = 0;
    if (zeros) *zeros = 0;
  } else {
    /* If the shift is zero, perturb it to a very small positive value.
       The goal is that the nonzero pattern is the same in all cases and reuse
       the symbolic factorizations */
    nzshift = (shift==0.0)? 10.0/PETSC_MAX_REAL: shift;
    PetscCall(STSetShift(eps->st,nzshift));
    PetscCall(STGetKSP(eps->st,&ksp));
    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCREDUNDANT,&flg));
    if (flg) {
      PetscCall(PCRedundantGetKSP(pc,&kspr));
      PetscCall(KSPGetPC(kspr,&pc));
    }
    PetscCall(PCFactorGetMatrix(pc,&F));
    PetscCall(MatGetInertia(F,inertia,zeros,NULL));
  }
  if (inertia) PetscCall(PetscInfo(eps,"Computed inertia at shift %g: %" PetscInt_FMT "\n",(double)nzshift,*inertia));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Dummy backtransform operation
 */
static PetscErrorCode EPSBackTransform_Skip(EPS eps)
{
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSetUp_KrylovSchur_Slice(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data,*ctx_glob;
  EPS_SR          sr,sr_loc,sr_glob;
  PetscInt        nEigs,dssz=1,i,zeros=0,off=0,method,hiteig=0;
  PetscMPIInt     nproc,rank=0,aux;
  PetscReal       r;
  MPI_Request     req;
  Mat             A,B=NULL;
  DSParallelType  ptype;
  MPI_Comm        child;

  PetscFunctionBegin;
  if (eps->nev==0) eps->nev = 1;
  if (ctx->global) {
    EPSCheckHermitianDefiniteCondition(eps,PETSC_TRUE," with spectrum slicing");
    EPSCheckSinvertCayleyCondition(eps,PETSC_TRUE," with spectrum slicing");
    PetscCheck(eps->inta!=eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues unless you provide a computational interval with EPSSetInterval()");
    PetscCheck(eps->intb<PETSC_MAX_REAL || eps->inta>PETSC_MIN_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The defined computational interval should have at least one of their sides bounded");
    PetscCheck(eps->nds==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Spectrum slicing not supported in combination with deflation space");
    EPSCheckUnsupportedCondition(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_STOPPING,PETSC_TRUE," with spectrum slicing");
    EPSCheckIgnoredCondition(eps,EPS_FEATURE_BALANCE,PETSC_TRUE," with spectrum slicing");
    if (eps->tol==(PetscReal)PETSC_DETERMINE) {
#if defined(PETSC_USE_REAL_SINGLE)
      eps->tol = SLEPC_DEFAULT_TOL;
#else
      /* use tighter tolerance */
      eps->tol = SLEPC_DEFAULT_TOL*1e-2;
#endif
    }
    if (eps->max_it==PETSC_DETERMINE) eps->max_it = 100;
    if (ctx->nev==1) ctx->nev = PetscMin(40,eps->n);  /* nev not set, use default value */
    PetscCheck(eps->n<=10 || ctx->nev>=10,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"nev cannot be less than 10 in spectrum slicing runs");
  }
  eps->ops->backtransform = EPSBackTransform_Skip;

  /* create spectrum slicing context and initialize it */
  PetscCall(EPSSliceResetSR(eps));
  PetscCall(PetscNew(&sr));
  ctx->sr = sr;
  sr->itsKs = 0;
  sr->nleap = 0;
  sr->nMAXCompl = eps->nev/4;
  sr->iterCompl = eps->max_it/4;
  sr->sPres = NULL;
  sr->nS = 0;

  if (ctx->npart==1 || ctx->global) {
    /* check presence of ends and finding direction */
    if ((eps->inta > PETSC_MIN_REAL && !(ctx->subintervals && ctx->subintervals[0]==ctx->subintervals[1])) || eps->intb >= PETSC_MAX_REAL) {
      sr->int0 = eps->inta;
      sr->int1 = eps->intb;
      sr->dir = 1;
      if (eps->intb >= PETSC_MAX_REAL) { /* Right-open interval */
        sr->hasEnd = PETSC_FALSE;
      } else sr->hasEnd = PETSC_TRUE;
    } else {
      sr->int0 = eps->intb;
      sr->int1 = eps->inta;
      sr->dir = -1;
      sr->hasEnd = PetscNot(eps->inta <= PETSC_MIN_REAL);
    }
  }
  if (ctx->global) {
    PetscCall(EPSSetDimensions_Default(eps,&ctx->nev,&ctx->ncv,&ctx->mpd));
    /* create subintervals and initialize auxiliary eps for slicing runs */
    PetscCall(EPSKrylovSchurGetChildEPS(eps,&ctx->eps));
    /* prevent computation of factorization in global eps */
    PetscCall(STSetTransform(eps->st,PETSC_FALSE));
    PetscCall(EPSSliceGetEPS(eps));
    sr_loc = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;
    if (ctx->npart>1) {
      PetscCall(PetscSubcommGetChild(ctx->subc,&child));
      if ((sr->dir>0&&ctx->subc->color==0)||(sr->dir<0&&ctx->subc->color==ctx->npart-1)) sr->inertia0 = sr_loc->inertia0;
      PetscCallMPI(MPI_Comm_rank(child,&rank));
      if (!rank) {
        PetscCall(PetscMPIIntCast((sr->dir>0)?0:ctx->npart-1,&aux));
        PetscCallMPI(MPI_Bcast(&sr->inertia0,1,MPIU_INT,aux,ctx->commrank));
      }
      PetscCallMPI(MPI_Bcast(&sr->inertia0,1,MPIU_INT,0,child));
      PetscCall(PetscFree(ctx->nconv_loc));
      PetscCall(PetscMalloc1(ctx->npart,&ctx->nconv_loc));
      PetscCallMPI(MPI_Comm_size(((PetscObject)eps)->comm,&nproc));
      if (sr->dir<0) off = 1;
      if (nproc%ctx->npart==0) { /* subcommunicators with the same size */
        PetscCall(PetscMPIIntCast(sr_loc->numEigs,&aux));
        PetscCallMPI(MPI_Allgather(&aux,1,MPI_INT,ctx->nconv_loc,1,MPI_INT,ctx->commrank));
        PetscCallMPI(MPI_Allgather(sr_loc->dir==sr->dir?&sr_loc->int0:&sr_loc->int1,1,MPIU_REAL,ctx->subintervals+off,1,MPIU_REAL,ctx->commrank));
      } else {
        PetscCallMPI(MPI_Comm_rank(child,&rank));
        if (!rank) {
          PetscCall(PetscMPIIntCast(sr_loc->numEigs,&aux));
          PetscCallMPI(MPI_Allgather(&aux,1,MPI_INT,ctx->nconv_loc,1,MPI_INT,ctx->commrank));
          PetscCallMPI(MPI_Allgather(sr_loc->dir==sr->dir?&sr_loc->int0:&sr_loc->int1,1,MPIU_REAL,ctx->subintervals+off,1,MPIU_REAL,ctx->commrank));
        }
        PetscCall(PetscMPIIntCast(ctx->npart,&aux));
        PetscCallMPI(MPI_Bcast(ctx->nconv_loc,aux,MPI_INT,0,child));
        PetscCallMPI(MPI_Bcast(ctx->subintervals+off,aux,MPIU_REAL,0,child));
      }
      nEigs = 0;
      for (i=0;i<ctx->npart;i++) nEigs += ctx->nconv_loc[i];
    } else {
      nEigs = sr_loc->numEigs;
      sr->inertia0 = sr_loc->inertia0;
      sr->dir = sr_loc->dir;
    }
    sr->inertia1 = sr->inertia0+sr->dir*nEigs;
    sr->numEigs = nEigs;
    eps->nev = nEigs;
    eps->ncv = nEigs;
    eps->mpd = nEigs;
  } else {
    ctx_glob = (EPS_KRYLOVSCHUR*)ctx->eps->data;
    sr_glob = ctx_glob->sr;
    if (ctx->npart>1) {
      sr->dir = sr_glob->dir;
      sr->int0 = (sr->dir==1)?eps->inta:eps->intb;
      sr->int1 = (sr->dir==1)?eps->intb:eps->inta;
      if ((sr->dir>0&&ctx->subc->color==ctx->npart-1)||(sr->dir<0&&ctx->subc->color==0)) sr->hasEnd = sr_glob->hasEnd;
      else sr->hasEnd = PETSC_TRUE;
    }
    /* sets first shift */
    PetscCall(STSetShift(eps->st,(sr->int0==0.0)?10.0/PETSC_MAX_REAL:sr->int0));
    PetscCall(STSetUp(eps->st));

    /* compute inertia0 */
    PetscCall(EPSSliceGetInertia(eps,sr->int0,&sr->inertia0,ctx->detect?&zeros:NULL));
    /* undocumented option to control what to do when an eigenvalue is found:
       - error out if it's the endpoint of the user-provided interval (or sub-interval)
       - if it's an endpoint computed internally:
          + if hiteig=0 error out
          + else if hiteig=1 the subgroup that hit the eigenvalue does nothing
          + otherwise the subgroup that hit the eigenvalue perturbs the shift and recomputes inertia
    */
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-eps_krylovschur_hiteigenvalue",&hiteig,NULL));
    if (zeros) { /* error in factorization */
      PetscCheck(sr->int0!=ctx->eps->inta && sr->int0!=ctx->eps->intb,((PetscObject)eps)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in the interval endpoint");
      PetscCheck(!ctx_glob->subintset || hiteig,((PetscObject)eps)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in an interval endpoint defined by user");
      if (hiteig==1) { /* idle subgroup */
        sr->inertia0 = -1;
      } else { /* perturb shift */
        sr->int0 *= (1.0+SLICE_PTOL);
        PetscCall(EPSSliceGetInertia(eps,sr->int0,&sr->inertia0,&zeros));
        PetscCheck(zeros==0,((PetscObject)eps)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",(double)sr->int1);
      }
    }
    if (ctx->npart>1) {
      PetscCall(PetscSubcommGetChild(ctx->subc,&child));
      /* inertia1 is received from neighbour */
      PetscCallMPI(MPI_Comm_rank(child,&rank));
      if (!rank) {
        if (sr->inertia0!=-1 && ((sr->dir>0 && ctx->subc->color>0) || (sr->dir<0 && ctx->subc->color<ctx->npart-1))) { /* send inertia0 to neighbour0 */
          PetscCall(PetscMPIIntCast(ctx->subc->color-sr->dir,&aux));
          PetscCallMPI(MPI_Isend(&sr->inertia0,1,MPIU_INT,aux,0,ctx->commrank,&req));
          PetscCallMPI(MPI_Isend(&sr->int0,1,MPIU_REAL,aux,0,ctx->commrank,&req));
        }
        if ((sr->dir>0 && ctx->subc->color<ctx->npart-1)|| (sr->dir<0 && ctx->subc->color>0)) { /* receive inertia1 from neighbour1 */
          PetscCall(PetscMPIIntCast(ctx->subc->color+sr->dir,&aux));
          PetscCallMPI(MPI_Recv(&sr->inertia1,1,MPIU_INT,aux,0,ctx->commrank,MPI_STATUS_IGNORE));
          PetscCallMPI(MPI_Recv(&sr->int1,1,MPIU_REAL,aux,0,ctx->commrank,MPI_STATUS_IGNORE));
        }
        if (sr->inertia0==-1 && !(sr->dir>0 && ctx->subc->color==ctx->npart-1) && !(sr->dir<0 && ctx->subc->color==0)) {
          sr->inertia0 = sr->inertia1; sr->int0 = sr->int1;
          PetscCall(PetscMPIIntCast(ctx->subc->color-sr->dir,&aux));
          PetscCallMPI(MPI_Isend(&sr->inertia0,1,MPIU_INT,aux,0,ctx->commrank,&req));
          PetscCallMPI(MPI_Isend(&sr->int0,1,MPIU_REAL,aux,0,ctx->commrank,&req));
        }
      }
      if ((sr->dir>0 && ctx->subc->color<ctx->npart-1)||(sr->dir<0 && ctx->subc->color>0)) {
        PetscCallMPI(MPI_Bcast(&sr->inertia1,1,MPIU_INT,0,child));
        PetscCallMPI(MPI_Bcast(&sr->int1,1,MPIU_REAL,0,child));
      } else sr_glob->inertia1 = sr->inertia1;
    }

    /* last process in eps comm computes inertia1 */
    if (ctx->npart==1 || ((sr->dir>0 && ctx->subc->color==ctx->npart-1) || (sr->dir<0 && ctx->subc->color==0))) {
      PetscCall(EPSSliceGetInertia(eps,sr->int1,&sr->inertia1,ctx->detect?&zeros:NULL));
      PetscCheck(zeros==0,((PetscObject)eps)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in an interval endpoint defined by user");
      if (!rank && sr->inertia0==-1) {
        sr->inertia0 = sr->inertia1; sr->int0 = sr->int1;
        PetscCall(PetscMPIIntCast(ctx->subc->color-sr->dir,&aux));
        PetscCallMPI(MPI_Isend(&sr->inertia0,1,MPIU_INT,aux,0,ctx->commrank,&req));
        PetscCallMPI(MPI_Isend(&sr->int0,1,MPIU_REAL,aux,0,ctx->commrank,&req));
      }
      if (sr->hasEnd) {
        sr->dir = -sr->dir; r = sr->int0; sr->int0 = sr->int1; sr->int1 = r;
        i = sr->inertia0; sr->inertia0 = sr->inertia1; sr->inertia1 = i;
      }
    }

    /* number of eigenvalues in interval */
    sr->numEigs = (sr->dir)*(sr->inertia1 - sr->inertia0);
    if (ctx->npart>1) {
      /* memory allocate for subinterval eigenpairs */
      PetscCall(EPSSliceAllocateSolution(eps,1));
    }
    dssz = eps->ncv+1;
    PetscCall(DSGetParallel(ctx->eps->ds,&ptype));
    PetscCall(DSSetParallel(eps->ds,ptype));
    PetscCall(DSGetMethod(ctx->eps->ds,&method));
    PetscCall(DSSetMethod(eps->ds,method));
  }
  PetscCall(DSSetType(eps->ds,DSHEP));
  PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
  PetscCall(DSAllocate(eps->ds,dssz));
  /* keep state of subcomm matrices to check that the user does not modify them */
  PetscCall(EPSGetOperators(eps,&A,&B));
  PetscCall(MatGetState(A,&ctx->Astate));
  PetscCall(PetscObjectGetId((PetscObject)A,&ctx->Aid));
  if (B) {
    PetscCall(MatGetState(B,&ctx->Bstate));
    PetscCall(PetscObjectGetId((PetscObject)B,&ctx->Bid));
  } else {
    ctx->Bstate=0;
    ctx->Bid=0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSliceGatherEigenVectors(EPS eps)
{
  Vec             v,vg,v_loc;
  IS              is1,is2;
  VecScatter      vec_sc;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        nloc,m0,n0,i,si,idx,*idx1,*idx2,j;
  PetscScalar     *array;
  EPS_SR          sr_loc;
  BV              V_loc;

  PetscFunctionBegin;
  sr_loc = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;
  V_loc = sr_loc->V;

  /* Gather parallel eigenvectors */
  PetscCall(BVGetColumn(eps->V,0,&v));
  PetscCall(VecGetOwnershipRange(v,&n0,&m0));
  PetscCall(BVRestoreColumn(eps->V,0,&v));
  PetscCall(BVGetColumn(ctx->eps->V,0,&v));
  PetscCall(VecGetLocalSize(v,&nloc));
  PetscCall(BVRestoreColumn(ctx->eps->V,0,&v));
  PetscCall(PetscMalloc2(m0-n0,&idx1,m0-n0,&idx2));
  PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)eps),nloc,PETSC_DECIDE,&vg));
  idx = -1;
  for (si=0;si<ctx->npart;si++) {
    j = 0;
    for (i=n0;i<m0;i++) {
      idx1[j]   = i;
      idx2[j++] = i+eps->n*si;
    }
    PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)eps),(m0-n0),idx1,PETSC_COPY_VALUES,&is1));
    PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)eps),(m0-n0),idx2,PETSC_COPY_VALUES,&is2));
    PetscCall(BVGetColumn(eps->V,0,&v));
    PetscCall(VecScatterCreate(v,is1,vg,is2,&vec_sc));
    PetscCall(BVRestoreColumn(eps->V,0,&v));
    PetscCall(ISDestroy(&is1));
    PetscCall(ISDestroy(&is2));
    for (i=0;i<ctx->nconv_loc[si];i++) {
      PetscCall(BVGetColumn(eps->V,++idx,&v));
      if (ctx->subc->color==si) {
        PetscCall(BVGetColumn(V_loc,i,&v_loc));
        PetscCall(VecGetArray(v_loc,&array));
        PetscCall(VecPlaceArray(vg,array));
      }
      PetscCall(VecScatterBegin(vec_sc,vg,v,INSERT_VALUES,SCATTER_REVERSE));
      PetscCall(VecScatterEnd(vec_sc,vg,v,INSERT_VALUES,SCATTER_REVERSE));
      if (ctx->subc->color==si) {
        PetscCall(VecResetArray(vg));
        PetscCall(VecRestoreArray(v_loc,&array));
        PetscCall(BVRestoreColumn(V_loc,i,&v_loc));
      }
      PetscCall(BVRestoreColumn(eps->V,idx,&v));
    }
    PetscCall(VecScatterDestroy(&vec_sc));
  }
  PetscCall(PetscFree2(idx1,idx2));
  PetscCall(VecDestroy(&vg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPSComputeVectors_Slice - Recover Eigenvectors from subcomunicators
 */
PetscErrorCode EPSComputeVectors_Slice(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (ctx->global && ctx->npart>1) {
    PetscCall(EPSComputeVectors(ctx->eps));
    PetscCall(EPSSliceGatherEigenVectors(eps));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSliceGetInertias(EPS eps,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i=0,j,tmpi;
  PetscReal       v,tmpr;
  EPS_shift       s;

  PetscFunctionBegin;
  PetscCheck(eps->state,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  PetscCheck(ctx->sr,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (!ctx->sr->s0) {  /* EPSSolve not called yet */
    *n = 2;
  } else {
    *n = 1;
    s = ctx->sr->s0;
    while (s) {
      (*n)++;
      s = s->neighb[1];
    }
  }
  PetscCall(PetscMalloc1(*n,shifts));
  PetscCall(PetscMalloc1(*n,inertias));
  if (!ctx->sr->s0) {  /* EPSSolve not called yet */
    (*shifts)[0]   = ctx->sr->int0;
    (*shifts)[1]   = ctx->sr->int1;
    (*inertias)[0] = ctx->sr->inertia0;
    (*inertias)[1] = ctx->sr->inertia1;
  } else {
    s = ctx->sr->s0;
    while (s) {
      (*shifts)[i]     = s->value;
      (*inertias)[i++] = s->inertia;
      s = s->neighb[1];
    }
    (*shifts)[i]   = ctx->sr->int1;
    (*inertias)[i] = ctx->sr->inertia1;
  }
  /* remove possible duplicate in last position */
  if ((*shifts)[(*n)-1]==(*shifts)[(*n)-2]) (*n)--;
  /* sort result */
  for (i=0;i<*n;i++) {
    v = (*shifts)[i];
    for (j=i+1;j<*n;j++) {
      if (v > (*shifts)[j]) {
        SlepcSwap((*shifts)[i],(*shifts)[j],tmpr);
        SlepcSwap((*inertias)[i],(*inertias)[j],tmpi);
        v = (*shifts)[i];
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSliceGatherSolution(EPS eps)
{
  PetscMPIInt     rank,nproc,*disp,*ns_loc,aux;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,idx,j,*perm_loc,off=0,*inertias_loc,ns;
  PetscScalar     *eigr_loc;
  EPS_SR          sr_loc;
  PetscReal       *shifts_loc;
  MPI_Comm        child;

  PetscFunctionBegin;
  eps->nconv = 0;
  for (i=0;i<ctx->npart;i++) eps->nconv += ctx->nconv_loc[i];
  sr_loc = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;

  /* Gather the shifts used and the inertias computed */
  PetscCall(EPSSliceGetInertias(ctx->eps,&ns,&shifts_loc,&inertias_loc));
  if (ctx->sr->dir>0 && shifts_loc[ns-1]==sr_loc->int1 && ctx->subc->color<ctx->npart-1) ns--;
  if (ctx->sr->dir<0 && shifts_loc[ns-1]==sr_loc->int0 && ctx->subc->color>0) {
    ns--;
    for (i=0;i<ns;i++) {
      inertias_loc[i] = inertias_loc[i+1];
      shifts_loc[i] = shifts_loc[i+1];
    }
  }
  PetscCall(PetscMalloc1(ctx->npart,&ns_loc));
  PetscCall(PetscSubcommGetChild(ctx->subc,&child));
  PetscCallMPI(MPI_Comm_rank(child,&rank));
  PetscCall(PetscMPIIntCast(ns,&aux));
  if (!rank) PetscCallMPI(MPI_Allgather(&aux,1,MPI_INT,ns_loc,1,MPI_INT,ctx->commrank));
  PetscCall(PetscMPIIntCast(ctx->npart,&aux));
  PetscCallMPI(MPI_Bcast(ns_loc,aux,MPI_INT,0,child));
  ctx->nshifts = 0;
  for (i=0;i<ctx->npart;i++) ctx->nshifts += ns_loc[i];
  PetscCall(PetscFree(ctx->inertias));
  PetscCall(PetscFree(ctx->shifts));
  PetscCall(PetscMalloc1(ctx->nshifts,&ctx->inertias));
  PetscCall(PetscMalloc1(ctx->nshifts,&ctx->shifts));

  /* Gather eigenvalues (same ranks have fully set of eigenvalues)*/
  eigr_loc = sr_loc->eigr;
  perm_loc = sr_loc->perm;
  PetscCallMPI(MPI_Comm_size(((PetscObject)eps)->comm,&nproc));
  PetscCall(PetscMalloc1(ctx->npart,&disp));
  disp[0] = 0;
  for (i=1;i<ctx->npart;i++) disp[i] = disp[i-1]+ctx->nconv_loc[i-1];
  if (nproc%ctx->npart==0) { /* subcommunicators with the same size */
    PetscCall(PetscMPIIntCast(sr_loc->numEigs,&aux));
    PetscCallMPI(MPI_Allgatherv(eigr_loc,aux,MPIU_SCALAR,eps->eigr,ctx->nconv_loc,disp,MPIU_SCALAR,ctx->commrank)); /* eigenvalues */
    PetscCallMPI(MPI_Allgatherv(perm_loc,aux,MPIU_INT,eps->perm,ctx->nconv_loc,disp,MPIU_INT,ctx->commrank)); /* perm */
    for (i=1;i<ctx->npart;i++) disp[i] = disp[i-1]+ns_loc[i-1];
    PetscCall(PetscMPIIntCast(ns,&aux));
    PetscCallMPI(MPI_Allgatherv(shifts_loc,aux,MPIU_REAL,ctx->shifts,ns_loc,disp,MPIU_REAL,ctx->commrank)); /* shifts */
    PetscCallMPI(MPI_Allgatherv(inertias_loc,aux,MPIU_INT,ctx->inertias,ns_loc,disp,MPIU_INT,ctx->commrank)); /* inertias */
    PetscCallMPI(MPIU_Allreduce(&sr_loc->itsKs,&eps->its,1,MPIU_INT,MPI_SUM,ctx->commrank));
  } else { /* subcommunicators with different size */
    if (!rank) {
      PetscCall(PetscMPIIntCast(sr_loc->numEigs,&aux));
      PetscCallMPI(MPI_Allgatherv(eigr_loc,aux,MPIU_SCALAR,eps->eigr,ctx->nconv_loc,disp,MPIU_SCALAR,ctx->commrank)); /* eigenvalues */
      PetscCallMPI(MPI_Allgatherv(perm_loc,aux,MPIU_INT,eps->perm,ctx->nconv_loc,disp,MPIU_INT,ctx->commrank)); /* perm */
      for (i=1;i<ctx->npart;i++) disp[i] = disp[i-1]+ns_loc[i-1];
      PetscCall(PetscMPIIntCast(ns,&aux));
      PetscCallMPI(MPI_Allgatherv(shifts_loc,aux,MPIU_REAL,ctx->shifts,ns_loc,disp,MPIU_REAL,ctx->commrank)); /* shifts */
      PetscCallMPI(MPI_Allgatherv(inertias_loc,aux,MPIU_INT,ctx->inertias,ns_loc,disp,MPIU_INT,ctx->commrank)); /* inertias */
      PetscCallMPI(MPIU_Allreduce(&sr_loc->itsKs,&eps->its,1,MPIU_INT,MPI_SUM,ctx->commrank));
    }
    PetscCall(PetscMPIIntCast(eps->nconv,&aux));
    PetscCallMPI(MPI_Bcast(eps->eigr,aux,MPIU_SCALAR,0,child));
    PetscCallMPI(MPI_Bcast(eps->perm,aux,MPIU_INT,0,child));
    PetscCall(PetscMPIIntCast(ctx->nshifts,&aux));
    PetscCallMPI(MPI_Bcast(ctx->shifts,aux,MPIU_REAL,0,child));
    PetscCallMPI(MPI_Bcast(ctx->inertias,aux,MPIU_INT,0,child));
    PetscCallMPI(MPI_Bcast(&eps->its,1,MPIU_INT,0,child));
  }
  /* Update global array eps->perm */
  idx = ctx->nconv_loc[0];
  for (i=1;i<ctx->npart;i++) {
    off += ctx->nconv_loc[i-1];
    for (j=0;j<ctx->nconv_loc[i];j++) eps->perm[idx++] += off;
  }

  /* Gather parallel eigenvectors */
  PetscCall(PetscFree(ns_loc));
  PetscCall(PetscFree(disp));
  PetscCall(PetscFree(shifts_loc));
  PetscCall(PetscFree(inertias_loc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Fills the fields of a shift structure
*/
static PetscErrorCode EPSCreateShift(EPS eps,PetscReal val,EPS_shift neighb0,EPS_shift neighb1)
{
  EPS_shift       s,*pending2;
  PetscInt        i;
  EPS_SR          sr;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  if ((neighb0 && val==neighb0->value) || (neighb1 && val==neighb1->value)) {
    sr->nPend++;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(PetscNew(&s));
  s->value = val;
  s->neighb[0] = neighb0;
  if (neighb0) neighb0->neighb[1] = s;
  s->neighb[1] = neighb1;
  if (neighb1) neighb1->neighb[0] = s;
  s->comp[0] = PETSC_FALSE;
  s->comp[1] = PETSC_FALSE;
  s->index = -1;
  s->neigs = 0;
  s->nconv[0] = s->nconv[1] = 0;
  s->nsch[0] = s->nsch[1]=0;
  /* Inserts in the stack of pending shifts */
  /* If needed, the array is resized */
  if (sr->nPend >= sr->maxPend) {
    sr->maxPend *= 2;
    PetscCall(PetscMalloc1(sr->maxPend,&pending2));
    for (i=0;i<sr->nPend;i++) pending2[i] = sr->pending[i];
    PetscCall(PetscFree(sr->pending));
    sr->pending = pending2;
  }
  sr->pending[sr->nPend++]=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Prepare for Rational Krylov update */
static PetscErrorCode EPSPrepareRational(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        dir,i,k,ld,nv;
  PetscScalar     *A;
  EPS_SR          sr = ctx->sr;
  Vec             v;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  dir = (sr->sPres->neighb[0] == sr->sPrev)?1:-1;
  dir*=sr->dir;
  k = 0;
  for (i=0;i<sr->nS;i++) {
    if (dir*PetscRealPart(sr->S[i])>0.0) {
      sr->S[k] = sr->S[i];
      sr->S[sr->nS+k] = sr->S[sr->nS+i];
      PetscCall(BVGetColumn(sr->Vnext,k,&v));
      PetscCall(BVCopyVec(eps->V,eps->nconv+i,v));
      PetscCall(BVRestoreColumn(sr->Vnext,k,&v));
      k++;
      if (k>=sr->nS/2) break;
    }
  }
  /* Copy to DS */
  PetscCall(DSSetCompact(eps->ds,PETSC_FALSE));  /* make sure DS_MAT_A is allocated */
  PetscCall(DSGetArray(eps->ds,DS_MAT_A,&A));
  PetscCall(PetscArrayzero(A,ld*ld));
  for (i=0;i<k;i++) {
    A[i*(1+ld)] = sr->S[i];
    A[k+i*ld] = sr->S[sr->nS+i];
  }
  sr->nS = k;
  PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&A));
  PetscCall(DSGetDimensions(eps->ds,&nv,NULL,NULL,NULL));
  PetscCall(DSSetDimensions(eps->ds,nv,0,k));
  /* Append u to V */
  PetscCall(BVGetColumn(sr->Vnext,sr->nS,&v));
  PetscCall(BVCopyVec(eps->V,sr->nv,v));
  PetscCall(BVRestoreColumn(sr->Vnext,sr->nS,&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Provides next shift to be computed */
static PetscErrorCode EPSExtractShift(EPS eps)
{
  PetscInt        iner,zeros=0;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr;
  PetscReal       newShift,diam,ptol;
  EPS_shift       sPres;

  PetscFunctionBegin;
  sr = ctx->sr;
  if (sr->nPend > 0) {
    if (sr->sPres==sr->pending[sr->nPend-1]) {
      eps->reason = EPS_CONVERGED_ITERATING;
      eps->its = 0;
      sr->nPend--;
      sr->sPres->rep = PETSC_TRUE;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
    sr->sPrev = sr->sPres;
    sr->sPres = sr->pending[--sr->nPend];
    sPres = sr->sPres;
    PetscCall(EPSSliceGetInertia(eps,sPres->value,&iner,ctx->detect?&zeros:NULL));
    if (zeros) {
      diam = PetscMin(PetscAbsReal(sPres->neighb[0]->value-sPres->value),PetscAbsReal(sPres->neighb[1]->value-sPres->value));
      ptol = PetscMin(SLICE_PTOL,diam/2);
      newShift = sPres->value*(1.0+ptol);
      if (sr->dir*(sPres->neighb[0] && newShift-sPres->neighb[0]->value) < 0) newShift = (sPres->value+sPres->neighb[0]->value)/2;
      else if (sPres->neighb[1] && sr->dir*(sPres->neighb[1]->value-newShift) < 0) newShift = (sPres->value+sPres->neighb[1]->value)/2;
      PetscCall(EPSSliceGetInertia(eps,newShift,&iner,&zeros));
      PetscCheck(zeros==0,((PetscObject)eps)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",(double)newShift);
      sPres->value = newShift;
    }
    sr->sPres->inertia = iner;
    eps->target = sr->sPres->value;
    eps->reason = EPS_CONVERGED_ITERATING;
    eps->its = 0;
  } else sr->sPres = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Symmetric KrylovSchur adapted to spectrum slicing:
   Allows searching an specific amount of eigenvalues in the subintervals left and right.
   Returns whether the search has succeeded
*/
static PetscErrorCode EPSKrylovSchur_Slice(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,k,l,ld,nv,*iwork,j,count0,count1,iterCompl=0,n0,n1;
  Mat             U,Op,T;
  PetscScalar     *Q,*A;
  PetscReal       *a,*b,beta,lambda;
  EPS_shift       sPres;
  PetscBool       breakdown,complIterating,sch0,sch1;
  EPS_SR          sr = ctx->sr;

  PetscFunctionBegin;
  /* Spectrum slicing data */
  sPres = sr->sPres;
  complIterating =PETSC_FALSE;
  sch1 = sch0 = PETSC_TRUE;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(PetscMalloc1(2*ld,&iwork));
  count0=0;count1=0; /* Found on both sides */
  if (!sPres->rep && sr->nS > 0 && (sPres->neighb[0] == sr->sPrev || sPres->neighb[1] == sr->sPrev)) {
    /* Rational Krylov */
    PetscCall(DSTranslateRKS(eps->ds,sr->sPrev->value-sPres->value));
    PetscCall(DSGetDimensions(eps->ds,NULL,NULL,&l,NULL));
    PetscCall(DSSetDimensions(eps->ds,l+1,0,0));
    PetscCall(BVSetActiveColumns(eps->V,0,l+1));
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(BVMultInPlace(eps->V,U,0,l+1));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));
  } else {
    /* Get the starting Lanczos vector */
    PetscCall(EPSGetStartVector(eps,0,NULL));
    l = 0;
  }
  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++; sr->itsKs++;
    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    PetscCall(DSGetMat(eps->ds,DS_MAT_T,&T));
    PetscCall(STGetOperator(eps->st,&Op));
    PetscCall(BVMatLanczos(eps->V,Op,T,eps->nconv+l,&nv,&beta,&breakdown));
    PetscCall(STRestoreOperator(eps->st,&Op));
    sr->nv = nv;
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_T,&T));
    PetscCall(DSSetDimensions(eps->ds,nv,eps->nconv,eps->nconv+l));
    if (l==0) PetscCall(DSSetState(eps->ds,DS_STATE_INTERMEDIATE));
    else PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
    PetscCall(BVSetActiveColumns(eps->V,eps->nconv,nv));

    /* Solve projected problem and compute residual norm estimates */
    if (eps->its == 1 && l > 0) { /* After rational update, DS_MAT_A is available */
      PetscCall(DSGetArray(eps->ds,DS_MAT_A,&A));
      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
      b = a + ld;
      k = eps->nconv+l;
      A[k*ld+k-1] = A[(k-1)*ld+k];
      A[k*ld+k] = a[k];
      for (j=k+1; j< nv; j++) {
        A[j*ld+j] = a[j];
        A[j*ld+j-1] = b[j-1] ;
        A[(j-1)*ld+j] = b[j-1];
      }
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_A,&A));
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
      PetscCall(DSSolve(eps->ds,eps->eigr,NULL));
      PetscCall(DSSort(eps->ds,eps->eigr,NULL,NULL,NULL,NULL));
      PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
    } else { /* Restart */
      PetscCall(DSSolve(eps->ds,eps->eigr,NULL));
      PetscCall(DSSort(eps->ds,eps->eigr,NULL,NULL,NULL,NULL));
    }
    PetscCall(DSSynchronize(eps->ds,eps->eigr,NULL));

    /* Residual */
    PetscCall(EPSKrylovConvergence(eps,PETSC_TRUE,eps->nconv,nv-eps->nconv,beta,0.0,1.0,&k));
    /* Checking values obtained for completing */
    for (i=0;i<k;i++) {
      sr->back[i]=eps->eigr[i];
    }
    PetscCall(STBackTransform(eps->st,k,sr->back,eps->eigi));
    count0=count1=0;
    for (i=0;i<k;i++) {
      lambda = PetscRealPart(sr->back[i]);
      if ((sr->dir*(sPres->value - lambda) > 0) && (sr->dir*(lambda - sPres->ext[0]) > 0)) count0++;
      if ((sr->dir*(lambda - sPres->value) > 0) && (sr->dir*(sPres->ext[1] - lambda) > 0)) count1++;
    }
    if (k>eps->nev && eps->ncv-k<5) eps->reason = EPS_CONVERGED_TOL;
    else {
      /* Checks completion */
      if ((!sch0||count0 >= sPres->nsch[0]) && (!sch1 ||count1 >= sPres->nsch[1])) {
        eps->reason = EPS_CONVERGED_TOL;
      } else {
        if (!complIterating && eps->its >= eps->max_it) eps->reason = EPS_DIVERGED_ITS;
        if (complIterating) {
          if (--iterCompl <= 0) eps->reason = EPS_DIVERGED_ITS;
        } else if (k >= eps->nev) {
          n0 = sPres->nsch[0]-count0;
          n1 = sPres->nsch[1]-count1;
          if (sr->iterCompl>0 && ((n0>0 && n0<= sr->nMAXCompl)||(n1>0&&n1<=sr->nMAXCompl))) {
            /* Iterating for completion*/
            complIterating = PETSC_TRUE;
            if (n0 >sr->nMAXCompl)sch0 = PETSC_FALSE;
            if (n1 >sr->nMAXCompl)sch1 = PETSC_FALSE;
            iterCompl = sr->iterCompl;
          } else eps->reason = EPS_CONVERGED_TOL;
        }
      }
    }
    /* Update l */
    if (eps->reason == EPS_CONVERGED_ITERATING) l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
    else l = nv-k;
    if (breakdown) l=0;
    if (!ctx->lock && l>0 && eps->reason == EPS_CONVERGED_ITERATING) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      if (breakdown) {
        /* Start a new Lanczos factorization */
        PetscCall(PetscInfo(eps,"Breakdown in Krylov-Schur method (it=%" PetscInt_FMT " norm=%g)\n",eps->its,(double)beta));
        PetscCall(EPSGetStartVector(eps,k,&breakdown));
        if (breakdown) {
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
        }
      } else {
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&a));
        PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Q));
        b = a + ld;
        for (i=k;i<k+l;i++) {
          a[i] = PetscRealPart(eps->eigr[i]);
          b[i] = PetscRealPart(Q[nv-1+i*ld]*beta);
        }
        PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&a));
        PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Q));
      }
    }
    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(BVMultInPlace(eps->V,U,eps->nconv,k+l));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));

    /* Normalize u and append it to V */
    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) PetscCall(BVCopyColumn(eps->V,nv,k+l));
    eps->nconv = k;
    if (eps->reason != EPS_CONVERGED_ITERATING) {
      /* Store approximated values for next shift */
      PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Q));
      sr->nS = l;
      for (i=0;i<l;i++) {
        sr->S[i] = eps->eigr[i+k];/* Diagonal elements */
        sr->S[i+l] = Q[nv-1+(i+k)*ld]*beta; /* Out of diagonal elements */
      }
      PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Q));
    }
  }
  /* Check for completion */
  for (i=0;i< eps->nconv; i++) {
    if (sr->dir*PetscRealPart(eps->eigr[i])>0) sPres->nconv[1]++;
    else sPres->nconv[0]++;
  }
  sPres->comp[0] = PetscNot(count0 < sPres->nsch[0]);
  sPres->comp[1] = PetscNot(count1 < sPres->nsch[1]);
  PetscCall(PetscInfo(eps,"Lanczos: %" PetscInt_FMT " evals in [%g,%g]%s and %" PetscInt_FMT " evals in [%g,%g]%s\n",count0,(double)(sr->dir==1?sPres->ext[0]:sPres->value),(double)(sr->dir==1?sPres->value:sPres->ext[0]),sPres->comp[0]?"*":"",count1,(double)(sr->dir==1?sPres->value:sPres->ext[1]),(double)(sr->dir==1?sPres->ext[1]:sPres->value),sPres->comp[1]?"*":""));
  PetscCheck(count0<=sPres->nsch[0] && count1<=sPres->nsch[1],PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia%s",ctx->detect?"":", consider using EPSKrylovSchurSetDetectZeros()");
  PetscCall(PetscFree(iwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Obtains value of subsequent shift
*/
static PetscErrorCode EPSGetNewShiftValue(EPS eps,PetscInt side,PetscReal *newS)
{
  PetscReal       lambda,d_prev;
  PetscInt        i,idxP;
  EPS_SR          sr;
  EPS_shift       sPres,s;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;
  if (sPres->neighb[side]) {
    /* Completing a previous interval */
    *newS = (sPres->value + sPres->neighb[side]->value)/2;
    if (PetscAbsReal(sPres->value - *newS)/PetscAbsReal(sPres->value)<=100*PETSC_SQRT_MACHINE_EPSILON) *newS = sPres->value;
  } else { /* (Only for side=1). Creating a new interval. */
    if (sPres->neigs==0) {/* No value has been accepted*/
      if (sPres->neighb[0]) {
        /* Multiplying by 10 the previous distance */
        *newS = sPres->value + 10*sr->dir*PetscAbsReal(sPres->value - sPres->neighb[0]->value);
        sr->nleap++;
        /* Stops when the interval is open and no values are found in the last 5 shifts (there might be infinite eigenvalues) */
        PetscCheck(sr->hasEnd || sr->nleap<=5,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Unable to compute the wanted eigenvalues with open interval");
      } else { /* First shift */
        PetscCheck(eps->nconv!=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"First shift renders no information");
        /* Unaccepted values give information for next shift */
        idxP=0;/* Number of values left from shift */
        for (i=0;i<eps->nconv;i++) {
          lambda = PetscRealPart(eps->eigr[i]);
          if (sr->dir*(lambda - sPres->value) <0) idxP++;
          else break;
        }
        /* Avoiding subtraction of eigenvalues (might be the same).*/
        if (idxP>0) {
          d_prev = PetscAbsReal(sPres->value - PetscRealPart(eps->eigr[0]))/(idxP+0.3);
        } else {
          d_prev = PetscAbsReal(sPres->value - PetscRealPart(eps->eigr[eps->nconv-1]))/(eps->nconv+0.3);
        }
        *newS = sPres->value + (sr->dir*d_prev*eps->nev)/2;
      }
    } else { /* Accepted values found */
      sr->nleap = 0;
      /* Average distance of values in previous subinterval */
      s = sPres->neighb[0];
      while (s && PetscAbs(s->inertia - sPres->inertia)==0) {
        s = s->neighb[0];/* Looking for previous shifts with eigenvalues within */
      }
      if (s) {
        d_prev = PetscAbsReal((sPres->value - s->value)/(sPres->inertia - s->inertia));
      } else { /* First shift. Average distance obtained with values in this shift */
        /* first shift might be too far from first wanted eigenvalue (no values found outside the interval)*/
        if (sr->dir*(PetscRealPart(sr->eigr[0])-sPres->value)>0 && PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0]))/PetscRealPart(sr->eigr[0])) > PetscSqrtReal(eps->tol)) {
          d_prev =  PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0])))/(sPres->neigs+0.3);
        } else {
          d_prev = PetscAbsReal(PetscRealPart(sr->eigr[sr->indexEig-1]) - sPres->value)/(sPres->neigs+0.3);
        }
      }
      /* Average distance is used for next shift by adding it to value on the right or to shift */
      if (sr->dir*(PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1]) - sPres->value)>0) {
        *newS = PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1])+ (sr->dir*d_prev*eps->nev)/2;
      } else { /* Last accepted value is on the left of shift. Adding to shift */
        *newS = sPres->value + (sr->dir*d_prev*eps->nev)/2;
      }
    }
    /* End of interval can not be surpassed */
    if (sr->dir*(sr->int1 - *newS) < 0) *newS = sr->int1;
  }/* of neighb[side]==null */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Function for sorting an array of real values
*/
static PetscErrorCode sortRealEigenvalues(PetscScalar *r,PetscInt *perm,PetscInt nr,PetscBool prev,PetscInt dir)
{
  PetscReal re;
  PetscInt  i,j,tmp;

  PetscFunctionBegin;
  if (!prev) for (i=0;i<nr;i++) perm[i] = i;
  /* Insertion sort */
  for (i=1;i<nr;i++) {
    re = PetscRealPart(r[perm[i]]);
    j = i-1;
    while (j>=0 && dir*(re - PetscRealPart(r[perm[j]])) <= 0) {
      tmp = perm[j]; perm[j] = perm[j+1]; perm[j+1] = tmp; j--;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Stores the pairs obtained since the last shift in the global arrays */
static PetscErrorCode EPSStoreEigenpairs(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  PetscReal       lambda,err,norm;
  PetscInt        i,count;
  PetscBool       iscayley;
  EPS_SR          sr = ctx->sr;
  EPS_shift       sPres;
  Vec             v,w;

  PetscFunctionBegin;
  sPres = sr->sPres;
  sPres->index = sr->indexEig;
  count = sr->indexEig;
  /* Back-transform */
  PetscCall(STBackTransform(eps->st,eps->nconv,eps->eigr,eps->eigi));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&iscayley));
  /* Sort eigenvalues */
  PetscCall(sortRealEigenvalues(eps->eigr,eps->perm,eps->nconv,PETSC_FALSE,sr->dir));
  /* Values stored in global array */
  for (i=0;i<eps->nconv;i++) {
    lambda = PetscRealPart(eps->eigr[eps->perm[i]]);
    err = eps->errest[eps->perm[i]];

    if (sr->dir*(lambda - sPres->ext[0]) > 0 && (sr->dir)*(sPres->ext[1] - lambda) > 0) {/* Valid value */
      PetscCheck(count<sr->numEigs,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Unexpected error in Spectrum Slicing");
      sr->eigr[count] = lambda;
      sr->errest[count] = err;
      /* Explicit purification */
      PetscCall(BVGetColumn(eps->V,eps->perm[i],&w));
      if (eps->purify) {
        PetscCall(BVGetColumn(sr->V,count,&v));
        PetscCall(STApply(eps->st,w,v));
        PetscCall(BVRestoreColumn(sr->V,count,&v));
      } else PetscCall(BVInsertVec(sr->V,count,w));
      PetscCall(BVRestoreColumn(eps->V,eps->perm[i],&w));
      PetscCall(BVNormColumn(sr->V,count,NORM_2,&norm));
      PetscCall(BVScaleColumn(sr->V,count,1.0/norm));
      count++;
    }
  }
  sPres->neigs = count - sr->indexEig;
  sr->indexEig = count;
  /* Global ordering array updating */
  PetscCall(sortRealEigenvalues(sr->eigr,sr->perm,count,PETSC_TRUE,sr->dir));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLookForDeflation(EPS eps)
{
  PetscReal       val;
  PetscInt        i,count0=0,count1=0;
  EPS_shift       sPres;
  PetscInt        ini,fin,k,idx0,idx1;
  EPS_SR          sr;
  Vec             v;
  EPS_KRYLOVSCHUR *ctx=(EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;

  if (sPres->neighb[0]) ini = (sr->dir)*(sPres->neighb[0]->inertia - sr->inertia0);
  else ini = 0;
  fin = sr->indexEig;
  /* Selection of ends for searching new values */
  if (!sPres->neighb[0]) sPres->ext[0] = sr->int0;/* First shift */
  else sPres->ext[0] = sPres->neighb[0]->value;
  if (!sPres->neighb[1]) {
    if (sr->hasEnd) sPres->ext[1] = sr->int1;
    else sPres->ext[1] = (sr->dir > 0)?PETSC_MAX_REAL:PETSC_MIN_REAL;
  } else sPres->ext[1] = sPres->neighb[1]->value;
  /* Selection of values between right and left ends */
  for (i=ini;i<fin;i++) {
    val=PetscRealPart(sr->eigr[sr->perm[i]]);
    /* Values to the right of left shift */
    if (sr->dir*(val - sPres->ext[1]) < 0) {
      if (sr->dir*(val - sPres->value) < 0) count0++;
      else count1++;
    } else break;
  }
  /* The number of values on each side are found */
  if (sPres->neighb[0]) {
    sPres->nsch[0] = (sr->dir)*(sPres->inertia - sPres->neighb[0]->inertia)-count0;
    PetscCheck(sPres->nsch[0]>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia%s",ctx->detect?"":", consider using EPSKrylovSchurSetDetectZeros()");
  } else sPres->nsch[0] = 0;

  if (sPres->neighb[1]) {
    sPres->nsch[1] = (sr->dir)*(sPres->neighb[1]->inertia - sPres->inertia) - count1;
    PetscCheck(sPres->nsch[1]>=0,PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia%s",ctx->detect?"":", consider using EPSKrylovSchurSetDetectZeros()");
  } else sPres->nsch[1] = (sr->dir)*(sr->inertia1 - sPres->inertia);

  /* Completing vector of indexes for deflation */
  idx0 = ini;
  idx1 = ini+count0+count1;
  k=0;
  for (i=idx0;i<idx1;i++) sr->idxDef[k++]=sr->perm[i];
  PetscCall(BVDuplicateResize(eps->V,k+eps->ncv+1,&sr->Vnext));
  PetscCall(BVSetNumConstraints(sr->Vnext,k));
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(sr->Vnext,-i-1,&v));
    PetscCall(BVCopyVec(sr->V,sr->idxDef[i],v));
    PetscCall(BVRestoreColumn(sr->Vnext,-i-1,&v));
  }

  /* For rational Krylov */
  if (!sr->sPres->rep && sr->nS>0 && (sr->sPrev == sr->sPres->neighb[0] || sr->sPrev == sr->sPres->neighb[1])) PetscCall(EPSPrepareRational(eps));
  eps->nconv = 0;
  /* Get rid of temporary Vnext */
  PetscCall(BVDestroy(&eps->V));
  eps->V = sr->Vnext;
  sr->Vnext = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EPSSolve_KrylovSchur_Slice(EPS eps)
{
  PetscInt         i,lds,ti;
  PetscReal        newS;
  EPS_KRYLOVSCHUR  *ctx=(EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR           sr=ctx->sr;
  Mat              A,B=NULL;
  PetscObjectState Astate,Bstate=0;
  PetscObjectId    Aid,Bid=0;

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));
  if (ctx->global) {
    PetscCall(EPSSolve_KrylovSchur_Slice(ctx->eps));
    ctx->eps->state = EPS_STATE_SOLVED;
    eps->reason = EPS_CONVERGED_TOL;
    if (ctx->npart>1) {
      /* Gather solution from subsolvers */
      PetscCall(EPSSliceGatherSolution(eps));
    } else {
      eps->nconv = sr->numEigs;
      eps->its   = ctx->eps->its;
      PetscCall(PetscFree(ctx->inertias));
      PetscCall(PetscFree(ctx->shifts));
      PetscCall(EPSSliceGetInertias(ctx->eps,&ctx->nshifts,&ctx->shifts,&ctx->inertias));
    }
  } else {
    if (ctx->npart==1) {
      sr->eigr   = ctx->eps->eigr;
      sr->eigi   = ctx->eps->eigi;
      sr->perm   = ctx->eps->perm;
      sr->errest = ctx->eps->errest;
      sr->V      = ctx->eps->V;
    }
    /* Check that the user did not modify subcomm matrices */
    PetscCall(EPSGetOperators(eps,&A,&B));
    PetscCall(MatGetState(A,&Astate));
    PetscCall(PetscObjectGetId((PetscObject)A,&Aid));
    if (B) {
      PetscCall(MatGetState(B,&Bstate));
      PetscCall(PetscObjectGetId((PetscObject)B,&Bid));
    }
    PetscCheck(Astate==ctx->Astate && (!B || Bstate==ctx->Bstate) && Aid==ctx->Aid && (!B || Bid==ctx->Bid),PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Subcomm matrices have been modified by user");
    /* Only with eigenvalues present in the interval ...*/
    if (sr->numEigs==0) {
      eps->reason = EPS_CONVERGED_TOL;
      PetscFunctionReturn(PETSC_SUCCESS);
    }
    /* Array of pending shifts */
    sr->maxPend = 100; /* Initial size */
    sr->nPend = 0;
    PetscCall(PetscMalloc1(sr->maxPend,&sr->pending));
    PetscCall(EPSCreateShift(eps,sr->int0,NULL,NULL));
    /* extract first shift */
    sr->sPrev = NULL;
    sr->sPres = sr->pending[--sr->nPend];
    sr->sPres->inertia = sr->inertia0;
    eps->target = sr->sPres->value;
    sr->s0 = sr->sPres;
    sr->indexEig = 0;
    /* Memory reservation for auxiliary variables */
    lds = PetscMin(eps->mpd,eps->ncv);
    PetscCall(PetscCalloc1(lds*lds,&sr->S));
    PetscCall(PetscMalloc1(eps->ncv,&sr->back));
    for (i=0;i<sr->numEigs;i++) {
      sr->eigr[i]   = 0.0;
      sr->eigi[i]   = 0.0;
      sr->errest[i] = 0.0;
      sr->perm[i]   = i;
    }
    /* Vectors for deflation */
    PetscCall(PetscMalloc1(sr->numEigs,&sr->idxDef));
    sr->indexEig = 0;
    /* Main loop */
    while (sr->sPres) {
      /* Search for deflation */
      PetscCall(EPSLookForDeflation(eps));
      /* KrylovSchur */
      PetscCall(EPSKrylovSchur_Slice(eps));

      PetscCall(EPSStoreEigenpairs(eps));
      /* Select new shift */
      if (!sr->sPres->comp[1]) {
        PetscCall(EPSGetNewShiftValue(eps,1,&newS));
        PetscCall(EPSCreateShift(eps,newS,sr->sPres,sr->sPres->neighb[1]));
      }
      if (!sr->sPres->comp[0]) {
        /* Completing earlier interval */
        PetscCall(EPSGetNewShiftValue(eps,0,&newS));
        PetscCall(EPSCreateShift(eps,newS,sr->sPres->neighb[0],sr->sPres));
      }
      /* Preparing for a new search of values */
      PetscCall(EPSExtractShift(eps));
    }

    /* Updating eps values prior to exit */
    PetscCall(PetscFree(sr->S));
    PetscCall(PetscFree(sr->idxDef));
    PetscCall(PetscFree(sr->pending));
    PetscCall(PetscFree(sr->back));
    PetscCall(BVDuplicateResize(eps->V,eps->ncv+1,&sr->Vnext));
    PetscCall(BVSetNumConstraints(sr->Vnext,0));
    PetscCall(BVDestroy(&eps->V));
    eps->V      = sr->Vnext;
    eps->nconv  = sr->indexEig;
    eps->reason = EPS_CONVERGED_TOL;
    eps->its    = sr->itsKs;
    eps->nds    = 0;
    if (sr->dir<0) {
      for (i=0;i<eps->nconv/2;i++) {
        ti = sr->perm[i]; sr->perm[i] = sr->perm[eps->nconv-1-i]; sr->perm[eps->nconv-1-i] = ti;
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
