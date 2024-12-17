/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to eigensolvers in ChASE.
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include <petsc/private/petscscalapack.h>

#if !defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
/* s */
#define CHASEchase_init_blockcyclic pschase_init_blockcyclic_
#define CHASEchase                  pschase_
#define CHASEchase_finalize         pschase_finalize_
#elif defined(PETSC_USE_REAL_DOUBLE)
/* d */
#define CHASEchase_init_blockcyclic pdchase_init_blockcyclic_
#define CHASEchase                  pdchase_
#define CHASEchase_finalize         pdchase_finalize_
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
/* c */
#define CHASEchase_init_blockcyclic pcchase_init_blockcyclic_
#define CHASEchase                  pcchase_
#define CHASEchase_finalize         pcchase_finalize_
#elif defined(PETSC_USE_REAL_DOUBLE)
/* z */
#define CHASEchase_init_blockcyclic pzchase_init_blockcyclic_
#define CHASEchase                  pzchase_
#define CHASEchase_finalize         pzchase_finalize_
#endif
#endif

SLEPC_EXTERN void CHASEchase_init_blockcyclic(PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscInt*,PetscScalar*,PetscReal*,PetscInt*,PetscInt*,char*,PetscInt*,PetscInt*,MPI_Comm*,PetscInt*);
SLEPC_EXTERN void CHASEchase(PetscInt*,PetscReal*,char*,char*,char*);
SLEPC_EXTERN void CHASEchase_finalize(PetscInt*);

typedef struct {
  Mat       As;    /* converted matrix */
  PetscInt  nex;   /* extra searching space size */
  PetscInt  deg;   /* initial degree of Chebyshev polynomial filter */
  PetscBool opt;   /* internal optimization of polynomial degree */
} EPS_ChASE;

static PetscErrorCode EPSSetUp_ChASE(EPS eps)
{
  EPS_ChASE      *ctx = (EPS_ChASE*)eps->data;
  Mat            A;
  PetscBool      isshift;
  PetscScalar    shift;

  PetscFunctionBegin;
  EPSCheckStandard(eps);
  EPSCheckHermitian(eps);
  EPSCheckNotStructured(eps);
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  PetscCheck(isshift,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support spectral transformations");
  if (eps->nev==0) eps->nev = 1;
  if (eps->ncv==PETSC_DETERMINE) eps->ncv = PetscMin(eps->n,PetscMax(2*eps->nev,eps->nev+15));
  else PetscCheck(eps->ncv>=eps->nev+1 || (eps->ncv==eps->nev && eps->ncv==eps->n),PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev+1");
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  PetscCheck(eps->which==EPS_SMALLEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver only supports computation of leftmost eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  PetscCall(EPSAllocateSolution(eps,0));

  ctx->nex = eps->ncv-eps->nev;
  if (ctx->deg==PETSC_DEFAULT) ctx->deg = 10;

  /* convert matrices */
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatConvert(A,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->As));
  PetscCall(STGetShift(eps->st,&shift));
  if (shift != 0.0) PetscCall(MatShift(ctx->As,-shift));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_ChASE(EPS eps)
{
  EPS_ChASE      *ctx = (EPS_ChASE*)eps->data;
  Mat            As = ctx->As,Vs,V;
  Mat_ScaLAPACK  *a = (Mat_ScaLAPACK*)As->data,*v;
  MPI_Comm       comm;
  PetscReal      *w = eps->errest;  /* used to store real eigenvalues */
  PetscInt       i,init=0,flg=1;
  char           grid_major='R',mode='X',opt=ctx->opt?'S':'X',qr='C';

  PetscFunctionBegin;
  PetscCheck(a->grid->npcol==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Current implementation of the ChASE interface only supports process grids with one column; use -mat_scalapack_grid_height <p> where <p> is the number of processes");
  PetscCall(BVGetMat(eps->V,&V));
  PetscCall(MatConvert(V,MATSCALAPACK,MAT_INITIAL_MATRIX,&Vs));
  PetscCall(BVRestoreMat(eps->V,&V));
  v = (Mat_ScaLAPACK*)Vs->data;

  /* initialization */
  PetscCall(PetscObjectGetComm((PetscObject)As,&comm));
  CHASEchase_init_blockcyclic(&a->N,&eps->nev,&ctx->nex,&a->mb,&a->nb,a->loc,&a->lld,v->loc,w,&a->grid->nprow,&a->grid->npcol,&grid_major,&a->rsrc,&a->csrc,&comm,&init);
  PetscCheck(init==1,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Problem initializing ChASE");

  /* solve */
  CHASEchase(&ctx->deg,&eps->tol,&mode,&opt,&qr);
  CHASEchase_finalize(&flg);
  PetscCheck(flg==0,PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Problem solving with ChASE");

  for (i=0;i<eps->ncv;i++) {
    eps->eigr[i]   = eps->errest[i];
    eps->errest[i] = eps->tol;
  }

  PetscCall(BVGetMat(eps->V,&V));
  PetscCall(MatConvert(Vs,MATDENSE,MAT_REUSE_MATRIX,&V));
  PetscCall(BVRestoreMat(eps->V,&V));
  PetscCall(MatDestroy(&Vs));

  eps->nconv  = eps->nev;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCHASESetDegree_ChASE(EPS eps,PetscInt deg,PetscBool opt)
{
  EPS_ChASE *ctx = (EPS_ChASE*)eps->data;

  PetscFunctionBegin;
  if (deg==PETSC_DEFAULT || deg==PETSC_DECIDE) {
    ctx->deg   = PETSC_DEFAULT;
    eps->state = EPS_STATE_INITIAL;
  } else {
    PetscCheck(deg>0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"The degree must be >0");
    ctx->deg = deg;
  }
  ctx->opt = opt;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCHASESetDegree - Sets the degree of the Chebyshev polynomial filter in the ChASE solver.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  deg - initial degree of Chebyshev polynomial filter
-  opt - internal optimization of polynomial degree

   Options Database Keys:
+  -eps_chase_degree - Sets the initial degree
-  -eps_chase_degree_opt - Enables/disables the optimization

   Level: advanced

.seealso: EPSCHASEGetDegree()
@*/
PetscErrorCode EPSCHASESetDegree(EPS eps,PetscInt deg,PetscBool opt)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,deg,2);
  PetscValidLogicalCollectiveBool(eps,opt,3);
  PetscTryMethod(eps,"EPSCHASESetDegree_C",(EPS,PetscInt,PetscBool),(eps,deg,opt));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSCHASEGetDegree_ChASE(EPS eps,PetscInt *deg,PetscBool *opt)
{
  EPS_ChASE *ctx = (EPS_ChASE*)eps->data;

  PetscFunctionBegin;
  if (deg) *deg = ctx->deg;
  if (opt) *opt = ctx->opt;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSCHASEGetDegree - Gets the degree of the Chebyshev polynomial filter used in the ChASE solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  deg - initial degree of Chebyshev polynomial filter
-  opt - internal optimization of polynomial degree

   Level: advanced

.seealso: EPSCHASESetDegree()
@*/
PetscErrorCode EPSCHASEGetDegree(EPS eps,PetscInt *deg,PetscBool *opt)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSCHASEGetDegree_C",(EPS,PetscInt*,PetscBool*),(eps,deg,opt));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_ChASE(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_ChASE   *ctx = (EPS_ChASE*)eps->data;
  PetscBool   flg1,flg2;
  PetscInt    deg;
  PetscBool   opt;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS ChASE Options");

    deg = ctx->deg;
    PetscCall(PetscOptionsInt("-eps_chase_degree","Initial degree of Chebyshev polynomial filter","EPSCHASESetDegree",deg,&deg,&flg1));
    opt = ctx->opt;
    PetscCall(PetscOptionsBool("-eps_chase_degree_opt","Internal optimization of polynomial degree","EPSCHASESetDegree",opt,&opt,&flg2));
    if (flg1 || flg2) PetscCall(EPSCHASESetDegree(eps,deg,opt));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_ChASE(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCHASESetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCHASEGetDegree_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_ChASE(EPS eps)
{
  EPS_ChASE   *ctx = (EPS_ChASE*)eps->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->As));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_ChASE(EPS eps,PetscViewer viewer)
{
  EPS_ChASE   *ctx = (EPS_ChASE*)eps->data;
  PetscBool   isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  initial degree of Chebyshev polynomial filter: %" PetscInt_FMT "\n",ctx->deg));
    if (ctx->opt) PetscCall(PetscViewerASCIIPrintf(viewer,"  internal optimization of polynomial degree enabled\n"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_ChASE(EPS eps)
{
  EPS_ChASE   *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;
  ctx->deg  = PETSC_DEFAULT;
  ctx->opt  = PETSC_TRUE;

  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_ChASE;
  eps->ops->setup          = EPSSetUp_ChASE;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->setfromoptions = EPSSetFromOptions_ChASE;
  eps->ops->destroy        = EPSDestroy_ChASE;
  eps->ops->reset          = EPSReset_ChASE;
  eps->ops->view           = EPSView_ChASE;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_NoFactor;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCHASESetDegree_C",EPSCHASESetDegree_ChASE));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSCHASEGetDegree_C",EPSCHASEGetDegree_ChASE));
  PetscFunctionReturn(PETSC_SUCCESS);
}
