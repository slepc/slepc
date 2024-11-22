/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lobpcg"

   Method: Locally Optimal Block Preconditioned Conjugate Gradient

   Algorithm:

       LOBPCG with soft and hard locking. Follows the implementation
       in BLOPEX [2].

   References:

       [1] A. V. Knyazev, "Toward the optimal preconditioned eigensolver:
           locally optimal block preconditioned conjugate gradient method",
           SIAM J. Sci. Comput. 23(2):517-541, 2001.

       [2] A. V. Knyazev et al., "Block Locally Optimal Preconditioned
           Eigenvalue Xolvers (BLOPEX) in Hypre and PETSc", SIAM J. Sci.
           Comput. 29(5):2224-2239, 2007.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

typedef struct {
  PetscInt  bs;        /* block size */
  PetscBool lock;      /* soft locking active/inactive */
  PetscReal restart;   /* restart parameter */
  PetscInt  guard;     /* number of guard vectors */
} EPS_LOBPCG;

static PetscErrorCode EPSSetDimensions_LOBPCG(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;
  PetscInt   k;

  PetscFunctionBegin;
  if (*nev==0) *nev = 1;
  k = PetscMax(3*ctx->bs,((*nev-1)/ctx->bs+3)*ctx->bs);
  if (*ncv!=PETSC_DETERMINE) { /* ncv set */
    PetscCheck(*ncv>=k,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv is not sufficiently large");
  } else *ncv = k;
  if (*mpd==PETSC_DETERMINE) *mpd = 3*ctx->bs;
  else PetscCheck(*mpd==3*ctx->bs,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"This solver does not allow a value of mpd different from 3*blocksize");
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetUp_LOBPCG(EPS eps)
{
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  if (!ctx->bs) ctx->bs = PetscMin(16,eps->nev?eps->nev:1);
  PetscCheck(eps->n-eps->nds>=5*ctx->bs,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The problem size is too small relative to the block size");
  PetscCall(EPSSetDimensions_LOBPCG(eps,&eps->nev,&eps->ncv,&eps->mpd));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  PetscCheck(eps->which==EPS_SMALLEST_REAL || eps->which==EPS_LARGEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only smallest real or largest real eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_THRESHOLD);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE);

  if (!ctx->restart) ctx->restart = 0.9;

  /* number of guard vectors */
  if (ctx->bs==1) ctx->guard = 0;
  else ctx->guard = PetscMin((PetscInt)((1.0-ctx->restart)*ctx->bs+0.45),ctx->bs-1);

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(EPS_SetInnerProduct(eps));
  PetscCall(DSSetType(eps->ds,DSGHEP));
  PetscCall(DSAllocate(eps->ds,eps->mpd));
  PetscCall(EPSSetWorkVecs(eps,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_LOBPCG(EPS eps)
{
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;
  PetscInt       i,j,k,nv,ini,nmat,nc,nconv,locked,its,prev=0;
  PetscReal      norm;
  PetscScalar    *eigr,dot;
  PetscBool      breakdown,countc,flip=PETSC_FALSE,checkprecond=PETSC_FALSE;
  Mat            A,B,M,V=NULL,W=NULL;
  Vec            v,z,w=eps->work[0];
  BV             X,Y=NULL,Z,R,P,AX,BX;
  SlepcSC        sc;

  PetscFunctionBegin;
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  if (nmat>1) PetscCall(STGetMatrix(eps->st,1,&B));
  else B = NULL;

  if (eps->which==EPS_LARGEST_REAL) {  /* flip spectrum */
    flip = PETSC_TRUE;
    PetscCall(DSGetSlepcSC(eps->ds,&sc));
    sc->comparison = SlepcCompareSmallestReal;
  }

  /* undocumented option to check for a positive-definite preconditioner (turn-off by default) */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-eps_lobpcg_checkprecond",&checkprecond,NULL));

  /* 1. Allocate memory */
  PetscCall(PetscCalloc1(3*ctx->bs,&eigr));
  PetscCall(BVDuplicateResize(eps->V,3*ctx->bs,&Z));
  PetscCall(BVDuplicateResize(eps->V,ctx->bs,&X));
  PetscCall(BVDuplicateResize(eps->V,ctx->bs,&R));
  PetscCall(BVDuplicateResize(eps->V,ctx->bs,&P));
  PetscCall(BVDuplicateResize(eps->V,ctx->bs,&AX));
  if (B) PetscCall(BVDuplicateResize(eps->V,ctx->bs,&BX));
  nc = eps->nds;
  if (nc>0 || eps->nev>ctx->bs-ctx->guard) PetscCall(BVDuplicateResize(eps->V,nc+eps->nev,&Y));
  if (nc>0) {
    for (j=0;j<nc;j++) {
      PetscCall(BVGetColumn(eps->V,-nc+j,&v));
      PetscCall(BVInsertVec(Y,j,v));
      PetscCall(BVRestoreColumn(eps->V,-nc+j,&v));
    }
    PetscCall(BVSetActiveColumns(Y,0,nc));
  }

  /* 2. Apply the constraints to the initial vectors */
  /* 3. B-orthogonalize initial vectors */
  for (k=eps->nini;k<eps->ncv-ctx->bs;k++) { /* Generate more initial vectors if necessary */
    PetscCall(BVSetRandomColumn(eps->V,k));
    PetscCall(BVOrthonormalizeColumn(eps->V,k,PETSC_TRUE,NULL,NULL));
  }
  nv = ctx->bs;
  PetscCall(BVSetActiveColumns(eps->V,0,nv));
  PetscCall(BVSetActiveColumns(Z,0,nv));
  PetscCall(BVCopy(eps->V,Z));
  PetscCall(BVCopy(Z,X));

  /* 4. Compute initial Ritz vectors */
  PetscCall(BVMatMult(X,A,AX));
  PetscCall(DSSetDimensions(eps->ds,nv,0,0));
  PetscCall(DSGetMat(eps->ds,DS_MAT_A,&M));
  PetscCall(BVMatProject(AX,NULL,X,M));
  if (flip) PetscCall(MatScale(M,-1.0));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&M));
  PetscCall(DSSetIdentity(eps->ds,DS_MAT_B));
  PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
  PetscCall(DSSolve(eps->ds,eigr,NULL));
  PetscCall(DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL));
  PetscCall(DSSynchronize(eps->ds,eigr,NULL));
  for (j=0;j<nv;j++) eps->eigr[j] = flip? -eigr[j]: eigr[j];
  PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
  PetscCall(DSGetMat(eps->ds,DS_MAT_X,&M));
  PetscCall(BVMultInPlace(X,M,0,nv));
  PetscCall(BVMultInPlace(AX,M,0,nv));
  PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&M));

  /* 5. Initialize range of active iterates */
  locked = 0;  /* hard-locked vectors, the leading locked columns of V are eigenvectors */
  nconv  = 0;  /* number of converged eigenvalues in the current block */
  its    = 0;  /* iterations for the current block */

  /* 6. Main loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {

    if (ctx->lock) {
      PetscCall(BVSetActiveColumns(R,nconv,ctx->bs));
      PetscCall(BVSetActiveColumns(AX,nconv,ctx->bs));
      if (B) PetscCall(BVSetActiveColumns(BX,nconv,ctx->bs));
    }

    /* 7. Compute residuals */
    ini = (ctx->lock)? nconv: 0;
    PetscCall(BVCopy(AX,R));
    if (B) PetscCall(BVMatMult(X,B,BX));
    for (j=ini;j<ctx->bs;j++) {
      PetscCall(BVGetColumn(R,j,&v));
      PetscCall(BVGetColumn(B?BX:X,j,&z));
      PetscCall(VecAXPY(v,-eps->eigr[locked+j],z));
      PetscCall(BVRestoreColumn(R,j,&v));
      PetscCall(BVRestoreColumn(B?BX:X,j,&z));
    }

    /* 8. Compute residual norms and update index set of active iterates */
    k = ini;
    countc = PETSC_TRUE;
    for (j=ini;j<ctx->bs;j++) {
      i = locked+j;
      PetscCall(BVGetColumn(R,j,&v));
      PetscCall(VecNorm(v,NORM_2,&norm));
      PetscCall(BVRestoreColumn(R,j,&v));
      PetscCall((*eps->converged)(eps,eps->eigr[i],eps->eigi[i],norm,&eps->errest[i],eps->convergedctx));
      if (countc) {
        if (eps->errest[i] < eps->tol) k++;
        else countc = PETSC_FALSE;
      }
      if (!countc && !eps->trackall) break;
    }
    nconv = k;
    eps->nconv = locked + nconv;
    if (its) PetscCall(EPSMonitor(eps,eps->its+its,eps->nconv,eps->eigr,eps->eigi,eps->errest,locked+ctx->bs));
    PetscCall((*eps->stopping)(eps,eps->its+its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx));
    if (eps->reason != EPS_CONVERGED_ITERATING || nconv >= ctx->bs-ctx->guard) {
      PetscCall(BVSetActiveColumns(eps->V,locked,eps->nconv));
      PetscCall(BVSetActiveColumns(X,0,nconv));
      PetscCall(BVCopy(X,eps->V));
    }
    if (eps->reason != EPS_CONVERGED_ITERATING) {
      break;
    } else if (nconv >= ctx->bs-ctx->guard) {
      eps->its += its-1;
      its = 0;
    } else its++;

    if (nconv >= ctx->bs-ctx->guard) {  /* force hard locking of vectors and compute new R */

      /* extend constraints */
      PetscCall(BVSetActiveColumns(Y,nc+locked,nc+locked+nconv));
      PetscCall(BVCopy(X,Y));
      PetscCall(BVSetActiveColumns(Y,0,nc+locked+nconv));

      /* shift work BV's */
      for (j=nconv;j<ctx->bs;j++) {
        PetscCall(BVCopyColumn(X,j,j-nconv));
        PetscCall(BVCopyColumn(R,j,j-nconv));
        PetscCall(BVCopyColumn(P,j,j-nconv));
        PetscCall(BVCopyColumn(AX,j,j-nconv));
        if (B) PetscCall(BVCopyColumn(BX,j,j-nconv));
      }

      /* set new initial vectors */
      PetscCall(BVSetActiveColumns(eps->V,locked+ctx->bs,locked+ctx->bs+nconv));
      PetscCall(BVSetActiveColumns(X,ctx->bs-nconv,ctx->bs));
      PetscCall(BVCopy(eps->V,X));
      for (j=ctx->bs-nconv;j<ctx->bs;j++) {
        PetscCall(BVGetColumn(X,j,&v));
        PetscCall(BVOrthogonalizeVec(Y,v,NULL,&norm,&breakdown));
        if (norm>0.0 && !breakdown) PetscCall(VecScale(v,1.0/norm));
        else {
          PetscCall(PetscInfo(eps,"Orthogonalization of initial vector failed\n"));
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
        PetscCall(BVRestoreColumn(X,j,&v));
      }
      locked += nconv;
      nconv = 0;
      PetscCall(BVSetActiveColumns(X,nconv,ctx->bs));

      /* B-orthogonalize initial vectors */
      PetscCall(BVOrthogonalize(X,NULL));
      PetscCall(BVSetActiveColumns(Z,nconv,ctx->bs));
      PetscCall(BVSetActiveColumns(AX,nconv,ctx->bs));
      PetscCall(BVCopy(X,Z));

      /* compute initial Ritz vectors */
      nv = ctx->bs;
      PetscCall(BVMatMult(X,A,AX));
      PetscCall(DSSetDimensions(eps->ds,nv,0,0));
      PetscCall(DSGetMat(eps->ds,DS_MAT_A,&M));
      PetscCall(BVMatProject(AX,NULL,X,M));
      if (flip) PetscCall(MatScale(M,-1.0));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&M));
      PetscCall(DSSetIdentity(eps->ds,DS_MAT_B));
      PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
      PetscCall(DSSolve(eps->ds,eigr,NULL));
      PetscCall(DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL));
      PetscCall(DSSynchronize(eps->ds,eigr,NULL));
      for (j=0;j<nv;j++) if (locked+j<eps->ncv) eps->eigr[locked+j] = flip? -eigr[j]: eigr[j];
      PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));
      PetscCall(DSGetMat(eps->ds,DS_MAT_X,&M));
      PetscCall(BVMultInPlace(X,M,0,nv));
      PetscCall(BVMultInPlace(AX,M,0,nv));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&M));

      continue;   /* skip the rest of the iteration */
    }

    ini = (ctx->lock)? nconv: 0;
    if (ctx->lock) {
      PetscCall(BVSetActiveColumns(R,nconv,ctx->bs));
      PetscCall(BVSetActiveColumns(P,nconv,ctx->bs));
      PetscCall(BVSetActiveColumns(AX,nconv,ctx->bs));
      if (B) PetscCall(BVSetActiveColumns(BX,nconv,ctx->bs));
    }

    /* 9. Apply preconditioner to the residuals */
    PetscCall(BVGetMat(R,&V));
    if (prev != ctx->bs-ini) {
      prev = ctx->bs-ini;
      PetscCall(MatDestroy(&W));
      PetscCall(MatDuplicate(V,MAT_SHARE_NONZERO_PATTERN,&W));
    }
    PetscCall(STApplyMat(eps->st,V,W));
    if (checkprecond) {
      for (j=ini;j<ctx->bs;j++) {
        PetscCall(MatDenseGetColumnVecRead(V,j-ini,&v));
        PetscCall(MatDenseGetColumnVecRead(W,j-ini,&w));
        PetscCall(VecDot(v,w,&dot));
        PetscCall(MatDenseRestoreColumnVecRead(W,j-ini,&w));
        PetscCall(MatDenseRestoreColumnVecRead(V,j-ini,&v));
        if (PetscRealPart(dot)<0.0) {
          PetscCall(PetscInfo(eps,"The preconditioner is not positive-definite\n"));
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
      }
    }
    if (nc+locked>0) {
      for (j=ini;j<ctx->bs;j++) {
        PetscCall(MatDenseGetColumnVecWrite(W,j-ini,&w));
        PetscCall(BVOrthogonalizeVec(Y,w,NULL,&norm,&breakdown));
        if (norm>0.0 && !breakdown) PetscCall(VecScale(w,1.0/norm));
        PetscCall(MatDenseRestoreColumnVecWrite(W,j-ini,&w));
        if (norm<=0.0 || breakdown) {
          PetscCall(PetscInfo(eps,"Orthogonalization of preconditioned residual failed\n"));
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
      }
    }
    PetscCall(MatCopy(W,V,SAME_NONZERO_PATTERN));
    PetscCall(BVRestoreMat(R,&V));

    /* 11. B-orthonormalize preconditioned residuals */
    PetscCall(BVOrthogonalize(R,NULL));

    /* 13-16. B-orthonormalize conjugate directions */
    if (its>1) PetscCall(BVOrthogonalize(P,NULL));

    /* 17-23. Compute symmetric Gram matrices */
    PetscCall(BVSetActiveColumns(Z,0,ctx->bs));
    PetscCall(BVSetActiveColumns(X,0,ctx->bs));
    PetscCall(BVCopy(X,Z));
    PetscCall(BVSetActiveColumns(Z,ctx->bs,2*ctx->bs-ini));
    PetscCall(BVCopy(R,Z));
    if (its>1) {
      PetscCall(BVSetActiveColumns(Z,2*ctx->bs-ini,3*ctx->bs-2*ini));
      PetscCall(BVCopy(P,Z));
    }

    if (its>1) nv = 3*ctx->bs-2*ini;
    else nv = 2*ctx->bs-ini;

    PetscCall(BVSetActiveColumns(Z,0,nv));
    PetscCall(DSSetDimensions(eps->ds,nv,0,0));
    PetscCall(DSGetMat(eps->ds,DS_MAT_A,&M));
    PetscCall(BVMatProject(Z,A,Z,M));
    if (flip) PetscCall(MatScale(M,-1.0));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&M));
    PetscCall(DSGetMat(eps->ds,DS_MAT_B,&M));
    PetscCall(BVMatProject(Z,B,Z,M)); /* covers also the case B=NULL */
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_B,&M));

    /* 24. Solve the generalized eigenvalue problem */
    PetscCall(DSSetState(eps->ds,DS_STATE_RAW));
    PetscCall(DSSolve(eps->ds,eigr,NULL));
    PetscCall(DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL));
    PetscCall(DSSynchronize(eps->ds,eigr,NULL));
    for (j=0;j<nv;j++) if (locked+j<eps->ncv) eps->eigr[locked+j] = flip? -eigr[j]: eigr[j];
    PetscCall(DSVectors(eps->ds,DS_MAT_X,NULL,NULL));

    /* 25-33. Compute Ritz vectors */
    PetscCall(DSGetMat(eps->ds,DS_MAT_X,&M));
    PetscCall(BVSetActiveColumns(Z,ctx->bs,nv));
    if (ctx->lock) PetscCall(BVSetActiveColumns(P,0,ctx->bs));
    PetscCall(BVMult(P,1.0,0.0,Z,M));
    PetscCall(BVCopy(P,X));
    if (ctx->lock) PetscCall(BVSetActiveColumns(P,nconv,ctx->bs));
    PetscCall(BVSetActiveColumns(Z,0,ctx->bs));
    PetscCall(BVMult(X,1.0,1.0,Z,M));
    if (ctx->lock) PetscCall(BVSetActiveColumns(X,nconv,ctx->bs));
    PetscCall(BVMatMult(X,A,AX));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_X,&M));
  }

diverged:
  eps->its += its;

  if (flip) sc->comparison = SlepcCompareLargestReal;
  PetscCall(PetscFree(eigr));
  PetscCall(MatDestroy(&W));
  if (V) PetscCall(BVRestoreMat(R,&V)); /* only needed when goto diverged is reached */
  PetscCall(BVDestroy(&Z));
  PetscCall(BVDestroy(&X));
  PetscCall(BVDestroy(&R));
  PetscCall(BVDestroy(&P));
  PetscCall(BVDestroy(&AX));
  if (B) PetscCall(BVDestroy(&BX));
  if (nc>0 || eps->nev>ctx->bs-ctx->guard) PetscCall(BVDestroy(&Y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGSetBlockSize_LOBPCG(EPS eps,PetscInt bs)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT || bs == PETSC_DECIDE) bs = 0;
  else PetscCheck(bs>0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid block size %" PetscInt_FMT,bs);
  if (ctx->bs != bs) {
    ctx->bs = bs;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGSetBlockSize - Sets the block size of the LOBPCG method.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs  - the block size

   Options Database Key:
.  -eps_lobpcg_blocksize - Sets the block size

   Level: advanced

.seealso: EPSLOBPCGGetBlockSize()
@*/
PetscErrorCode EPSLOBPCGSetBlockSize(EPS eps,PetscInt bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  PetscTryMethod(eps,"EPSLOBPCGSetBlockSize_C",(EPS,PetscInt),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGGetBlockSize_LOBPCG(EPS eps,PetscInt *bs)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *bs = ctx->bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGGetBlockSize - Gets the block size used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - the block size

   Level: advanced

.seealso: EPSLOBPCGSetBlockSize()
@*/
PetscErrorCode EPSLOBPCGGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(bs,2);
  PetscUseMethod(eps,"EPSLOBPCGGetBlockSize_C",(EPS,PetscInt*),(eps,bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGSetRestart_LOBPCG(EPS eps,PetscReal restart)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  if (restart==(PetscReal)PETSC_DEFAULT || restart==(PetscReal)PETSC_DECIDE) restart = 0.9;
  PetscCheck(restart>=0.1 && restart<=1.0,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The restart argument %g must be in the range [0.1,1.0]",(double)restart);
  if (restart != ctx->restart) {
    ctx->restart = restart;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGSetRestart - Sets the restart parameter for the LOBPCG method.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  restart - the percentage of the block of vectors to force a restart

   Options Database Key:
.  -eps_lobpcg_restart - Sets the restart parameter

   Notes:
   The meaning of this parameter is the proportion of vectors within the
   current block iterate that must have converged in order to force a
   restart with hard locking.
   Allowed values are in the range [0.1,1.0]. The default is 0.9.

   Level: advanced

.seealso: EPSLOBPCGGetRestart()
@*/
PetscErrorCode EPSLOBPCGSetRestart(EPS eps,PetscReal restart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,restart,2);
  PetscTryMethod(eps,"EPSLOBPCGSetRestart_C",(EPS,PetscReal),(eps,restart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGGetRestart_LOBPCG(EPS eps,PetscReal *restart)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *restart = ctx->restart;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGGetRestart - Gets the restart parameter used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  restart - the restart parameter

   Level: advanced

.seealso: EPSLOBPCGSetRestart()
@*/
PetscErrorCode EPSLOBPCGGetRestart(EPS eps,PetscReal *restart)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(restart,2);
  PetscUseMethod(eps,"EPSLOBPCGGetRestart_C",(EPS,PetscReal*),(eps,restart));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGSetLocking_LOBPCG(EPS eps,PetscBool lock)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGSetLocking - Choose between locking and non-locking variants of
   the LOBPCG method.

   Logically Collective

   Input Parameters:
+  eps  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -eps_lobpcg_locking - Sets the locking flag

   Notes:
   This flag refers to soft locking (converged vectors within the current
   block iterate), since hard locking is always used (when nev is larger
   than the block size).

   Level: advanced

.seealso: EPSLOBPCGGetLocking()
@*/
PetscErrorCode EPSLOBPCGSetLocking(EPS eps,PetscBool lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,lock,2);
  PetscTryMethod(eps,"EPSLOBPCGSetLocking_C",(EPS,PetscBool),(eps,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLOBPCGGetLocking_LOBPCG(EPS eps,PetscBool *lock)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLOBPCGGetLocking - Gets the locking flag used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: EPSLOBPCGSetLocking()
@*/
PetscErrorCode EPSLOBPCGGetLocking(EPS eps,PetscBool *lock)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(lock,2);
  PetscUseMethod(eps,"EPSLOBPCGGetLocking_C",(EPS,PetscBool*),(eps,lock));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_LOBPCG(EPS eps,PetscViewer viewer)
{
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  block size %" PetscInt_FMT "\n",ctx->bs));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  restart parameter=%g (using %" PetscInt_FMT " guard vectors)\n",(double)ctx->restart,ctx->guard));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  soft locking %sactivated\n",ctx->lock?"":"de"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_LOBPCG(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      lock,flg;
  PetscInt       bs;
  PetscReal      restart;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS LOBPCG Options");

    PetscCall(PetscOptionsInt("-eps_lobpcg_blocksize","Block size","EPSLOBPCGSetBlockSize",20,&bs,&flg));
    if (flg) PetscCall(EPSLOBPCGSetBlockSize(eps,bs));

    PetscCall(PetscOptionsReal("-eps_lobpcg_restart","Percentage of the block of vectors to force a restart","EPSLOBPCGSetRestart",0.5,&restart,&flg));
    if (flg) PetscCall(EPSLOBPCGSetRestart(eps,restart));

    PetscCall(PetscOptionsBool("-eps_lobpcg_locking","Choose between locking and non-locking variants","EPSLOBPCGSetLocking",PETSC_TRUE,&lock,&flg));
    if (flg) PetscCall(EPSLOBPCGSetLocking(eps,lock));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_LOBPCG(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetBlockSize_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetLocking_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetLocking_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_LOBPCG(EPS eps)
{
  EPS_LOBPCG     *lobpcg;

  PetscFunctionBegin;
  PetscCall(PetscNew(&lobpcg));
  eps->data = (void*)lobpcg;
  lobpcg->lock = PETSC_TRUE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_LOBPCG;
  eps->ops->setup          = EPSSetUp_LOBPCG;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_LOBPCG;
  eps->ops->destroy        = EPSDestroy_LOBPCG;
  eps->ops->view           = EPSView_LOBPCG;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetBlockSize_C",EPSLOBPCGSetBlockSize_LOBPCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetBlockSize_C",EPSLOBPCGGetBlockSize_LOBPCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetRestart_C",EPSLOBPCGSetRestart_LOBPCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetRestart_C",EPSLOBPCGGetRestart_LOBPCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetLocking_C",EPSLOBPCGSetLocking_LOBPCG));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetLocking_C",EPSLOBPCGGetLocking_LOBPCG));
  PetscFunctionReturn(PETSC_SUCCESS);
}
