/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cross"

   Method: Uses a Hermitian eigensolver for A^T*A
*/

#include <slepc/private/svdimpl.h>                /*I "slepcsvd.h" I*/

typedef struct {
  PetscBool explicitmatrix;
  EPS       eps;
  PetscBool usereps;
  Mat       C,D;
} SVD_CROSS;

typedef struct {
  Mat       A,AT;
  Vec       w,diag,omega;
  PetscBool swapped;
} SVD_CROSS_SHELL;

static PetscErrorCode MatMult_Cross(Mat B,Vec x,Vec y)
{
  SVD_CROSS_SHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatMult(ctx->A,x,ctx->w));
  if (ctx->omega && !ctx->swapped) PetscCall(VecPointwiseMult(ctx->w,ctx->w,ctx->omega));
  PetscCall(MatMult(ctx->AT,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatGetDiagonal_Cross(Mat B,Vec d)
{
  SVD_CROSS_SHELL   *ctx;
  PetscMPIInt       len;
  PetscInt          N,n,i,j,start,end,ncols;
  PetscScalar       *work1,*work2,*diag;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  if (!ctx->diag) {
    /* compute diagonal from rows and store in ctx->diag */
    PetscCall(VecDuplicate(d,&ctx->diag));
    PetscCall(MatGetSize(ctx->A,NULL,&N));
    PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
    PetscCall(PetscCalloc2(N,&work1,N,&work2));
    if (ctx->swapped) {
      PetscCall(MatGetOwnershipRange(ctx->AT,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->AT,i,&ncols,NULL,&vals));
        for (j=0;j<ncols;j++) work1[i] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->AT,i,&ncols,NULL,&vals));
      }
    } else {
      PetscCall(MatGetOwnershipRange(ctx->A,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->A,i,&ncols,&cols,&vals));
        for (j=0;j<ncols;j++) work1[cols[j]] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->A,i,&ncols,&cols,&vals));
      }
    }
    PetscCall(PetscMPIIntCast(N,&len));
    PetscCallMPI(MPIU_Allreduce(work1,work2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)B)));
    PetscCall(VecGetOwnershipRange(ctx->diag,&start,&end));
    PetscCall(VecGetArrayWrite(ctx->diag,&diag));
    for (i=start;i<end;i++) diag[i-start] = work2[i];
    PetscCall(VecRestoreArrayWrite(ctx->diag,&diag));
    PetscCall(PetscFree2(work1,work2));
  }
  PetscCall(VecCopy(ctx->diag,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_Cross(Mat B)
{
  SVD_CROSS_SHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(VecDestroy(&ctx->diag));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCrossGetProductMat(SVD svd,Mat A,Mat AT,Mat *C)
{
  SVD_CROSS       *cross = (SVD_CROSS*)svd->data;
  SVD_CROSS_SHELL *ctx;
  PetscInt        n;
  VecType         vtype;
  Mat             B;

  PetscFunctionBegin;
  if (cross->explicitmatrix) {
    if (!svd->ishyperbolic || svd->swapped) B = (!svd->expltrans && svd->swapped)? AT: A;
    else {  /* duplicate A and scale by signature */
      PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&B));
      PetscCall(MatDiagonalScale(B,svd->omega,NULL));
    }
    if (svd->expltrans) {  /* explicit transpose */
      PetscCall(MatProductCreate(AT,B,NULL,C));
      PetscCall(MatProductSetType(*C,MATPRODUCT_AB));
    } else {  /* implicit transpose */
#if defined(PETSC_USE_COMPLEX)
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Must use explicit transpose with complex scalars");
#else
      if (!svd->swapped) {
        PetscCall(MatProductCreate(A,B,NULL,C));
        PetscCall(MatProductSetType(*C,MATPRODUCT_AtB));
      } else {
        PetscCall(MatProductCreate(B,AT,NULL,C));
        PetscCall(MatProductSetType(*C,MATPRODUCT_ABt));
      }
#endif
    }
    PetscCall(MatProductSetFromOptions(*C));
    PetscCall(MatProductSymbolic(*C));
    PetscCall(MatProductNumeric(*C));
    if (svd->ishyperbolic && !svd->swapped) PetscCall(MatDestroy(&B));
  } else {
    PetscCall(PetscNew(&ctx));
    ctx->A       = A;
    ctx->AT      = AT;
    ctx->omega   = svd->omega;
    ctx->swapped = svd->swapped;
    PetscCall(MatCreateVecs(A,NULL,&ctx->w));
    PetscCall(MatGetLocalSize(A,NULL,&n));
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)svd),n,n,PETSC_DETERMINE,PETSC_DETERMINE,(void*)ctx,C));
    PetscCall(MatShellSetOperation(*C,MATOP_MULT,(void(*)(void))MatMult_Cross));
    if (!svd->ishyperbolic || svd->swapped) PetscCall(MatShellSetOperation(*C,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Cross));
    PetscCall(MatShellSetOperation(*C,MATOP_DESTROY,(void(*)(void))MatDestroy_Cross));
    PetscCall(MatGetVecType(A,&vtype));
    PetscCall(MatSetVecType(*C,vtype));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Convergence test relative to the norm of R (used in GSVD only) */
static PetscErrorCode EPSConv_Cross(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  SVD svd = (SVD)ctx;

  PetscFunctionBegin;
  *errest = res/PetscMax(svd->nrma,svd->nrmb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetUp_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  ST             st;
  PetscBool      trackall,issinv,isks;
  EPSProblemType ptype;
  EPSWhich       which;
  Mat            Omega;
  MatType        Atype;
  PetscInt       n,N;

  PetscFunctionBegin;
  if (svd->nsv==0 && svd->stop!=SVD_STOP_THRESHOLD) svd->nsv = 1;
  if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
  PetscCall(MatDestroy(&cross->C));
  PetscCall(MatDestroy(&cross->D));
  PetscCall(SVDCrossGetProductMat(svd,svd->A,svd->AT,&cross->C));
  if (svd->isgeneralized) {
    PetscCall(SVDCrossGetProductMat(svd,svd->B,svd->BT,&cross->D));
    PetscCall(EPSSetOperators(cross->eps,cross->C,cross->D));
    PetscCall(EPSGetProblemType(cross->eps,&ptype));
    if (!ptype) PetscCall(EPSSetProblemType(cross->eps,EPS_GHEP));
  } else if (svd->ishyperbolic && svd->swapped) {
    PetscCall(MatGetType(svd->OP,&Atype));
    PetscCall(MatGetSize(svd->A,NULL,&N));
    PetscCall(MatGetLocalSize(svd->A,NULL,&n));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Omega));
    PetscCall(MatSetSizes(Omega,n,n,N,N));
    PetscCall(MatSetType(Omega,Atype));
    PetscCall(MatDiagonalSet(Omega,svd->omega,INSERT_VALUES));
    PetscCall(EPSSetOperators(cross->eps,cross->C,Omega));
    PetscCall(EPSSetProblemType(cross->eps,EPS_GHIEP));
    PetscCall(MatDestroy(&Omega));
  } else {
    PetscCall(EPSSetOperators(cross->eps,cross->C,NULL));
    PetscCall(EPSSetProblemType(cross->eps,EPS_HEP));
  }
  if (!cross->usereps) {
    PetscCall(EPSGetST(cross->eps,&st));
    PetscCall(PetscObjectTypeCompare((PetscObject)st,STSINVERT,&issinv));
    PetscCall(PetscObjectTypeCompare((PetscObject)cross->eps,EPSKRYLOVSCHUR,&isks));
    if (svd->isgeneralized && svd->which==SVD_SMALLEST) {
      if (cross->explicitmatrix && isks && !issinv) {  /* default to shift-and-invert */
        PetscCall(STSetType(st,STSINVERT));
        PetscCall(EPSSetTarget(cross->eps,0.0));
        which = EPS_TARGET_REAL;
      } else which = issinv?EPS_TARGET_REAL:EPS_SMALLEST_REAL;
    } else {
      if (issinv) which = EPS_TARGET_MAGNITUDE;
      else if (svd->ishyperbolic) which = svd->which==SVD_LARGEST?EPS_LARGEST_MAGNITUDE:EPS_SMALLEST_MAGNITUDE;
      else which = svd->which==SVD_LARGEST?EPS_LARGEST_MAGNITUDE:EPS_SMALLEST_MAGNITUDE;
    }
    PetscCall(EPSSetWhichEigenpairs(cross->eps,which));
    PetscCall(EPSSetDimensions(cross->eps,svd->nsv?svd->nsv:PETSC_CURRENT,svd->ncv,svd->mpd));
    if (svd->stop==SVD_STOP_THRESHOLD) PetscCall(EPSSetThreshold(cross->eps,svd->thres*svd->thres,svd->threlative));
    PetscCall(EPSSetTolerances(cross->eps,svd->tol==(PetscReal)PETSC_DETERMINE?SLEPC_DEFAULT_TOL/10.0:svd->tol,svd->max_it));
    switch (svd->conv) {
    case SVD_CONV_ABS:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_ABS));break;
    case SVD_CONV_REL:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_REL));break;
    case SVD_CONV_NORM:
      if (svd->isgeneralized) {
        if (!svd->nrma) PetscCall(MatNorm(svd->OP,NORM_INFINITY,&svd->nrma));
        if (!svd->nrmb) PetscCall(MatNorm(svd->OPb,NORM_INFINITY,&svd->nrmb));
        PetscCall(EPSSetConvergenceTestFunction(cross->eps,EPSConv_Cross,svd,NULL));
      } else {
        PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_NORM));break;
      }
      break;
    case SVD_CONV_MAXIT:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Maxit convergence test not supported in this solver");
    case SVD_CONV_USER:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined convergence test not supported in this solver");
    }
  }
  /* Transfer the trackall option from svd to eps */
  PetscCall(SVDGetTrackAll(svd,&trackall));
  PetscCall(EPSSetTrackAll(cross->eps,trackall));
  /* Transfer the initial space from svd to eps */
  if (svd->nini<0) {
    PetscCall(EPSSetInitialSpace(cross->eps,-svd->nini,svd->IS));
    PetscCall(SlepcBasisDestroy_Private(&svd->nini,&svd->IS));
  }
  PetscCall(EPSSetUp(cross->eps));
  PetscCall(EPSGetDimensions(cross->eps,NULL,&svd->ncv,&svd->mpd));
  PetscCall(EPSGetTolerances(cross->eps,NULL,&svd->max_it));
  if (svd->tol==(PetscReal)PETSC_DETERMINE) svd->tol = SLEPC_DEFAULT_TOL;

  svd->leftbasis = PETSC_FALSE;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSolve_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscInt       i;
  PetscScalar    lambda;
  PetscReal      sigma;

  PetscFunctionBegin;
  PetscCall(EPSSolve(cross->eps));
  PetscCall(EPSGetConverged(cross->eps,&svd->nconv));
  PetscCall(EPSGetIterationNumber(cross->eps,&svd->its));
  PetscCall(EPSGetConvergedReason(cross->eps,(EPSConvergedReason*)&svd->reason));
  for (i=0;i<svd->nconv;i++) {
    PetscCall(EPSGetEigenvalue(cross->eps,i,&lambda,NULL));
    sigma = PetscRealPart(lambda);
    if (svd->ishyperbolic) svd->sigma[i] = PetscSqrtReal(PetscAbsReal(sigma));
    else {
      PetscCheck(sigma>-10*PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)svd),PETSC_ERR_FP,"Negative eigenvalue computed by EPS: %g",(double)sigma);
      if (sigma<0.0) {
        PetscCall(PetscInfo(svd,"Negative eigenvalue computed by EPS: %g, resetting to 0\n",(double)sigma));
        sigma = 0.0;
      }
      svd->sigma[i] = PetscSqrtReal(sigma);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDComputeVectors_Cross(SVD svd)
{
  SVD_CROSS         *cross = (SVD_CROSS*)svd->data;
  PetscInt          i,mloc,ploc;
  Vec               u,v,x,uv,w,omega2=NULL;
  Mat               Omega;
  PetscScalar       *dst,alpha,lambda,*varray;
  const PetscScalar *src;
  PetscReal         nrm;

  PetscFunctionBegin;
  if (svd->isgeneralized) {
    PetscCall(MatCreateVecs(svd->A,NULL,&u));
    PetscCall(VecGetLocalSize(u,&mloc));
    PetscCall(MatCreateVecs(svd->B,NULL,&v));
    PetscCall(VecGetLocalSize(v,&ploc));
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&x));
      PetscCall(EPSGetEigenpair(cross->eps,i,&lambda,NULL,x,NULL));
      PetscCall(MatMult(svd->A,x,u));     /* u_i*c_i/alpha = A*x_i */
      PetscCall(VecNormalize(u,NULL));
      PetscCall(MatMult(svd->B,x,v));     /* v_i*s_i/alpha = B*x_i */
      PetscCall(VecNormalize(v,&nrm));    /* ||v||_2 = s_i/alpha   */
      alpha = 1.0/(PetscSqrtReal(1.0+PetscRealPart(lambda))*nrm);    /* alpha=s_i/||v||_2 */
      PetscCall(VecScale(x,alpha));
      PetscCall(BVRestoreColumn(svd->V,i,&x));
      /* copy [u;v] to U[i] */
      PetscCall(BVGetColumn(svd->U,i,&uv));
      PetscCall(VecGetArrayWrite(uv,&dst));
      PetscCall(VecGetArrayRead(u,&src));
      PetscCall(PetscArraycpy(dst,src,mloc));
      PetscCall(VecRestoreArrayRead(u,&src));
      PetscCall(VecGetArrayRead(v,&src));
      PetscCall(PetscArraycpy(dst+mloc,src,ploc));
      PetscCall(VecRestoreArrayRead(v,&src));
      PetscCall(VecRestoreArrayWrite(uv,&dst));
      PetscCall(BVRestoreColumn(svd->U,i,&uv));
    }
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&u));
  } else if (svd->ishyperbolic && svd->swapped) {  /* was solved as GHIEP, set u=Omega*u and normalize */
    PetscCall(EPSGetOperators(cross->eps,NULL,&Omega));
    PetscCall(MatCreateVecs(Omega,&w,NULL));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,svd->ncv,&omega2));
    PetscCall(VecGetArrayWrite(omega2,&varray));
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&v));
      PetscCall(EPSGetEigenvector(cross->eps,i,v,NULL));
      PetscCall(MatMult(Omega,v,w));
      PetscCall(VecDot(v,w,&alpha));
      svd->sign[i] = PetscSign(PetscRealPart(alpha));
      varray[i] = svd->sign[i];
      alpha = 1.0/PetscSqrtScalar(PetscAbsScalar(alpha));
      PetscCall(VecScale(w,alpha));
      PetscCall(VecCopy(w,v));
      PetscCall(BVRestoreColumn(svd->V,i,&v));
    }
    PetscCall(BVSetSignature(svd->V,omega2));
    PetscCall(VecRestoreArrayWrite(omega2,&varray));
    PetscCall(VecDestroy(&omega2));
    PetscCall(VecDestroy(&w));
    PetscCall(SVDComputeVectors_Left(svd));
  } else {
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&v));
      PetscCall(EPSGetEigenvector(cross->eps,i,v,NULL));
      PetscCall(BVRestoreColumn(svd->V,i,&v));
    }
    PetscCall(SVDComputeVectors_Left(svd));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSMonitor_Cross(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i,ncv;
  SVD            svd = (SVD)ctx;
  SVD_CROSS      *cross;
  PetscScalar    er,ei;
  ST             st;

  PetscFunctionBegin;
  if (svd->stop==SVD_STOP_THRESHOLD) {
    cross = (SVD_CROSS*)svd->data;
    PetscCall(EPSGetDimensions(cross->eps,NULL,&ncv,NULL));
    if (ncv!=svd->ncv) {  /* reallocate */
      PetscCall(SVDReallocateSolution(svd,ncv));
      for (i=svd->ncv;i<ncv;i++) svd->perm[i] = i;
      svd->ncv = ncv;
    }
  }
  PetscCall(EPSGetST(eps,&st));
  for (i=0;i<PetscMin(nest,svd->ncv);i++) {
    er = eigr[i]; ei = eigi[i];
    PetscCall(STBackTransform(st,1,&er,&ei));
    svd->sigma[i] = PetscSqrtReal(PetscAbsReal(PetscRealPart(er)));
    svd->errest[i] = errest[i];
  }
  PetscCall(SVDMonitor(svd,its,nconv,svd->sigma,svd->errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDSetFromOptions_Cross(SVD svd,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      set,val;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  ST             st;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD Cross Options");

    PetscCall(PetscOptionsBool("-svd_cross_explicitmatrix","Use cross explicit matrix","SVDCrossSetExplicitMatrix",cross->explicitmatrix,&val,&set));
    if (set) PetscCall(SVDCrossSetExplicitMatrix(svd,val));

  PetscOptionsHeadEnd();

  if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
  if (!cross->explicitmatrix && !cross->usereps) {
    /* use as default an ST with shell matrix and Jacobi */
    PetscCall(EPSGetST(cross->eps,&st));
    PetscCall(STSetMatMode(st,ST_MATMODE_SHELL));
  }
  PetscCall(EPSSetFromOptions(cross->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCrossSetExplicitMatrix_Cross(SVD svd,PetscBool explicitmatrix)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  if (cross->explicitmatrix != explicitmatrix) {
    cross->explicitmatrix = explicitmatrix;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCrossSetExplicitMatrix - Indicate if the eigensolver operator A^T*A must
   be computed explicitly.

   Logically Collective

   Input Parameters:
+  svd         - singular value solver
-  explicitmat - boolean flag indicating if A^T*A is built explicitly

   Options Database Key:
.  -svd_cross_explicitmatrix <boolean> - Indicates the boolean flag

   Level: advanced

.seealso: SVDCrossGetExplicitMatrix()
@*/
PetscErrorCode SVDCrossSetExplicitMatrix(SVD svd,PetscBool explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,explicitmat,2);
  PetscTryMethod(svd,"SVDCrossSetExplicitMatrix_C",(SVD,PetscBool),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCrossGetExplicitMatrix_Cross(SVD svd,PetscBool *explicitmat)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  *explicitmat = cross->explicitmatrix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCrossGetExplicitMatrix - Returns the flag indicating if A^T*A is built explicitly.

   Not Collective

   Input Parameter:
.  svd  - singular value solver

   Output Parameter:
.  explicitmat - the mode flag

   Level: advanced

.seealso: SVDCrossSetExplicitMatrix()
@*/
PetscErrorCode SVDCrossGetExplicitMatrix(SVD svd,PetscBool *explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(explicitmat,2);
  PetscUseMethod(svd,"SVDCrossGetExplicitMatrix_C",(SVD,PetscBool*),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCrossSetEPS_Cross(SVD svd,EPS eps)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&cross->eps));
  cross->eps     = eps;
  cross->usereps = PETSC_TRUE;
  svd->state     = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCrossSetEPS - Associate an eigensolver object (EPS) to the
   singular value solver.

   Collective

   Input Parameters:
+  svd - singular value solver
-  eps - the eigensolver object

   Level: advanced

.seealso: SVDCrossGetEPS()
@*/
PetscErrorCode SVDCrossSetEPS(SVD svd,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(svd,1,eps,2);
  PetscTryMethod(svd,"SVDCrossSetEPS_C",(SVD,EPS),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDCrossGetEPS_Cross(SVD svd,EPS *eps)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  if (!cross->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)svd),&cross->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)cross->eps,(PetscObject)svd,1));
    PetscCall(EPSSetOptionsPrefix(cross->eps,((PetscObject)svd)->prefix));
    PetscCall(EPSAppendOptionsPrefix(cross->eps,"svd_cross_"));
    PetscCall(PetscObjectSetOptions((PetscObject)cross->eps,((PetscObject)svd)->options));
    PetscCall(EPSSetWhichEigenpairs(cross->eps,EPS_LARGEST_MAGNITUDE));
    PetscCall(EPSMonitorSet(cross->eps,EPSMonitor_Cross,svd,NULL));
  }
  *eps = cross->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   SVDCrossGetEPS - Retrieve the eigensolver object (EPS) associated
   to the singular value solver.

   Collective

   Input Parameter:
.  svd - singular value solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: SVDCrossSetEPS()
@*/
PetscErrorCode SVDCrossGetEPS(SVD svd,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(svd,"SVDCrossGetEPS_C",(SVD,EPS*),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDView_Cross(SVD svd,PetscViewer viewer)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %s matrix\n",cross->explicitmatrix?"explicit":"implicit"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(cross->eps,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDReset_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSReset(cross->eps));
  PetscCall(MatDestroy(&cross->C));
  PetscCall(MatDestroy(&cross->D));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SVDDestroy_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSDestroy(&cross->eps));
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Cross(SVD svd)
{
  SVD_CROSS      *cross;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cross));
  svd->data = (void*)cross;

  svd->ops->solve          = SVDSolve_Cross;
  svd->ops->solveg         = SVDSolve_Cross;
  svd->ops->solveh         = SVDSolve_Cross;
  svd->ops->setup          = SVDSetUp_Cross;
  svd->ops->setfromoptions = SVDSetFromOptions_Cross;
  svd->ops->destroy        = SVDDestroy_Cross;
  svd->ops->reset          = SVDReset_Cross;
  svd->ops->view           = SVDView_Cross;
  svd->ops->computevectors = SVDComputeVectors_Cross;
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",SVDCrossSetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",SVDCrossGetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",SVDCrossSetExplicitMatrix_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",SVDCrossGetExplicitMatrix_Cross));
  PetscFunctionReturn(PETSC_SUCCESS);
}
