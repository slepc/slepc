/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lyapii"

   Method: Lyapunov inverse iteration

   Algorithm:

       Lyapunov inverse iteration using LME solvers

   References:

       [1] H.C. Elman and M. Wu, "Lyapunov inverse iteration for computing a
           few rightmost eigenvalues of large generalized eigenvalue problems",
           SIAM J. Matrix Anal. Appl. 34(4):1685-1707, 2013.

       [2] K. Meerbergen and A. Spence, "Inverse iteration for purely imaginary
           eigenvalues with application to the detection of Hopf bifurcations in
           large-scale problems", SIAM J. Matrix Anal. Appl. 31:1982-1999, 2010.
*/

#include <slepc/private/epsimpl.h>          /*I "slepceps.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  LME      lme;      /* Lyapunov solver */
  DS       ds;       /* used to compute the SVD for compression */
  PetscInt rkl;      /* prescribed rank for the Lyapunov solver */
  PetscInt rkc;      /* the compressed rank, cannot be larger than rkl */
} EPS_LYAPII;

typedef struct {
  Mat      S;        /* the operator matrix, S=A^{-1}*B */
  BV       Q;        /* orthogonal basis of converged eigenvectors */
} EPS_LYAPII_MATSHELL;

typedef struct {
  Mat      S;        /* the matrix from which the implicit operator is built */
  PetscInt n;        /* the size of matrix S, the operator is nxn */
  LME      lme;      /* dummy LME object */
#if defined(PETSC_USE_COMPLEX)
  Mat      A,B,F;
  Vec      w;
#endif
} EPS_EIG_MATSHELL;

static PetscErrorCode EPSSetUp_LyapII(EPS eps)
{
  PetscRandom    rand;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  EPSCheckSinvert(eps);
  EPSCheckNotStructured(eps);
  if (eps->nev==0) eps->nev = 1;
  if (eps->ncv!=PETSC_DETERMINE) {
    PetscCheck(eps->ncv>=eps->nev+1,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must be at least nev+1");
  } else eps->ncv = eps->nev+1;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(1000*eps->nev,100*eps->n);
  if (!eps->which) eps->which=EPS_LARGEST_REAL;
  PetscCheck(eps->which==EPS_LARGEST_REAL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver supports only largest real eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_THRESHOLD | EPS_FEATURE_TWOSIDED);

  if (!ctx->rkc) ctx->rkc = 10;
  if (!ctx->rkl) ctx->rkl = 3*ctx->rkc;
  if (!ctx->lme) PetscCall(EPSLyapIIGetLME(eps,&ctx->lme));
  PetscCall(LMESetProblemType(ctx->lme,LME_LYAPUNOV));
  PetscCall(LMESetErrorIfNotConverged(ctx->lme,PETSC_TRUE));

  if (!ctx->ds) {
    PetscCall(DSCreate(PetscObjectComm((PetscObject)eps),&ctx->ds));
    PetscCall(DSSetType(ctx->ds,DSSVD));
  }
  PetscCall(DSAllocate(ctx->ds,ctx->rkl));

  PetscCall(DSSetType(eps->ds,DSNHEP));
  PetscCall(DSAllocate(eps->ds,eps->ncv));

  PetscCall(EPSAllocateSolution(eps,0));
  PetscCall(BVGetRandomContext(eps->V,&rand));  /* make sure the random context is available when duplicating */
  PetscCall(EPSSetWorkVecs(eps,3));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_EPSLyapIIOperator(Mat M,Vec x,Vec r)
{
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  PetscCall(MatMult(matctx->S,x,r));
  PetscCall(BVOrthogonalizeVec(matctx->Q,r,NULL,NULL,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_EPSLyapIIOperator(Mat M)
{
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  PetscCall(MatDestroy(&matctx->S));
  PetscCall(PetscFree(matctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_EigOperator(Mat M,Vec x,Vec y)
{
  EPS_EIG_MATSHELL  *matctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt          n,lds;
  PetscScalar       *Y,*C,zero=0.0,done=1.0,dtwo=2.0;
  const PetscScalar *S,*X;
  PetscBLASInt      n_,lds_;
#endif

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));

#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatMult(matctx->B,x,matctx->w));
  PetscCall(MatSolve(matctx->F,matctx->w,y));
#else
  PetscCall(VecGetArrayRead(x,&X));
  PetscCall(VecGetArray(y,&Y));
  PetscCall(MatDenseGetArrayRead(matctx->S,&S));
  PetscCall(MatDenseGetLDA(matctx->S,&lds));

  n = matctx->n;
  PetscCall(PetscCalloc1(n*n,&C));
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(lds,&lds_));

  /* C = 2*S*X*S.' */
  PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&dtwo,S,&lds_,X,&n_,&zero,Y,&n_));
  PetscCallBLAS("BLASgemm",BLASgemm_("N","T",&n_,&n_,&n_,&done,Y,&n_,S,&lds_,&zero,C,&n_));

  /* Solve S*Y + Y*S' = -C */
  PetscCall(LMEDenseLyapunov(matctx->lme,n,(PetscScalar*)S,lds,C,n,Y,n));

  PetscCall(PetscFree(C));
  PetscCall(VecRestoreArrayRead(x,&X));
  PetscCall(VecRestoreArray(y,&Y));
  PetscCall(MatDenseRestoreArrayRead(matctx->S,&S));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_EigOperator(Mat M)
{
  EPS_EIG_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatDestroy(&matctx->A));
  PetscCall(MatDestroy(&matctx->B));
  PetscCall(MatDestroy(&matctx->F));
  PetscCall(VecDestroy(&matctx->w));
#else
  PetscCall(MatDestroy(&matctx->S));
#endif
  PetscCall(PetscFree(matctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EV2x2: solve the eigenproblem for a 2x2 matrix M
 */
static PetscErrorCode EV2x2(PetscScalar *M,PetscInt ld,PetscScalar *wr,PetscScalar *wi,PetscScalar *vec)
{
  PetscBLASInt   lwork=10,ld_;
  PetscScalar    work[10];
  PetscBLASInt   two=2,info;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      rwork[6];
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(ld,&ld_));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgeev",LAPACKgeev_("N","V",&two,M,&ld_,wr,wi,NULL,&ld_,vec,&ld_,work,&lwork,&info));
#else
  PetscCallBLAS("LAPACKgeev",LAPACKgeev_("N","V",&two,M,&ld_,wr,NULL,&ld_,vec,&ld_,work,&lwork,rwork,&info));
#endif
  SlepcCheckLapackInfo("geev",info);
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   LyapIIBuildRHS: prepare the right-hand side of the Lyapunov equation SY + YS' = -2*S*Z*S'
   in factored form:
      if (V)  U=sqrt(2)*S*V    (uses 1 work vector)
      else    U=sqrt(2)*S*U    (uses 2 work vectors)
   where U,V are assumed to have rk columns.
 */
static PetscErrorCode LyapIIBuildRHS(Mat S,PetscInt rk,Mat U,BV V,Vec *work)
{
  PetscScalar    *array,*uu;
  PetscInt       i,nloc;
  Vec            v,u=work[0];

  PetscFunctionBegin;
  PetscCall(MatGetLocalSize(U,&nloc,NULL));
  for (i=0;i<rk;i++) {
    PetscCall(MatDenseGetColumn(U,i,&array));
    if (V) PetscCall(BVGetColumn(V,i,&v));
    else {
      v = work[1];
      PetscCall(VecPlaceArray(v,array));
    }
    PetscCall(MatMult(S,v,u));
    if (V) PetscCall(BVRestoreColumn(V,i,&v));
    else PetscCall(VecResetArray(v));
    PetscCall(VecScale(u,PETSC_SQRT2));
    PetscCall(VecGetArray(u,&uu));
    PetscCall(PetscArraycpy(array,uu,nloc));
    PetscCall(VecRestoreArray(u,&uu));
    PetscCall(MatDenseRestoreColumn(U,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   LyapIIBuildEigenMat: create shell matrix Op=A\B with A = kron(I,S)+kron(S,I), B = -2*kron(S,S)
   where S is a sequential square dense matrix of order n.
   v0 is the initial vector, should have the form v0 = w*w' (for instance 1*1')
 */
static PetscErrorCode LyapIIBuildEigenMat(LME lme,Mat S,Mat *Op,Vec *v0)
{
  PetscInt          n,m;
  PetscBool         create=PETSC_FALSE;
  EPS_EIG_MATSHELL  *matctx;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar       theta,*aa,*bb;
  const PetscScalar *ss;
  PetscInt          i,j,f,c,off,ld,lds;
  IS                perm;
#endif

  PetscFunctionBegin;
  PetscCall(MatGetSize(S,&n,NULL));
  if (!*Op) create=PETSC_TRUE;
  else {
    PetscCall(MatGetSize(*Op,&m,NULL));
    if (m!=n*n) create=PETSC_TRUE;
  }
  if (create) {
    PetscCall(MatDestroy(Op));
    PetscCall(VecDestroy(v0));
    PetscCall(PetscNew(&matctx));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n*n,n*n,NULL,&matctx->A));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n*n,n*n,NULL,&matctx->B));
    PetscCall(MatCreateVecs(matctx->A,NULL,&matctx->w));
#else
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&matctx->S));
#endif
    PetscCall(MatCreateShell(PETSC_COMM_SELF,n*n,n*n,PETSC_DETERMINE,PETSC_DETERMINE,matctx,Op));
    PetscCall(MatShellSetOperation(*Op,MATOP_MULT,(void(*)(void))MatMult_EigOperator));
    PetscCall(MatShellSetOperation(*Op,MATOP_DESTROY,(void(*)(void))MatDestroy_EigOperator));
    PetscCall(MatCreateVecs(*Op,NULL,v0));
  } else {
    PetscCall(MatShellGetContext(*Op,&matctx));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(MatZeroEntries(matctx->A));
#endif
  }
#if defined(PETSC_USE_COMPLEX)
  PetscCall(MatDenseGetArray(matctx->A,&aa));
  PetscCall(MatDenseGetArray(matctx->B,&bb));
  PetscCall(MatDenseGetArrayRead(S,&ss));
  PetscCall(MatDenseGetLDA(S,&lds));
  ld = n*n;
  for (f=0;f<n;f++) {
    off = f*n+f*n*ld;
    for (i=0;i<n;i++) for (j=0;j<n;j++) aa[off+i+j*ld] = ss[i+j*lds];
    for (c=0;c<n;c++) {
      off = f*n+c*n*ld;
      theta = ss[f+c*lds];
      for (i=0;i<n;i++) aa[off+i+i*ld] += theta;
      for (i=0;i<n;i++) for (j=0;j<n;j++) bb[off+i+j*ld] = -2*theta*ss[i+j*lds];
    }
  }
  PetscCall(MatDenseRestoreArray(matctx->A,&aa));
  PetscCall(MatDenseRestoreArray(matctx->B,&bb));
  PetscCall(MatDenseRestoreArrayRead(S,&ss));
  PetscCall(ISCreateStride(PETSC_COMM_SELF,n*n,0,1,&perm));
  PetscCall(MatDestroy(&matctx->F));
  PetscCall(MatDuplicate(matctx->A,MAT_COPY_VALUES,&matctx->F));
  PetscCall(MatLUFactor(matctx->F,perm,perm,NULL));
  PetscCall(ISDestroy(&perm));
#else
  PetscCall(MatCopy(S,matctx->S,SAME_NONZERO_PATTERN));
#endif
  matctx->lme = lme;
  matctx->n = n;
  PetscCall(VecSet(*v0,1.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_LyapII(EPS eps)
{
  EPS_LYAPII          *ctx = (EPS_LYAPII*)eps->data;
  PetscInt            i,ldds,rk,nloc,mloc,nv,idx,k;
  Vec                 v,w,z=eps->work[0],v0=NULL;
  Mat                 S,C,Ux[2],Y,Y1,R,U,W,X,Op=NULL;
  BV                  V;
  BVOrthogType        type;
  BVOrthogRefineType  refine;
  PetscScalar         eigr[2],eigi[2],*array,er,ei,*uu,*s,*xx,*aa,pM[4],vec[4];
  PetscReal           eta;
  EPS                 epsrr;
  PetscReal           norm;
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(ctx->ds,&ldds));

  /* Operator for the Lyapunov equation */
  PetscCall(PetscNew(&matctx));
  PetscCall(STGetOperator(eps->st,&matctx->S));
  PetscCall(MatGetLocalSize(matctx->S,&mloc,&nloc));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)eps),mloc,nloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&S));
  matctx->Q = eps->V;
  PetscCall(MatShellSetOperation(S,MATOP_MULT,(void(*)(void))MatMult_EPSLyapIIOperator));
  PetscCall(MatShellSetOperation(S,MATOP_DESTROY,(void(*)(void))MatDestroy_EPSLyapIIOperator));
  PetscCall(LMESetCoefficients(ctx->lme,S,NULL,NULL,NULL));

  /* Right-hand side */
  PetscCall(BVDuplicateResize(eps->V,ctx->rkl,&V));
  PetscCall(BVGetOrthogonalization(V,&type,&refine,&eta,NULL));
  PetscCall(BVSetOrthogonalization(V,type,refine,eta,BV_ORTHOG_BLOCK_TSQR));
  PetscCall(MatCreateDense(PetscObjectComm((PetscObject)eps),eps->nloc,PETSC_DECIDE,PETSC_DECIDE,1,NULL,&Ux[0]));
  PetscCall(MatCreateDense(PetscObjectComm((PetscObject)eps),eps->nloc,PETSC_DECIDE,PETSC_DECIDE,2,NULL,&Ux[1]));
  nv = ctx->rkl;
  PetscCall(PetscMalloc1(nv,&s));

  /* Initialize first column */
  PetscCall(EPSGetStartVector(eps,0,NULL));
  PetscCall(BVGetColumn(eps->V,0,&v));
  PetscCall(BVInsertVec(V,0,v));
  PetscCall(BVRestoreColumn(eps->V,0,&v));
  PetscCall(BVSetActiveColumns(eps->V,0,0));  /* no deflation at the beginning */
  PetscCall(LyapIIBuildRHS(S,1,Ux[0],V,eps->work));
  idx = 0;

  /* EPS for rank reduction */
  PetscCall(EPSCreate(PETSC_COMM_SELF,&epsrr));
  PetscCall(EPSSetOptionsPrefix(epsrr,((PetscObject)eps)->prefix));
  PetscCall(EPSAppendOptionsPrefix(epsrr,"eps_lyapii_"));
  PetscCall(EPSSetDimensions(epsrr,1,PETSC_CURRENT,PETSC_CURRENT));
  PetscCall(EPSSetTolerances(epsrr,PETSC_MACHINE_EPSILON*100,PETSC_CURRENT));

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Matrix for placing the solution of the Lyapunov equation (an alias of V) */
    PetscCall(BVSetActiveColumns(V,0,nv));
    PetscCall(BVGetMat(V,&Y1));
    PetscCall(MatZeroEntries(Y1));
    PetscCall(MatCreateLRC(NULL,Y1,NULL,NULL,&Y));
    PetscCall(LMESetSolution(ctx->lme,Y));

    /* Solve the Lyapunov equation SY + YS' = -2*S*Z*S' */
    PetscCall(MatCreateLRC(NULL,Ux[idx],NULL,NULL,&C));
    PetscCall(LMESetRHS(ctx->lme,C));
    PetscCall(MatDestroy(&C));
    PetscCall(LMESolve(ctx->lme));
    PetscCall(BVRestoreMat(V,&Y1));
    PetscCall(MatDestroy(&Y));

    /* SVD of the solution: [Q,R]=qr(V); [U,Sigma,~]=svd(R) */
    PetscCall(DSSetDimensions(ctx->ds,nv,0,0));
    PetscCall(DSSVDSetDimensions(ctx->ds,nv));
    PetscCall(DSGetMat(ctx->ds,DS_MAT_A,&R));
    PetscCall(BVOrthogonalize(V,R));
    PetscCall(DSRestoreMat(ctx->ds,DS_MAT_A,&R));
    PetscCall(DSSetState(ctx->ds,DS_STATE_RAW));
    PetscCall(DSSolve(ctx->ds,s,NULL));

    /* Determine rank */
    rk = nv;
    for (i=1;i<nv;i++) if (PetscAbsScalar(s[i]/s[0])<PETSC_SQRT_MACHINE_EPSILON) {rk=i; break;}
    PetscCall(PetscInfo(eps,"The computed solution of the Lyapunov equation has rank %" PetscInt_FMT "\n",rk));
    rk = PetscMin(rk,ctx->rkc);
    PetscCall(DSGetMat(ctx->ds,DS_MAT_U,&U));
    PetscCall(BVMultInPlace(V,U,0,rk));
    PetscCall(DSRestoreMat(ctx->ds,DS_MAT_U,&U));
    PetscCall(BVSetActiveColumns(V,0,rk));

    /* Rank reduction */
    PetscCall(DSSetDimensions(ctx->ds,rk,0,0));
    PetscCall(DSSVDSetDimensions(ctx->ds,rk));
    PetscCall(DSGetMat(ctx->ds,DS_MAT_A,&W));
    PetscCall(BVMatProject(V,S,V,W));
    PetscCall(LyapIIBuildEigenMat(ctx->lme,W,&Op,&v0)); /* Op=A\B, A=kron(I,S)+kron(S,I), B=-2*kron(S,S) */
    PetscCall(DSRestoreMat(ctx->ds,DS_MAT_A,&W));
    PetscCall(EPSSetOperators(epsrr,Op,NULL));
    PetscCall(EPSSetInitialSpace(epsrr,1,&v0));
    PetscCall(EPSSolve(epsrr));
    PetscCall(EPSComputeVectors(epsrr));
    /* Copy first eigenvector, vec(A)=x */
    PetscCall(BVGetArray(epsrr->V,&xx));
    PetscCall(DSGetArray(ctx->ds,DS_MAT_A,&aa));
    for (i=0;i<rk;i++) PetscCall(PetscArraycpy(aa+i*ldds,xx+i*rk,rk));
    PetscCall(DSRestoreArray(ctx->ds,DS_MAT_A,&aa));
    PetscCall(BVRestoreArray(epsrr->V,&xx));
    PetscCall(DSSetState(ctx->ds,DS_STATE_RAW));
    /* Compute [U,Sigma,~] = svd(A), its rank should be 1 or 2 */
    PetscCall(DSSolve(ctx->ds,s,NULL));
    if (PetscAbsScalar(s[1]/s[0])<PETSC_SQRT_MACHINE_EPSILON) rk=1;
    else rk = 2;
    PetscCall(PetscInfo(eps,"The eigenvector has rank %" PetscInt_FMT "\n",rk));
    PetscCall(DSGetMat(ctx->ds,DS_MAT_U,&U));
    PetscCall(BVMultInPlace(V,U,0,rk));
    PetscCall(DSRestoreMat(ctx->ds,DS_MAT_U,&U));

    /* Save V in Ux */
    idx = (rk==2)?1:0;
    for (i=0;i<rk;i++) {
      PetscCall(BVGetColumn(V,i,&v));
      PetscCall(VecGetArray(v,&uu));
      PetscCall(MatDenseGetColumn(Ux[idx],i,&array));
      PetscCall(PetscArraycpy(array,uu,eps->nloc));
      PetscCall(MatDenseRestoreColumn(Ux[idx],&array));
      PetscCall(VecRestoreArray(v,&uu));
      PetscCall(BVRestoreColumn(V,i,&v));
    }

    /* Eigenpair approximation */
    PetscCall(BVGetColumn(V,0,&v));
    PetscCall(MatMult(S,v,z));
    PetscCall(VecDot(z,v,pM));
    PetscCall(BVRestoreColumn(V,0,&v));
    if (rk>1) {
      PetscCall(BVGetColumn(V,1,&w));
      PetscCall(VecDot(z,w,pM+1));
      PetscCall(MatMult(S,w,z));
      PetscCall(VecDot(z,w,pM+3));
      PetscCall(BVGetColumn(V,0,&v));
      PetscCall(VecDot(z,v,pM+2));
      PetscCall(BVRestoreColumn(V,0,&v));
      PetscCall(BVRestoreColumn(V,1,&w));
      PetscCall(EV2x2(pM,2,eigr,eigi,vec));
      PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,2,2,vec,&X));
      PetscCall(BVSetActiveColumns(V,0,rk));
      PetscCall(BVMultInPlace(V,X,0,rk));
      PetscCall(MatDestroy(&X));
#if !defined(PETSC_USE_COMPLEX)
      norm = eigr[0]*eigr[0]+eigi[0]*eigi[0];
      er = eigr[0]/norm; ei = -eigi[0]/norm;
#else
      er =1.0/eigr[0]; ei = 0.0;
#endif
    } else {
      eigr[0] = pM[0]; eigi[0] = 0.0;
      er = 1.0/eigr[0]; ei = 0.0;
    }
    PetscCall(BVGetColumn(V,0,&v));
    if (eigi[0]!=0.0) PetscCall(BVGetColumn(V,1,&w));
    else w = NULL;
    eps->eigr[eps->nconv] = eigr[0]; eps->eigi[eps->nconv] = eigi[0];
    PetscCall(EPSComputeResidualNorm_Private(eps,PETSC_FALSE,er,ei,v,w,eps->work,&norm));
    PetscCall(BVRestoreColumn(V,0,&v));
    if (w) PetscCall(BVRestoreColumn(V,1,&w));
    PetscCall((*eps->converged)(eps,er,ei,norm,&eps->errest[eps->nconv],eps->convergedctx));
    k = 0;
    if (eps->errest[eps->nconv]<eps->tol) {
      k++;
      if (rk==2) {
#if !defined (PETSC_USE_COMPLEX)
        eps->eigr[eps->nconv+k] = eigr[0]; eps->eigi[eps->nconv+k] = -eigi[0];
#else
        eps->eigr[eps->nconv+k] = PetscConj(eps->eigr[eps->nconv]);
#endif
        k++;
      }
      /* Store converged eigenpairs and vectors for deflation */
      for (i=0;i<k;i++) {
        PetscCall(BVGetColumn(V,i,&v));
        PetscCall(BVInsertVec(eps->V,eps->nconv+i,v));
        PetscCall(BVRestoreColumn(V,i,&v));
      }
      eps->nconv += k;
      PetscCall(BVSetActiveColumns(eps->V,eps->nconv-rk,eps->nconv));
      PetscCall(BVOrthogonalize(eps->V,NULL));
      PetscCall(DSSetDimensions(eps->ds,eps->nconv,0,0));
      PetscCall(DSGetMat(eps->ds,DS_MAT_A,&W));
      PetscCall(BVMatProject(eps->V,matctx->S,eps->V,W));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_A,&W));
      if (eps->nconv<eps->nev) {
        idx = 0;
        PetscCall(BVSetRandomColumn(V,0));
        PetscCall(BVNormColumn(V,0,NORM_2,&norm));
        PetscCall(BVScaleColumn(V,0,1.0/norm));
        PetscCall(LyapIIBuildRHS(S,1,Ux[idx],V,eps->work));
      }
    } else {
      /* Prepare right-hand side */
      PetscCall(LyapIIBuildRHS(S,rk,Ux[idx],NULL,eps->work));
    }
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx));
    PetscCall(EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,eps->nconv+1));
  }
  PetscCall(STRestoreOperator(eps->st,&matctx->S));
  PetscCall(MatDestroy(&S));
  PetscCall(MatDestroy(&Ux[0]));
  PetscCall(MatDestroy(&Ux[1]));
  PetscCall(MatDestroy(&Op));
  PetscCall(VecDestroy(&v0));
  PetscCall(BVDestroy(&V));
  PetscCall(EPSDestroy(&epsrr));
  PetscCall(PetscFree(s));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_LyapII(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;
  PetscInt       k,array[2]={PETSC_DETERMINE,PETSC_DETERMINE};
  PetscBool      flg;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Lyapunov Inverse Iteration Options");

    k = 2;
    PetscCall(PetscOptionsIntArray("-eps_lyapii_ranks","Ranks for Lyapunov equation (one or two comma-separated integers)","EPSLyapIISetRanks",array,&k,&flg));
    if (flg) PetscCall(EPSLyapIISetRanks(eps,array[0],array[1]));

  PetscOptionsHeadEnd();

  if (!ctx->lme) PetscCall(EPSLyapIIGetLME(eps,&ctx->lme));
  PetscCall(LMESetFromOptions(ctx->lme));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLyapIISetRanks_LyapII(EPS eps,PetscInt rkc,PetscInt rkl)
{
  EPS_LYAPII *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (rkc==PETSC_DETERMINE) {
    if (ctx->rkc != 10) eps->state = EPS_STATE_INITIAL;
    ctx->rkc = 10;
  } else if (rkc!=PETSC_CURRENT) {
    PetscCheck(rkc>1,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The compressed rank %" PetscInt_FMT " must be larger than 1",rkc);
    if (ctx->rkc != rkc) eps->state = EPS_STATE_INITIAL;
    ctx->rkc = rkc;
  }
  if (rkl==PETSC_DETERMINE) {
    if (ctx->rkl != 3*rkc) eps->state = EPS_STATE_INITIAL;
    ctx->rkl = 3*rkc;
  } else if (rkl!=PETSC_CURRENT) {
    PetscCheck(rkl>=rkc,PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The Lyapunov rank %" PetscInt_FMT " cannot be smaller than the compressed rank %" PetscInt_FMT,rkl,rkc);
    if (ctx->rkl != rkl) eps->state = EPS_STATE_INITIAL;
    ctx->rkl = rkl;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLyapIISetRanks - Set the ranks used in the solution of the Lyapunov equation.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
.  rkc - the compressed rank
-  rkl - the Lyapunov rank

   Options Database Key:
.  -eps_lyapii_ranks <rkc,rkl> - Sets the rank parameters

   Notes:
   PETSC_CURRENT can be used to preserve the current value of any of the
   arguments, and PETSC_DETERMINE to set them to a default value.

   Lyapunov inverse iteration needs to solve a large-scale Lyapunov equation
   at each iteration of the eigensolver. For this, an iterative solver (LME)
   is used, which requires to prescribe the rank of the solution matrix X. This
   is the meaning of parameter rkl. Later, this matrix is compressed into
   another matrix of rank rkc. If not provided, rkl is a small multiple of rkc.

   Level: intermediate

.seealso: EPSLyapIIGetRanks()
@*/
PetscErrorCode EPSLyapIISetRanks(EPS eps,PetscInt rkc,PetscInt rkl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,rkc,2);
  PetscValidLogicalCollectiveInt(eps,rkl,3);
  PetscTryMethod(eps,"EPSLyapIISetRanks_C",(EPS,PetscInt,PetscInt),(eps,rkc,rkl));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLyapIIGetRanks_LyapII(EPS eps,PetscInt *rkc,PetscInt *rkl)
{
  EPS_LYAPII *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (rkc) *rkc = ctx->rkc;
  if (rkl) *rkl = ctx->rkl;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLyapIIGetRanks - Return the rank values used for the Lyapunov step.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  rkc - the compressed rank
-  rkl - the Lyapunov rank

   Level: intermediate

.seealso: EPSLyapIISetRanks()
@*/
PetscErrorCode EPSLyapIIGetRanks(EPS eps,PetscInt *rkc,PetscInt *rkl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscUseMethod(eps,"EPSLyapIIGetRanks_C",(EPS,PetscInt*,PetscInt*),(eps,rkc,rkl));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLyapIISetLME_LyapII(EPS eps,LME lme)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)lme));
  PetscCall(LMEDestroy(&ctx->lme));
  ctx->lme = lme;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLyapIISetLME - Associate a linear matrix equation solver object (LME) to the
   eigenvalue solver.

   Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  lme - the linear matrix equation solver object

   Level: advanced

.seealso: EPSLyapIIGetLME()
@*/
PetscErrorCode EPSLyapIISetLME(EPS eps,LME lme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(lme,LME_CLASSID,2);
  PetscCheckSameComm(eps,1,lme,2);
  PetscTryMethod(eps,"EPSLyapIISetLME_C",(EPS,LME),(eps,lme));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLyapIIGetLME_LyapII(EPS eps,LME *lme)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (!ctx->lme) {
    PetscCall(LMECreate(PetscObjectComm((PetscObject)eps),&ctx->lme));
    PetscCall(LMESetOptionsPrefix(ctx->lme,((PetscObject)eps)->prefix));
    PetscCall(LMEAppendOptionsPrefix(ctx->lme,"eps_lyapii_"));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)ctx->lme,(PetscObject)eps,1));
  }
  *lme = ctx->lme;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLyapIIGetLME - Retrieve the linear matrix equation solver object (LME)
   associated with the eigenvalue solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lme - the linear matrix equation solver object

   Level: advanced

.seealso: EPSLyapIISetLME()
@*/
PetscErrorCode EPSLyapIIGetLME(EPS eps,LME *lme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(lme,2);
  PetscUseMethod(eps,"EPSLyapIIGetLME_C",(EPS,LME*),(eps,lme));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_LyapII(EPS eps,PetscViewer viewer)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  ranks: for Lyapunov solver=%" PetscInt_FMT ", after compression=%" PetscInt_FMT "\n",ctx->rkl,ctx->rkc));
    if (!ctx->lme) PetscCall(EPSLyapIIGetLME(eps,&ctx->lme));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(LMEView(ctx->lme,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_LyapII(EPS eps)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (!ctx->lme) PetscCall(LMEReset(ctx->lme));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_LyapII(EPS eps)
{
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  PetscCall(LMEDestroy(&ctx->lme));
  PetscCall(DSDestroy(&ctx->ds));
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetLME_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetLME_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetRanks_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetRanks_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetDefaultST_LyapII(EPS eps)
{
  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) PetscCall(STSetType(eps->st,STSINVERT));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_LyapII(EPS eps)
{
  EPS_LYAPII     *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_LyapII;
  eps->ops->setup          = EPSSetUp_LyapII;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_LyapII;
  eps->ops->reset          = EPSReset_LyapII;
  eps->ops->destroy        = EPSDestroy_LyapII;
  eps->ops->view           = EPSView_LyapII;
  eps->ops->setdefaultst   = EPSSetDefaultST_LyapII;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Schur;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetLME_C",EPSLyapIISetLME_LyapII));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetLME_C",EPSLyapIIGetLME_LyapII));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetRanks_C",EPSLyapIISetRanks_LyapII));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetRanks_C",EPSLyapIIGetRanks_LyapII));
  PetscFunctionReturn(PETSC_SUCCESS);
}
