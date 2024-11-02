/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "jd"

   Method: Jacobi-Davidson

   Algorithm:

       Jacobi-Davidson for polynomial eigenvalue problems.

   References:

       [1] C. Campos and J.E. Roman, "A polynomial Jacobi-Davidson solver
           with support for non-monomial bases and deflation", BIT Numer.
           Math. 60:295-318, 2020.

       [2] G.L.G. Sleijpen et al., "Jacobi-Davidson type methods for
           generalized eigenproblems and polynomial eigenproblems", BIT
           36(3):595-633, 1996.

       [3] Feng-Nan Hwang, Zih-Hao Wei, Tsung-Ming Huang, Weichung Wang,
           "A Parallel Additive Schwarz Preconditioned Jacobi-Davidson
           Algorithm for Polynomial Eigenvalue Problems in Quantum Dot
           Simulation", J. Comput. Phys. 229(8):2932-2947, 2010.
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-slice-qep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"A polynomial {Jacobi-Davidson} solver with support for non-monomial bases and deflation\",\n"
  "   journal = \"{BIT} Numer. Math.\",\n"
  "   volume = \"60\",\n"
  "   pages = \"295--318\",\n"
  "   year = \"2020,\"\n"
  "   doi = \"https://doi.org/10.1007/s10543-019-00778-z\"\n"
  "}\n";

typedef struct {
  PetscReal   keep;          /* restart parameter */
  PetscReal   fix;           /* fix parameter */
  PetscBool   reusepc;       /* flag indicating whether pc is rebuilt or not */
  BV          V;             /* work basis vectors to store the search space */
  BV          W;             /* work basis vectors to store the test space */
  BV          *TV;           /* work basis vectors to store T*V (each TV[i] is the coefficient for \lambda^i of T*V for the extended T) */
  BV          *AX;           /* work basis vectors to store A_i*X for locked eigenvectors */
  BV          N[2];          /* auxiliary work BVs */
  BV          X;             /* locked eigenvectors */
  PetscScalar *T;            /* matrix of the invariant pair */
  PetscScalar *Tj;           /* matrix containing the powers of the invariant pair matrix */
  PetscScalar *XpX;          /* X^H*X */
  PetscInt    ld;            /* leading dimension for Tj and XpX */
  PC          pcshell;       /* preconditioner including basic precond+projector */
  Mat         Pshell;        /* auxiliary shell matrix */
  PetscInt    nlock;         /* number of locked vectors in the invariant pair */
  Vec         vtempl;        /* reference nested vector */
  PetscInt    midx;          /* minimality index */
  PetscInt    mmidx;         /* maximum allowed minimality index */
  PEPJDProjection proj;      /* projection type (orthogonal, harmonic) */
} PEP_JD;

typedef struct {
  PEP         pep;
  PC          pc;            /* basic preconditioner */
  Vec         Bp[2];         /* preconditioned residual of derivative polynomial, B\p */
  Vec         u[2];          /* Ritz vector */
  PetscScalar gamma[2];      /* precomputed scalar u'*B\p */
  PetscScalar theta;
  PetscScalar *M;
  PetscScalar *ps;
  PetscInt    ld;
  Vec         *work;
  Mat         PPr;
  BV          X;
  PetscInt    n;
} PEP_JD_PCSHELL;

typedef struct {
  Mat         Pr,Pi;         /* matrix polynomial evaluated at theta */
  PEP         pep;
  Vec         *work;
  PetscScalar theta[2];
} PEP_JD_MATSHELL;

/*
   Duplicate and resize auxiliary basis
*/
static PetscErrorCode PEPJDDuplicateBasis(PEP pep,BV *basis)
{
  PEP_JD             *pjd = (PEP_JD*)pep->data;
  PetscInt           nloc,m;
  BVType             type;
  BVOrthogType       otype;
  BVOrthogRefineType oref;
  PetscReal          oeta;
  BVOrthogBlockType  oblock;
  VecType            vtype;

  PetscFunctionBegin;
  if (pjd->ld>1) {
    PetscCall(BVCreate(PetscObjectComm((PetscObject)pep),basis));
    PetscCall(BVGetSizes(pep->V,&nloc,NULL,&m));
    nloc += pjd->ld-1;
    PetscCall(BVSetSizes(*basis,nloc,PETSC_DECIDE,m));
    PetscCall(BVGetType(pep->V,&type));
    PetscCall(BVSetType(*basis,type));
    PetscCall(BVGetVecType(pep->V,&vtype));
    PetscCall(BVSetVecType(*basis,vtype));
    PetscCall(BVGetOrthogonalization(pep->V,&otype,&oref,&oeta,&oblock));
    PetscCall(BVSetOrthogonalization(*basis,otype,oref,oeta,oblock));
    PetscCall(PetscObjectStateIncrease((PetscObject)*basis));
  } else PetscCall(BVDuplicate(pep->V,basis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetUp_JD(PEP pep)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBool      isprecond,flg;
  PetscRandom    rand;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd));
  if (pep->max_it==PETSC_DETERMINE) pep->max_it = PetscMax(100,2*pep->n/pep->ncv);
  if (!pep->which) pep->which = PEP_TARGET_MAGNITUDE;
  PetscCheck(pep->which==PEP_TARGET_MAGNITUDE || pep->which==PEP_TARGET_REAL || pep->which==PEP_TARGET_IMAGINARY,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"The JD solver supports only target which, see PEPSetWhichEigenpairs()");

  PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STPRECOND,&isprecond));
  PetscCheck(isprecond,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"The JD solver only works with PRECOND spectral transformation");

  PetscCall(STGetTransform(pep->st,&flg));
  PetscCheck(!flg,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"The JD solver requires the ST transform flag unset, see STSetTransform()");
  PEPCheckIgnored(pep,PEP_FEATURE_EXTRACT);

  if (!pjd->mmidx) pjd->mmidx = pep->nmat-1;
  pjd->mmidx = PetscMin(pjd->mmidx,pep->nmat-1);
  if (!pjd->keep) pjd->keep = 0.5;
  PetscCall(PEPBasisCoefficients(pep,pep->pbc));
  PetscCall(PEPAllocateSolution(pep,0));
  PetscCall(BVGetRandomContext(pep->V,&rand));  /* make sure the random context is available when duplicating */
  PetscCall(PEPSetWorkVecs(pep,5));
  pjd->ld = pep->nev;
#if !defined (PETSC_USE_COMPLEX)
  pjd->ld++;
#endif
  PetscCall(PetscMalloc2(pep->nmat,&pjd->TV,pep->nmat,&pjd->AX));
  for (i=0;i<pep->nmat;i++) PetscCall(PEPJDDuplicateBasis(pep,pjd->TV+i));
  if (pjd->ld>1) {
    PetscCall(PEPJDDuplicateBasis(pep,&pjd->V));
    PetscCall(BVSetFromOptions(pjd->V));
    for (i=0;i<pep->nmat;i++) PetscCall(BVDuplicateResize(pep->V,pjd->ld-1,pjd->AX+i));
    PetscCall(BVDuplicateResize(pep->V,pjd->ld-1,pjd->N));
    PetscCall(BVDuplicateResize(pep->V,pjd->ld-1,pjd->N+1));
    pjd->X = pep->V;
    PetscCall(PetscCalloc3(pjd->ld*pjd->ld,&pjd->XpX,pep->ncv*pep->ncv,&pjd->T,pjd->ld*pjd->ld*pep->nmat,&pjd->Tj));
  } else pjd->V = pep->V;
  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) PetscCall(PEPJDDuplicateBasis(pep,&pjd->W));
  else pjd->W = pjd->V;
  PetscCall(DSSetType(pep->ds,DSPEP));
  PetscCall(DSPEPSetDegree(pep->ds,pep->nmat-1));
  if (pep->basis!=PEP_BASIS_MONOMIAL) PetscCall(DSPEPSetCoefficients(pep->ds,pep->pbc));
  PetscCall(DSAllocate(pep->ds,pep->ncv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Updates columns (low to (high-1)) of TV[i]
*/
static PetscErrorCode PEPJDUpdateTV(PEP pep,PetscInt low,PetscInt high,Vec *w)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       pp,col,i,nloc,nconv;
  Vec            v1,v2,t1,t2;
  PetscScalar    *array1,*array2,*x2,*xx,*N,*Np,*y2=NULL,zero=0.0,sone=1.0,*pT,fact,*psc;
  PetscReal      *cg,*ca,*cb;
  PetscMPIInt    rk,np;
  PetscBLASInt   n_,ld_,one=1;
  Mat            T;
  BV             pbv;

  PetscFunctionBegin;
  ca = pep->pbc; cb = ca+pep->nmat; cg = cb + pep->nmat;
  nconv = pjd->nlock;
  PetscCall(PetscMalloc5(nconv,&x2,nconv,&xx,nconv*nconv,&pT,nconv*nconv,&N,nconv*nconv,&Np));
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
  PetscCall(BVGetSizes(pep->V,&nloc,NULL,NULL));
  t1 = w[0];
  t2 = w[1];
  PetscCall(PetscBLASIntCast(pjd->nlock,&n_));
  PetscCall(PetscBLASIntCast(pjd->ld,&ld_));
  if (nconv) {
    for (i=0;i<nconv;i++) PetscCall(PetscArraycpy(pT+i*nconv,pjd->T+i*pep->ncv,nconv));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nconv,nconv,pT,&T));
  }
  for (col=low;col<high;col++) {
    PetscCall(BVGetColumn(pjd->V,col,&v1));
    PetscCall(VecGetArray(v1,&array1));
    if (nconv>0) {
      for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    }
    PetscCall(VecPlaceArray(t1,array1));
    if (nconv) {
      PetscCall(BVSetActiveColumns(pjd->N[0],0,nconv));
      PetscCall(BVSetActiveColumns(pjd->N[1],0,nconv));
      PetscCall(BVDotVec(pjd->X,t1,xx));
    }
    for (pp=pep->nmat-1;pp>=0;pp--) {
      PetscCall(BVGetColumn(pjd->TV[pp],col,&v2));
      PetscCall(VecGetArray(v2,&array2));
      PetscCall(VecPlaceArray(t2,array2));
      PetscCall(MatMult(pep->A[pp],t1,t2));
      if (nconv) {
        if (pp<pep->nmat-3) {
          PetscCall(BVMult(pjd->N[0],1.0,-cg[pp+2],pjd->AX[pp+1],NULL));
          PetscCall(MatShift(T,-cb[pp+1]));
          PetscCall(BVMult(pjd->N[0],1.0/ca[pp],1.0/ca[pp],pjd->N[1],T));
          pbv = pjd->N[0]; pjd->N[0] = pjd->N[1]; pjd->N[1] = pbv;
          PetscCall(BVMultVec(pjd->N[1],1.0,1.0,t2,x2));
          PetscCall(MatShift(T,cb[pp+1]));
        } else if (pp==pep->nmat-3) {
          PetscCall(BVCopy(pjd->AX[pp+2],pjd->N[0]));
          PetscCall(BVScale(pjd->N[0],1/ca[pp+1]));
          PetscCall(BVCopy(pjd->AX[pp+1],pjd->N[1]));
          PetscCall(MatShift(T,-cb[pp+1]));
          PetscCall(BVMult(pjd->N[1],1.0/ca[pp],1.0/ca[pp],pjd->N[0],T));
          PetscCall(BVMultVec(pjd->N[1],1.0,1.0,t2,x2));
          PetscCall(MatShift(T,cb[pp+1]));
        } else if (pp==pep->nmat-2) PetscCall(BVMultVec(pjd->AX[pp+1],1.0/ca[pp],1.0,t2,x2));
        if (pp<pjd->midx) {
          y2 = array2+nloc;
          PetscCallBLAS("BLASgemv",BLASgemv_("C",&n_,&n_,&sone,pjd->Tj+pjd->ld*pjd->ld*pp,&ld_,xx,&one,&zero,y2,&one));
          if (pp<pjd->midx-2) {
            fact = -cg[pp+2];
            PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,pjd->Tj+(pp+1)*pjd->ld*pjd->ld,&ld_,pjd->XpX,&ld_,&fact,Np,&n_));
            fact = 1/ca[pp];
            PetscCall(MatShift(T,-cb[pp+1]));
            PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&fact,N,&n_,pT,&n_,&fact,Np,&n_));
            PetscCall(MatShift(T,cb[pp+1]));
            psc = Np; Np = N; N = psc;
            PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,N,&n_,x2,&one,&sone,y2,&one));
          } else if (pp==pjd->midx-2) {
            fact = 1/ca[pp];
            PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&fact,pjd->Tj+(pp+1)*pjd->ld*pjd->ld,&ld_,pjd->XpX,&ld_,&zero,N,&n_));
            PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,N,&n_,x2,&one,&sone,y2,&one));
          } else if (pp==pjd->midx-1) PetscCall(PetscArrayzero(Np,nconv*nconv));
        }
        for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
      }
      PetscCall(VecResetArray(t2));
      PetscCall(VecRestoreArray(v2,&array2));
      PetscCall(BVRestoreColumn(pjd->TV[pp],col,&v2));
    }
    PetscCall(VecResetArray(t1));
    PetscCall(VecRestoreArray(v1,&array1));
    PetscCall(BVRestoreColumn(pjd->V,col,&v1));
  }
  if (nconv) PetscCall(MatDestroy(&T));
  PetscCall(PetscFree5(x2,xx,pT,N,Np));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   RRQR of X. Xin*P=Xou*R. Rank of R is rk
*/
static PetscErrorCode PEPJDOrthogonalize(PetscInt row,PetscInt col,PetscScalar *X,PetscInt ldx,PetscInt *rk,PetscInt *P,PetscScalar *R,PetscInt ldr)
{
  PetscInt       i,j,n,r;
  PetscBLASInt   row_,col_,ldx_,*p,lwork,info,n_;
  PetscScalar    *tau,*work;
  PetscReal      tol,*rwork;

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(row,&row_));
  PetscCall(PetscBLASIntCast(col,&col_));
  PetscCall(PetscBLASIntCast(ldx,&ldx_));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  n = PetscMin(row,col);
  PetscCall(PetscBLASIntCast(n,&n_));
  lwork = 3*col_+1;
  PetscCall(PetscMalloc4(col,&p,n,&tau,lwork,&work,2*col,&rwork));
  for (i=1;i<col;i++) p[i] = 0;
  p[0] = 1;

  /* rank revealing QR */
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKgeqp3",LAPACKgeqp3_(&row_,&col_,X,&ldx_,p,tau,work,&lwork,rwork,&info));
#else
  PetscCallBLAS("LAPACKgeqp3",LAPACKgeqp3_(&row_,&col_,X,&ldx_,p,tau,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("geqp3",info);
  if (P) for (i=0;i<col;i++) P[i] = p[i]-1;

  /* rank computation */
  tol = PetscMax(row,col)*PETSC_MACHINE_EPSILON*PetscAbsScalar(X[0]);
  r = 1;
  for (i=1;i<n;i++) {
    if (PetscAbsScalar(X[i+ldx*i])>tol) r++;
    else break;
  }
  if (rk) *rk=r;

  /* copy upper triangular matrix if requested */
  if (R) {
     for (i=0;i<r;i++) {
       PetscCall(PetscArrayzero(R+i*ldr,r));
       for (j=0;j<=i;j++) R[i*ldr+j] = X[i*ldx+j];
     }
  }
  PetscCallBLAS("LAPACKorgqr",LAPACKorgqr_(&row_,&n_,&n_,X,&ldx_,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("orgqr",info);
  PetscCall(PetscFPTrapPop());
  PetscCall(PetscFree4(p,tau,work,rwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Application of extended preconditioner
*/
static PetscErrorCode PEPJDExtendedPCApply(PC pc,Vec x,Vec y)
{
  PetscInt          i,j,nloc,n,ld=0;
  PetscMPIInt       np;
  Vec               tx,ty;
  PEP_JD_PCSHELL    *ctx;
  const PetscScalar *array1;
  PetscScalar       *x2=NULL,*t=NULL,*ps=NULL,*array2,zero=0.0,sone=1.0;
  PetscBLASInt      one=1,ld_,n_,ncv_;
  PEP_JD            *pjd=NULL;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pc),&np));
  PetscCall(PCShellGetContext(pc,&ctx));
  n  = ctx->n;
  if (n) {
    pjd = (PEP_JD*)ctx->pep->data;
    ps = ctx->ps;
    ld = pjd->ld;
    PetscCall(PetscMalloc2(n,&x2,n,&t));
    PetscCall(VecGetLocalSize(ctx->work[0],&nloc));
    PetscCall(VecGetArrayRead(x,&array1));
    for (i=0;i<n;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    PetscCall(VecRestoreArrayRead(x,&array1));
  }

  /* y = B\x apply PC */
  tx = ctx->work[0];
  ty = ctx->work[1];
  PetscCall(VecGetArrayRead(x,&array1));
  PetscCall(VecPlaceArray(tx,array1));
  PetscCall(VecGetArray(y,&array2));
  PetscCall(VecPlaceArray(ty,array2));
  PetscCall(PCApply(ctx->pc,tx,ty));
  if (n) {
    PetscCall(PetscBLASIntCast(ld,&ld_));
    PetscCall(PetscBLASIntCast(n,&n_));
    for (i=0;i<n;i++) {
      t[i] = 0.0;
      for (j=0;j<n;j++) t[i] += ctx->M[i+j*ld]*x2[j];
    }
    if (pjd->midx==1) {
      PetscCall(PetscBLASIntCast(ctx->pep->ncv,&ncv_));
      for (i=0;i<n;i++) pjd->T[i*(1+ctx->pep->ncv)] -= ctx->theta;
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,pjd->T,&ncv_,t,&one,&zero,x2,&one));
      for (i=0;i<n;i++) pjd->T[i*(1+ctx->pep->ncv)] += ctx->theta;
      for (i=0;i<n;i++) array2[nloc+i] = x2[i];
      for (i=0;i<n;i++) x2[i] = -t[i];
    } else {
      for (i=0;i<n;i++) array2[nloc+i] = t[i];
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,ps,&ld_,t,&one,&zero,x2,&one));
    }
    for (i=0;i<n;i++) array2[nloc+i] /= PetscSqrtReal(np);
    PetscCall(BVSetActiveColumns(pjd->X,0,n));
    PetscCall(BVMultVec(pjd->X,-1.0,1.0,ty,x2));
    PetscCall(PetscFree2(x2,t));
  }
  PetscCall(VecResetArray(tx));
  PetscCall(VecResetArray(ty));
  PetscCall(VecRestoreArrayRead(x,&array1));
  PetscCall(VecRestoreArray(y,&array2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Application of shell preconditioner:
      y = B\x - eta*B\p,  with eta = (u'*B\x)/(u'*B\p)
*/
static PetscErrorCode PCShellApply_PEPJD(PC pc,Vec x,Vec y)
{
  PetscScalar    rr,eta;
  PEP_JD_PCSHELL *ctx;
  PetscInt       sz;
  const Vec      *xs,*ys;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    rx,xr,xx;
#endif

  PetscFunctionBegin;
  PetscCall(PCShellGetContext(pc,&ctx));
  PetscCall(VecCompGetSubVecs(x,&sz,&xs));
  PetscCall(VecCompGetSubVecs(y,NULL,&ys));
  /* y = B\x apply extended PC */
  PetscCall(PEPJDExtendedPCApply(pc,xs[0],ys[0]));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) PetscCall(PEPJDExtendedPCApply(pc,xs[1],ys[1]));
#endif

  /* Compute eta = u'*y / u'*Bp */
  PetscCall(VecDot(ys[0],ctx->u[0],&rr));
  eta  = -rr*ctx->gamma[0];
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    PetscCall(VecDot(ys[0],ctx->u[1],&xr));
    PetscCall(VecDot(ys[1],ctx->u[0],&rx));
    PetscCall(VecDot(ys[1],ctx->u[1],&xx));
    eta += -ctx->gamma[0]*xx-ctx->gamma[1]*(-xr+rx);
  }
#endif
  eta /= ctx->gamma[0]*ctx->gamma[0]+ctx->gamma[1]*ctx->gamma[1];

  /* y = y - eta*Bp */
  PetscCall(VecAXPY(ys[0],eta,ctx->Bp[0]));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    PetscCall(VecAXPY(ys[1],eta,ctx->Bp[1]));
    eta = -ctx->gamma[1]*(rr+xx)+ctx->gamma[0]*(-xr+rx);
    eta /= ctx->gamma[0]*ctx->gamma[0]+ctx->gamma[1]*ctx->gamma[1];
    PetscCall(VecAXPY(ys[0],eta,ctx->Bp[1]));
    PetscCall(VecAXPY(ys[1],-eta,ctx->Bp[0]));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDCopyToExtendedVec(PEP pep,Vec v,PetscScalar *a,PetscInt na,PetscInt off,Vec vex,PetscBool back)
{
  PetscMPIInt    np,rk,count;
  PetscScalar    *array1,*array2;
  PetscInt       nloc;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
  PetscCall(BVGetSizes(pep->V,&nloc,NULL,NULL));
  if (v) {
    PetscCall(VecGetArray(v,&array1));
    PetscCall(VecGetArray(vex,&array2));
    if (back) PetscCall(PetscArraycpy(array1,array2,nloc));
    else PetscCall(PetscArraycpy(array2,array1,nloc));
    PetscCall(VecRestoreArray(v,&array1));
    PetscCall(VecRestoreArray(vex,&array2));
  }
  if (a) {
    PetscCall(VecGetArray(vex,&array2));
    if (back) {
      PetscCall(PetscArraycpy(a,array2+nloc+off,na));
      PetscCall(PetscMPIIntCast(na,&count));
      PetscCallMPI(MPI_Bcast(a,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
    } else {
      PetscCall(PetscArraycpy(array2+nloc+off,a,na));
      PetscCall(PetscMPIIntCast(na,&count));
      PetscCallMPI(MPI_Bcast(array2+nloc+off,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
    }
    PetscCall(VecRestoreArray(vex,&array2));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Computes Phi^hat(lambda) times a vector or its derivative (depends on beval)
     if no vector is provided returns a matrix
 */
static PetscErrorCode PEPJDEvaluateHatBasis(PEP pep,PetscInt n,PetscScalar *H,PetscInt ldh,PetscScalar *beval,PetscScalar *t,PetscInt idx,PetscScalar *qpp,PetscScalar *qp,PetscScalar *q)
{
  PetscInt       j,i;
  PetscBLASInt   n_,ldh_,one=1;
  PetscReal      *a,*b,*g;
  PetscScalar    sone=1.0,zero=0.0;

  PetscFunctionBegin;
  a = pep->pbc; b=a+pep->nmat; g=b+pep->nmat;
  PetscCall(PetscBLASIntCast(n,&n_));
  PetscCall(PetscBLASIntCast(ldh,&ldh_));
  if (idx<1) PetscCall(PetscArrayzero(q,t?n:n*n));
  else if (idx==1) {
    if (t) {for (j=0;j<n;j++) q[j] = t[j]*beval[idx-1]/a[0];}
    else {
      PetscCall(PetscArrayzero(q,n*n));
      for (j=0;j<n;j++) q[(j+1)*n] = beval[idx-1]/a[0];
    }
  } else {
    if (t) {
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,H,&ldh_,qp,&one,&zero,q,&one));
      for (j=0;j<n;j++) {
        q[j] += beval[idx-1]*t[j]-b[idx-1]*qp[j]-g[idx-1]*qpp[j];
        q[j] /= a[idx-1];
      }
    } else {
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,H,&ldh_,qp,&n_,&zero,q,&n_));
      for (j=0;j<n;j++) {
        q[j+n*j] += beval[idx-1];
        for (i=0;i<n;i++) {
          q[i+n*j] += -b[idx-1]*qp[j*n+i]-g[idx-1]*qpp[j*n+i];
          q[i+n*j] /= a[idx-1];
        }
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDComputeResidual(PEP pep,PetscBool derivative,PetscInt sz,Vec *u,PetscScalar *theta,Vec *p,Vec *work)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscMPIInt    rk,np,count;
  Vec            tu,tp,w;
  PetscScalar    *dval,*dvali,*array1,*array2,*x2=NULL,*y2,*qj=NULL,*tt=NULL,*xx=NULL,*xxi=NULL,sone=1.0;
  PetscInt       i,j,nconv,nloc;
  PetscBLASInt   n,ld,one=1;
#if !defined(PETSC_USE_COMPLEX)
  Vec            tui=NULL,tpi=NULL;
  PetscScalar    *x2i=NULL,*qji=NULL,*qq,*y2i,*arrayi1,*arrayi2;
#endif

  PetscFunctionBegin;
  nconv = pjd->nlock;
  if (!nconv) PetscCall(PetscMalloc1(2*sz*pep->nmat,&dval));
  else {
    PetscCall(PetscMalloc5(2*pep->nmat,&dval,2*nconv,&xx,nconv,&tt,sz*nconv,&x2,(sz==2?3:1)*nconv*pep->nmat,&qj));
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
    PetscCall(BVGetSizes(pep->V,&nloc,NULL,NULL));
    PetscCall(VecGetArray(u[0],&array1));
    for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]*PetscSqrtReal(np);
    PetscCall(VecRestoreArray(u[0],&array1));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      x2i = x2+nconv;
      PetscCall(VecGetArray(u[1],&arrayi1));
      for (i=0;i<nconv;i++) x2i[i] = arrayi1[nloc+i]*PetscSqrtReal(np);
      PetscCall(VecRestoreArray(u[1],&arrayi1));
    }
#endif
  }
  dvali = dval+pep->nmat;
  tu = work[0];
  tp = work[1];
  w  = work[2];
  PetscCall(VecGetArray(u[0],&array1));
  PetscCall(VecPlaceArray(tu,array1));
  PetscCall(VecGetArray(p[0],&array2));
  PetscCall(VecPlaceArray(tp,array2));
  PetscCall(VecSet(tp,0.0));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    tui = work[3];
    tpi = work[4];
    PetscCall(VecGetArray(u[1],&arrayi1));
    PetscCall(VecPlaceArray(tui,arrayi1));
    PetscCall(VecGetArray(p[1],&arrayi2));
    PetscCall(VecPlaceArray(tpi,arrayi2));
    PetscCall(VecSet(tpi,0.0));
  }
#endif
  if (derivative) PetscCall(PEPEvaluateBasisDerivative(pep,theta[0],theta[1],dval,dvali));
  else PetscCall(PEPEvaluateBasis(pep,theta[0],theta[1],dval,dvali));
  for (i=derivative?1:0;i<pep->nmat;i++) {
    PetscCall(MatMult(pep->A[i],tu,w));
    PetscCall(VecAXPY(tp,dval[i],w));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      PetscCall(VecAXPY(tpi,dvali[i],w));
      PetscCall(MatMult(pep->A[i],tui,w));
      PetscCall(VecAXPY(tpi,dval[i],w));
      PetscCall(VecAXPY(tp,-dvali[i],w));
    }
#endif
  }
  if (nconv) {
    for (i=0;i<pep->nmat;i++) PetscCall(PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dval,x2,i,i>1?qj+(i-2)*nconv:NULL,i>0?qj+(i-1)*nconv:NULL,qj+i*nconv));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      qji = qj+nconv*pep->nmat;
      qq = qji+nconv*pep->nmat;
      for (i=0;i<pep->nmat;i++) PetscCall(PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dvali,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv));
      for (i=0;i<nconv*pep->nmat;i++) qj[i] -= qji[i];
      for (i=0;i<pep->nmat;i++) {
        PetscCall(PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dval,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv));
        PetscCall(PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dvali,x2,i,i>1?qq+(i-2)*nconv:NULL,i>0?qq+(i-1)*nconv:NULL,qq+i*nconv));
      }
      for (i=0;i<nconv*pep->nmat;i++) qji[i] += qq[i];
      for (i=derivative?2:1;i<pep->nmat;i++) PetscCall(BVMultVec(pjd->AX[i],1.0,1.0,tpi,qji+i*nconv));
    }
#endif
    for (i=derivative?2:1;i<pep->nmat;i++) PetscCall(BVMultVec(pjd->AX[i],1.0,1.0,tp,qj+i*nconv));

    /* extended vector part */
    PetscCall(BVSetActiveColumns(pjd->X,0,nconv));
    PetscCall(BVDotVec(pjd->X,tu,xx));
    xxi = xx+nconv;
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) PetscCall(BVDotVec(pjd->X,tui,xxi));
#endif
    if (sz==1) PetscCall(PetscArrayzero(xxi,nconv));
    if (rk==np-1) {
      PetscCall(PetscBLASIntCast(nconv,&n));
      PetscCall(PetscBLASIntCast(pjd->ld,&ld));
      y2  = array2+nloc;
      PetscCall(PetscArrayzero(y2,nconv));
      for (j=derivative?1:0;j<pjd->midx;j++) {
        for (i=0;i<nconv;i++) tt[i] = dval[j]*xx[i]-dvali[j]*xxi[i];
        PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qj+j*nconv,&one,&sone,tt,&one));
        PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2,&one));
      }
      for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) {
        y2i = arrayi2+nloc;
        PetscCall(PetscArrayzero(y2i,nconv));
        for (j=derivative?1:0;j<pjd->midx;j++) {
          for (i=0;i<nconv;i++) tt[i] = dval[j]*xxi[i]+dvali[j]*xx[i];
          PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qji+j*nconv,&one,&sone,tt,&one));
          PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2i,&one));
        }
        for (i=0;i<nconv;i++) arrayi2[nloc+i] /= PetscSqrtReal(np);
      }
#endif
    }
    PetscCall(PetscMPIIntCast(nconv,&count));
    PetscCallMPI(MPI_Bcast(array2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) PetscCallMPI(MPI_Bcast(arrayi2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
#endif
  }
  if (nconv) PetscCall(PetscFree5(dval,xx,tt,x2,qj));
  else PetscCall(PetscFree(dval));
  PetscCall(VecResetArray(tu));
  PetscCall(VecRestoreArray(u[0],&array1));
  PetscCall(VecResetArray(tp));
  PetscCall(VecRestoreArray(p[0],&array2));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    PetscCall(VecResetArray(tui));
    PetscCall(VecRestoreArray(u[1],&arrayi1));
    PetscCall(VecResetArray(tpi));
    PetscCall(VecRestoreArray(p[1],&arrayi2));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDProcessInitialSpace(PEP pep,Vec *w)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscScalar    *tt,target[2];
  Vec            vg,wg;
  PetscInt       i;
  PetscReal      norm;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(pjd->ld-1,&tt));
  PetscCheck(pep->nini==0,PETSC_COMM_SELF,PETSC_ERR_SUP,"Support for initial vectors not implemented yet");
  PetscCall(BVSetRandomColumn(pjd->V,0));
  for (i=0;i<pjd->ld-1;i++) tt[i] = 0.0;
  PetscCall(BVGetColumn(pjd->V,0,&vg));
  PetscCall(PEPJDCopyToExtendedVec(pep,NULL,tt,pjd->ld-1,0,vg,PETSC_FALSE));
  PetscCall(BVRestoreColumn(pjd->V,0,&vg));
  PetscCall(BVNormColumn(pjd->V,0,NORM_2,&norm));
  PetscCall(BVScaleColumn(pjd->V,0,1.0/norm));
  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
    PetscCall(BVGetColumn(pjd->V,0,&vg));
    PetscCall(BVGetColumn(pjd->W,0,&wg));
    PetscCall(VecSet(wg,0.0));
    target[0] = pep->target; target[1] = 0.0;
    PetscCall(PEPJDComputeResidual(pep,PETSC_TRUE,1,&vg,target,&wg,w));
    PetscCall(BVRestoreColumn(pjd->W,0,&wg));
    PetscCall(BVRestoreColumn(pjd->V,0,&vg));
    PetscCall(BVNormColumn(pjd->W,0,NORM_2,&norm));
    PetscCall(BVScaleColumn(pjd->W,0,1.0/norm));
  }
  PetscCall(PetscFree(tt));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_PEPJD(Mat P,Vec x,Vec y)
{
  PEP_JD_MATSHELL *matctx;
  PEP_JD          *pjd;
  PetscInt        i,j,nconv,nloc,nmat,ldt,ncv,sz;
  Vec             tx,ty;
  const Vec       *xs,*ys;
  PetscScalar     *array1,*array2,*x2=NULL,*y2,*tt=NULL,*xx=NULL,*xxi,theta[2],sone=1.0,*qj,*val,*vali=NULL;
  PetscBLASInt    n,ld,one=1;
  PetscMPIInt     np;
#if !defined(PETSC_USE_COMPLEX)
  Vec             txi=NULL,tyi=NULL;
  PetscScalar     *x2i=NULL,*qji=NULL,*qq,*y2i,*arrayi1,*arrayi2;
#endif

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)P),&np));
  PetscCall(MatShellGetContext(P,&matctx));
  pjd   = (PEP_JD*)matctx->pep->data;
  nconv = pjd->nlock;
  nmat  = matctx->pep->nmat;
  ncv   = matctx->pep->ncv;
  ldt   = pjd->ld;
  PetscCall(VecCompGetSubVecs(x,&sz,&xs));
  PetscCall(VecCompGetSubVecs(y,NULL,&ys));
  theta[0] = matctx->theta[0];
  theta[1] = (sz==2)?matctx->theta[1]:0.0;
  if (nconv>0) {
    PetscCall(PetscMalloc5(nconv,&tt,sz*nconv,&x2,(sz==2?3:1)*nconv*nmat,&qj,2*nconv,&xx,2*nmat,&val));
    PetscCall(BVGetSizes(matctx->pep->V,&nloc,NULL,NULL));
    PetscCall(VecGetArray(xs[0],&array1));
    for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    PetscCall(VecRestoreArray(xs[0],&array1));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      x2i = x2+nconv;
      PetscCall(VecGetArray(xs[1],&arrayi1));
      for (i=0;i<nconv;i++) x2i[i] = arrayi1[nloc+i]* PetscSqrtReal(np);
      PetscCall(VecRestoreArray(xs[1],&arrayi1));
    }
#endif
    vali = val+nmat;
  }
  tx = matctx->work[0];
  ty = matctx->work[1];
  PetscCall(VecGetArray(xs[0],&array1));
  PetscCall(VecPlaceArray(tx,array1));
  PetscCall(VecGetArray(ys[0],&array2));
  PetscCall(VecPlaceArray(ty,array2));
  PetscCall(MatMult(matctx->Pr,tx,ty));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    txi = matctx->work[2];
    tyi = matctx->work[3];
    PetscCall(VecGetArray(xs[1],&arrayi1));
    PetscCall(VecPlaceArray(txi,arrayi1));
    PetscCall(VecGetArray(ys[1],&arrayi2));
    PetscCall(VecPlaceArray(tyi,arrayi2));
    PetscCall(MatMult(matctx->Pr,txi,tyi));
    if (theta[1]!=0.0) {
      PetscCall(MatMult(matctx->Pi,txi,matctx->work[4]));
      PetscCall(VecAXPY(ty,-1.0,matctx->work[4]));
      PetscCall(MatMult(matctx->Pi,tx,matctx->work[4]));
      PetscCall(VecAXPY(tyi,1.0,matctx->work[4]));
    }
  }
#endif
  if (nconv>0) {
    PetscCall(PEPEvaluateBasis(matctx->pep,theta[0],theta[1],val,vali));
    for (i=0;i<nmat;i++) PetscCall(PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,ncv,val,x2,i,i>1?qj+(i-2)*nconv:NULL,i>0?qj+(i-1)*nconv:NULL,qj+i*nconv));
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      qji = qj+nconv*nmat;
      qq = qji+nconv*nmat;
      for (i=0;i<nmat;i++) PetscCall(PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,vali,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv));
      for (i=0;i<nconv*nmat;i++) qj[i] -= qji[i];
      for (i=0;i<nmat;i++) {
        PetscCall(PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,val,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv));
        PetscCall(PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,vali,x2,i,i>1?qq+(i-2)*nconv:NULL,i>0?qq+(i-1)*nconv:NULL,qq+i*nconv));
      }
      for (i=0;i<nconv*nmat;i++) qji[i] += qq[i];
      for (i=1;i<matctx->pep->nmat;i++) PetscCall(BVMultVec(pjd->AX[i],1.0,1.0,tyi,qji+i*nconv));
    }
#endif
    for (i=1;i<nmat;i++) PetscCall(BVMultVec(pjd->AX[i],1.0,1.0,ty,qj+i*nconv));

    /* extended vector part */
    PetscCall(BVSetActiveColumns(pjd->X,0,nconv));
    PetscCall(BVDotVec(pjd->X,tx,xx));
    xxi = xx+nconv;
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) PetscCall(BVDotVec(pjd->X,txi,xxi));
#endif
    if (sz==1) PetscCall(PetscArrayzero(xxi,nconv));
    PetscCall(PetscBLASIntCast(pjd->nlock,&n));
    PetscCall(PetscBLASIntCast(ldt,&ld));
    y2 = array2+nloc;
    PetscCall(PetscArrayzero(y2,nconv));
    for (j=0;j<pjd->midx;j++) {
      for (i=0;i<nconv;i++) tt[i] = val[j]*xx[i]-vali[j]*xxi[i];
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qj+j*nconv,&one,&sone,tt,&one));
      PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2,&one));
    }
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      y2i = arrayi2+nloc;
      PetscCall(PetscArrayzero(y2i,nconv));
      for (j=0;j<pjd->midx;j++) {
        for (i=0;i<nconv;i++) tt[i] = val[j]*xxi[i]+vali[j]*xx[i];
        PetscCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qji+j*nconv,&one,&sone,tt,&one));
        PetscCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2i,&one));
      }
      for (i=0;i<nconv;i++) arrayi2[nloc+i] /= PetscSqrtReal(np);
    }
#endif
    for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
    PetscCall(PetscFree5(tt,x2,qj,xx,val));
  }
  PetscCall(VecResetArray(tx));
  PetscCall(VecRestoreArray(xs[0],&array1));
  PetscCall(VecResetArray(ty));
  PetscCall(VecRestoreArray(ys[0],&array2));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    PetscCall(VecResetArray(txi));
    PetscCall(VecRestoreArray(xs[1],&arrayi1));
    PetscCall(VecResetArray(tyi));
    PetscCall(VecRestoreArray(ys[1],&arrayi2));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatCreateVecs_PEPJD(Mat A,Vec *right,Vec *left)
{
  PEP_JD_MATSHELL *matctx;
  PEP_JD          *pjd;
  PetscInt        kspsf=1,i;
  Vec             v[2];

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&matctx));
  pjd   = (PEP_JD*)matctx->pep->data;
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  for (i=0;i<kspsf;i++) PetscCall(BVCreateVec(pjd->V,v+i));
  if (right) PetscCall(VecCreateCompWithVecs(v,kspsf,pjd->vtempl,right));
  if (left) PetscCall(VecCreateCompWithVecs(v,kspsf,pjd->vtempl,left));
  for (i=0;i<kspsf;i++) PetscCall(VecDestroy(&v[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDUpdateExtendedPC(PEP pep,PetscScalar theta)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL *pcctx;
  PetscInt       i,j,k,n=pjd->nlock,ld=pjd->ld,deg=pep->nmat-1;
  PetscScalar    *M,*ps,*work,*U,*V,*S,*Sp,*Spp,snone=-1.0,sone=1.0,zero=0.0,*val;
  PetscReal      tol,maxeig=0.0,*sg,*rwork;
  PetscBLASInt   n_,info,ld_,*p,lw_,rk=0;

  PetscFunctionBegin;
  if (n) {
    PetscCall(PCShellGetContext(pjd->pcshell,&pcctx));
    pcctx->theta = theta;
    pcctx->n = n;
    M  = pcctx->M;
    PetscCall(PetscBLASIntCast(n,&n_));
    PetscCall(PetscBLASIntCast(ld,&ld_));
    PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
    if (pjd->midx==1) {
      PetscCall(PetscArraycpy(M,pjd->XpX,ld*ld));
      PetscCall(PetscCalloc2(10*n,&work,n,&p));
    } else {
      ps = pcctx->ps;
      PetscCall(PetscCalloc7(2*n*n,&U,3*n*n,&S,n,&sg,10*n,&work,5*n,&rwork,n,&p,deg+1,&val));
      V = U+n*n;
      /* pseudo-inverse */
      for (j=0;j<n;j++) {
        for (i=0;i<n;i++) S[n*j+i] = -pjd->T[pep->ncv*j+i];
        S[n*j+j] += theta;
      }
      lw_ = 10*n_;
#if !defined (PETSC_USE_COMPLEX)
      PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&n_,&n_,S,&n_,sg,U,&n_,V,&n_,work,&lw_,&info));
#else
      PetscCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&n_,&n_,S,&n_,sg,U,&n_,V,&n_,work,&lw_,rwork,&info));
#endif
      SlepcCheckLapackInfo("gesvd",info);
      for (i=0;i<n;i++) maxeig = PetscMax(maxeig,sg[i]);
      tol = 10*PETSC_MACHINE_EPSILON*n*maxeig;
      for (j=0;j<n;j++) {
        if (sg[j]>tol) {
          for (i=0;i<n;i++) U[j*n+i] /= sg[j];
          rk++;
        } else break;
      }
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&rk,&sone,U,&n_,V,&n_,&zero,ps,&ld_));

      /* compute M */
      PetscCall(PEPEvaluateBasis(pep,theta,0.0,val,NULL));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&snone,pjd->XpX,&ld_,ps,&ld_,&zero,M,&ld_));
      PetscCall(PetscArrayzero(S,2*n*n));
      Sp = S+n*n;
      for (j=0;j<n;j++) S[j*(n+1)] = 1.0;
      for (k=1;k<pjd->midx;k++) {
        for (j=0;j<n;j++) for (i=0;i<n;i++) V[j*n+i] = S[j*n+i] - ps[j*ld+i]*val[k];
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,pjd->XpX,&ld_,V,&n_,&zero,U,&n_));
        PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,pjd->Tj+k*ld*ld,&ld_,U,&n_,&sone,M,&ld_));
        Spp = Sp; Sp = S;
        PetscCall(PEPJDEvaluateHatBasis(pep,n,pjd->T,pep->ncv,val,NULL,k+1,Spp,Sp,S));
      }
    }
    /* inverse */
    PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,M,&ld_,p,&info));
    SlepcCheckLapackInfo("getrf",info);
    PetscCallBLAS("LAPACKgetri",LAPACKgetri_(&n_,M,&ld_,p,work,&n_,&info));
    SlepcCheckLapackInfo("getri",info);
    PetscCall(PetscFPTrapPop());
    if (pjd->midx==1) PetscCall(PetscFree2(work,p));
    else PetscCall(PetscFree7(U,S,sg,work,rwork,p,val));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDMatSetUp(PEP pep,PetscInt sz,PetscScalar *theta)
{
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PEP_JD_MATSHELL *matctx;
  PEP_JD_PCSHELL  *pcctx;
  MatStructure    str;
  PetscScalar     *vals,*valsi;
  PetscBool       skipmat=PETSC_FALSE;
  PetscInt        i;
  Mat             Pr=NULL;

  PetscFunctionBegin;
  if (sz==2 && theta[1]==0.0) sz = 1;
  PetscCall(MatShellGetContext(pjd->Pshell,&matctx));
  PetscCall(PCShellGetContext(pjd->pcshell,&pcctx));
  if (matctx->Pr && matctx->theta[0]==theta[0] && ((!matctx->Pi && sz==1) || (sz==2 && matctx->theta[1]==theta[1]))) {
    if (pcctx->n == pjd->nlock) PetscFunctionReturn(PETSC_SUCCESS);
    skipmat = PETSC_TRUE;
  }
  if (!skipmat) {
    PetscCall(PetscMalloc2(pep->nmat,&vals,pep->nmat,&valsi));
    PetscCall(STGetMatStructure(pep->st,&str));
    PetscCall(PEPEvaluateBasis(pep,theta[0],theta[1],vals,valsi));
    if (!matctx->Pr) PetscCall(MatDuplicate(pep->A[0],MAT_COPY_VALUES,&matctx->Pr));
    else PetscCall(MatCopy(pep->A[0],matctx->Pr,str));
    for (i=1;i<pep->nmat;i++) PetscCall(MatAXPY(matctx->Pr,vals[i],pep->A[i],str));
    if (!pjd->reusepc) {
      if (pcctx->PPr && sz==2) {
        PetscCall(MatCopy(matctx->Pr,pcctx->PPr,str));
        Pr = pcctx->PPr;
      } else Pr = matctx->Pr;
    }
    matctx->theta[0] = theta[0];
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      if (!matctx->Pi) PetscCall(MatDuplicate(pep->A[0],MAT_COPY_VALUES,&matctx->Pi));
      else PetscCall(MatCopy(pep->A[1],matctx->Pi,str));
      PetscCall(MatScale(matctx->Pi,valsi[1]));
      for (i=2;i<pep->nmat;i++) PetscCall(MatAXPY(matctx->Pi,valsi[i],pep->A[i],str));
      matctx->theta[1] = theta[1];
    }
#endif
    PetscCall(PetscFree2(vals,valsi));
  }
  if (!pjd->reusepc) {
    if (!skipmat) {
      PetscCall(PCSetOperators(pcctx->pc,Pr,Pr));
      PetscCall(PCSetUp(pcctx->pc));
    }
    PetscCall(PEPJDUpdateExtendedPC(pep,theta[0]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDCreateShellPC(PEP pep,Vec *ww)
{
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL  *pcctx;
  PEP_JD_MATSHELL *matctx;
  KSP             ksp;
  PetscInt        nloc,mloc,kspsf=1;
  Vec             v[2];
  PetscScalar     target[2];
  Mat             Pr;

  PetscFunctionBegin;
  /* Create the reference vector */
  PetscCall(BVGetColumn(pjd->V,0,&v[0]));
  v[1] = v[0];
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  PetscCall(VecCreateCompWithVecs(v,kspsf,NULL,&pjd->vtempl));
  PetscCall(BVRestoreColumn(pjd->V,0,&v[0]));

  /* Replace preconditioner with one containing projectors */
  PetscCall(PCCreate(PetscObjectComm((PetscObject)pep),&pjd->pcshell));
  PetscCall(PCSetType(pjd->pcshell,PCSHELL));
  PetscCall(PCShellSetName(pjd->pcshell,"PCPEPJD"));
  PetscCall(PCShellSetApply(pjd->pcshell,PCShellApply_PEPJD));
  PetscCall(PetscNew(&pcctx));
  PetscCall(PCShellSetContext(pjd->pcshell,pcctx));
  PetscCall(STGetKSP(pep->st,&ksp));
  PetscCall(BVCreateVec(pjd->V,&pcctx->Bp[0]));
  PetscCall(VecDuplicate(pcctx->Bp[0],&pcctx->Bp[1]));
  PetscCall(KSPGetPC(ksp,&pcctx->pc));
  PetscCall(PetscObjectReference((PetscObject)pcctx->pc));
  PetscCall(MatGetLocalSize(pep->A[0],&mloc,&nloc));
  if (pjd->ld>1) {
    nloc += pjd->ld-1; mloc += pjd->ld-1;
  }
  PetscCall(PetscNew(&matctx));
  PetscCall(MatCreateShell(PetscObjectComm((PetscObject)pep),kspsf*nloc,kspsf*mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&pjd->Pshell));
  PetscCall(MatShellSetOperation(pjd->Pshell,MATOP_MULT,(void(*)(void))MatMult_PEPJD));
  PetscCall(MatShellSetOperation(pjd->Pshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_PEPJD));
  matctx->pep = pep;
  target[0] = pep->target; target[1] = 0.0;
  PetscCall(PEPJDMatSetUp(pep,1,target));
  Pr = matctx->Pr;
  pcctx->PPr = NULL;
#if !defined(PETSC_USE_COMPLEX)
  if (!pjd->reusepc) {
    PetscCall(MatDuplicate(matctx->Pr,MAT_COPY_VALUES,&pcctx->PPr));
    Pr = pcctx->PPr;
  }
#endif
  PetscCall(PCSetOperators(pcctx->pc,Pr,Pr));
  PetscCall(PCSetErrorIfFailure(pcctx->pc,PETSC_TRUE));
  PetscCall(KSPSetPC(ksp,pjd->pcshell));
  if (pjd->reusepc) {
    PetscCall(PCSetReusePreconditioner(pcctx->pc,PETSC_TRUE));
    PetscCall(KSPSetReusePreconditioner(ksp,PETSC_TRUE));
  }
  PetscCall(PEP_KSPSetOperators(ksp,pjd->Pshell,pjd->Pshell));
  PetscCall(KSPSetUp(ksp));
  if (pjd->ld>1) {
    PetscCall(PetscMalloc2(pjd->ld*pjd->ld,&pcctx->M,pjd->ld*pjd->ld,&pcctx->ps));
    pcctx->pep = pep;
  }
  matctx->work = ww;
  pcctx->work  = ww;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDEigenvectors(PEP pep)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBLASInt   ld,nconv,info,nc;
  PetscScalar    *Z;
  PetscReal      *wr;
  Mat            U;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar    *w;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(pep->ncv,&ld));
  PetscCall(PetscBLASIntCast(pep->nconv,&nconv));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscMalloc2(pep->nconv*pep->nconv,&Z,3*pep->ncv,&wr));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","A",NULL,&nconv,pjd->T,&ld,NULL,&nconv,Z,&nconv,&nconv,&nc,wr,&info));
#else
  PetscCall(PetscMalloc3(pep->nconv*pep->nconv,&Z,3*pep->ncv,&wr,2*pep->ncv,&w));
  PetscCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","A",NULL,&nconv,pjd->T,&ld,NULL,&nconv,Z,&nconv,&nconv,&nc,w,wr,&info));
#endif
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("trevc",info);
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nconv,nconv,Z,&U));
  PetscCall(BVSetActiveColumns(pjd->X,0,pep->nconv));
  PetscCall(BVMultInPlace(pjd->X,U,0,pep->nconv));
  PetscCall(BVNormalize(pjd->X,pep->eigi));
  PetscCall(MatDestroy(&U));
#if !defined(PETSC_USE_COMPLEX)
  PetscCall(PetscFree2(Z,wr));
#else
  PetscCall(PetscFree3(Z,wr,w));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDLockConverged(PEP pep,PetscInt *nv,PetscInt sz)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       j,i,*P,ldds,rk=0,nvv=*nv;
  Vec            v,x,w;
  PetscScalar    *R,*r,*pX,target[2];
  Mat            X;
  PetscBLASInt   sz_,rk_,nv_,info;
  PetscMPIInt    np;

  PetscFunctionBegin;
  /* update AX and XpX */
  for (i=sz;i>0;i--) {
    PetscCall(BVGetColumn(pjd->X,pjd->nlock-i,&x));
    for (j=0;j<pep->nmat;j++) {
      PetscCall(BVGetColumn(pjd->AX[j],pjd->nlock-i,&v));
      PetscCall(MatMult(pep->A[j],x,v));
      PetscCall(BVRestoreColumn(pjd->AX[j],pjd->nlock-i,&v));
      PetscCall(BVSetActiveColumns(pjd->AX[j],0,pjd->nlock-i+1));
    }
    PetscCall(BVRestoreColumn(pjd->X,pjd->nlock-i,&x));
    PetscCall(BVDotColumn(pjd->X,(pjd->nlock-i),pjd->XpX+(pjd->nlock-i)*pjd->ld));
    pjd->XpX[(pjd->nlock-i)*(1+pjd->ld)] = 1.0;
    for (j=0;j<pjd->nlock-i;j++) pjd->XpX[j*pjd->ld+pjd->nlock-i] = PetscConj(pjd->XpX[(pjd->nlock-i)*pjd->ld+j]);
  }

  /* minimality index */
  pjd->midx = PetscMin(pjd->mmidx,pjd->nlock);

  /* evaluate the polynomial basis in T */
  PetscCall(PetscArrayzero(pjd->Tj,pjd->ld*pjd->ld*pep->nmat));
  for (j=0;j<pep->nmat;j++) PetscCall(PEPEvaluateBasisMat(pep,pjd->nlock,pjd->T,pep->ncv,j,(j>1)?pjd->Tj+(j-2)*pjd->ld*pjd->ld:NULL,pjd->ld,j?pjd->Tj+(j-1)*pjd->ld*pjd->ld:NULL,pjd->ld,pjd->Tj+j*pjd->ld*pjd->ld,pjd->ld));

  /* Extend search space */
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
  PetscCall(PetscCalloc3(nvv,&P,nvv*nvv,&R,nvv*sz,&r));
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));
  PetscCall(DSGetArray(pep->ds,DS_MAT_X,&pX));
  PetscCall(PEPJDOrthogonalize(nvv,nvv,pX,ldds,&rk,P,R,nvv));
  for (j=0;j<sz;j++) {
    for (i=0;i<rk;i++) r[i*sz+j] = PetscConj(R[nvv*i+j]*pep->eigr[P[i]]); /* first row scaled with permuted diagonal */
  }
  PetscCall(PetscBLASIntCast(rk,&rk_));
  PetscCall(PetscBLASIntCast(sz,&sz_));
  PetscCall(PetscBLASIntCast(nvv,&nv_));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&rk_,R,&nv_,&info));
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("trtri",info);
  for (i=0;i<sz;i++) PetscCallBLAS("BLAStrmv",BLAStrmv_("U","C","N",&rk_,R,&nv_,r+i,&sz_));
  for (i=0;i<sz*rk;i++) r[i] = PetscConj(r[i])/PetscSqrtReal(np); /* revert */
  PetscCall(BVSetActiveColumns(pjd->V,0,nvv));
  rk -= sz;
  for (j=0;j<rk;j++) PetscCall(PetscArraycpy(R+j*nvv,pX+(j+sz)*ldds,nvv));
  PetscCall(DSRestoreArray(pep->ds,DS_MAT_X,&pX));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,nvv,rk,R,&X));
  PetscCall(BVMultInPlace(pjd->V,X,0,rk));
  PetscCall(MatDestroy(&X));
  PetscCall(BVSetActiveColumns(pjd->V,0,rk));
  for (j=0;j<rk;j++) {
    PetscCall(BVGetColumn(pjd->V,j,&v));
    PetscCall(PEPJDCopyToExtendedVec(pep,NULL,r+sz*(j+sz),sz,pjd->nlock-sz,v,PETSC_FALSE));
    PetscCall(BVRestoreColumn(pjd->V,j,&v));
  }
  PetscCall(BVOrthogonalize(pjd->V,NULL));

  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
    for (j=0;j<rk;j++) {
      /* W = P(target)*V */
      PetscCall(BVGetColumn(pjd->W,j,&w));
      PetscCall(BVGetColumn(pjd->V,j,&v));
      target[0] = pep->target; target[1] = 0.0;
      PetscCall(PEPJDComputeResidual(pep,PETSC_FALSE,1,&v,target,&w,pep->work));
      PetscCall(BVRestoreColumn(pjd->V,j,&v));
      PetscCall(BVRestoreColumn(pjd->W,j,&w));
    }
    PetscCall(BVSetActiveColumns(pjd->W,0,rk));
    PetscCall(BVOrthogonalize(pjd->W,NULL));
  }
  *nv = rk;
  PetscCall(PetscFree3(P,R,r));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSystemSetUp(PEP pep,PetscInt sz,PetscScalar *theta,Vec *u,Vec *p,Vec *ww)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL *pcctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    s[2];
#endif

  PetscFunctionBegin;
  PetscCall(PCShellGetContext(pjd->pcshell,&pcctx));
  PetscCall(PEPJDMatSetUp(pep,sz,theta));
  pcctx->u[0] = u[0]; pcctx->u[1] = u[1];
  /* Compute r'. p is a work space vector */
  PetscCall(PEPJDComputeResidual(pep,PETSC_TRUE,sz,u,theta,p,ww));
  PetscCall(PEPJDExtendedPCApply(pjd->pcshell,p[0],pcctx->Bp[0]));
  PetscCall(VecDot(pcctx->Bp[0],u[0],pcctx->gamma));
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    PetscCall(PEPJDExtendedPCApply(pjd->pcshell,p[1],pcctx->Bp[1]));
    PetscCall(VecDot(pcctx->Bp[0],u[1],pcctx->gamma+1));
    PetscCall(VecMDot(pcctx->Bp[1],2,u,s));
    pcctx->gamma[0] += s[1];
    pcctx->gamma[1] = -pcctx->gamma[1]+s[0];
  }
#endif
  if (sz==1) {
    PetscCall(VecZeroEntries(pcctx->Bp[1]));
    pcctx->gamma[1] = 0.0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSolve_JD(PEP pep)
{
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PetscInt        k,nv,nvc,ld,minv,dim,bupdated=0,sz=1,kspsf=1,idx,off,maxits,nloc;
  PetscMPIInt     np,count;
  PetscScalar     theta[2]={0.0,0.0},ritz[2]={0.0,0.0},*pX,*eig,*eigi,*array;
  PetscReal       norm,*res,tol=0.0,rtol,abstol, dtol;
  PetscBool       lindep,ini=PETSC_TRUE;
  Vec             tc,t[2]={NULL,NULL},u[2]={NULL,NULL},p[2]={NULL,NULL};
  Vec             rc,rr[2],r[2]={NULL,NULL},*ww=pep->work,v[2];
  Mat             G,X,Y;
  KSP             ksp;
  PEP_JD_PCSHELL  *pcctx;
  PEP_JD_MATSHELL *matctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal       norm1;
#endif

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np));
  PetscCall(BVGetSizes(pep->V,&nloc,NULL,NULL));
  PetscCall(DSGetLeadingDimension(pep->ds,&ld));
  PetscCall(PetscCalloc3(pep->ncv+pep->nev,&eig,pep->ncv+pep->nev,&eigi,pep->ncv+pep->nev,&res));
  pjd->nlock = 0;
  PetscCall(STGetKSP(pep->st,&ksp));
  PetscCall(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  PetscCall(PEPJDProcessInitialSpace(pep,ww));
  nv = (pep->nini)?pep->nini:1;

  /* Replace preconditioner with one containing projectors */
  PetscCall(PEPJDCreateShellPC(pep,ww));
  PetscCall(PCShellGetContext(pjd->pcshell,&pcctx));

  /* Create auxiliary vectors */
  PetscCall(BVCreateVec(pjd->V,&u[0]));
  PetscCall(VecDuplicate(u[0],&p[0]));
  PetscCall(VecDuplicate(u[0],&r[0]));
#if !defined (PETSC_USE_COMPLEX)
  PetscCall(VecDuplicate(u[0],&u[1]));
  PetscCall(VecDuplicate(u[0],&p[1]));
  PetscCall(VecDuplicate(u[0],&r[1]));
#endif

  /* Restart loop */
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    PetscCall(DSSetDimensions(pep->ds,nv,0,0));
    PetscCall(BVSetActiveColumns(pjd->V,bupdated,nv));
    PetscCall(PEPJDUpdateTV(pep,bupdated,nv,ww));
    if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) PetscCall(BVSetActiveColumns(pjd->W,bupdated,nv));
    for (k=0;k<pep->nmat;k++) {
      PetscCall(BVSetActiveColumns(pjd->TV[k],bupdated,nv));
      PetscCall(DSGetMat(pep->ds,DSMatExtra[k],&G));
      PetscCall(BVMatProject(pjd->TV[k],NULL,pjd->W,G));
      PetscCall(DSRestoreMat(pep->ds,DSMatExtra[k],&G));
    }
    PetscCall(BVSetActiveColumns(pjd->V,0,nv));
    PetscCall(BVSetActiveColumns(pjd->W,0,nv));

    /* Solve projected problem */
    PetscCall(DSSetState(pep->ds,DS_STATE_RAW));
    PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
    PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
    PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));
    idx = 0;
    do {
      ritz[0] = pep->eigr[idx];
#if !defined(PETSC_USE_COMPLEX)
      ritz[1] = pep->eigi[idx];
      sz = (ritz[1]==0.0)?1:2;
#endif
      /* Compute Ritz vector u=V*X(:,1) */
      PetscCall(DSGetArray(pep->ds,DS_MAT_X,&pX));
      PetscCall(BVSetActiveColumns(pjd->V,0,nv));
      PetscCall(BVMultVec(pjd->V,1.0,0.0,u[0],pX+idx*ld));
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) PetscCall(BVMultVec(pjd->V,1.0,0.0,u[1],pX+(idx+1)*ld));
#endif
      PetscCall(DSRestoreArray(pep->ds,DS_MAT_X,&pX));
      PetscCall(PEPJDComputeResidual(pep,PETSC_FALSE,sz,u,ritz,r,ww));
      /* Check convergence */
      PetscCall(VecNorm(r[0],NORM_2,&norm));
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) {
        PetscCall(VecNorm(r[1],NORM_2,&norm1));
        norm = SlepcAbs(norm,norm1);
      }
#endif
      PetscCall((*pep->converged)(pep,ritz[0],ritz[1],norm,&pep->errest[pep->nconv],pep->convergedctx));
      if (sz==2) pep->errest[pep->nconv+1] = pep->errest[pep->nconv];
      if (ini) {
        tol = PetscMin(.1,pep->errest[pep->nconv]); ini = PETSC_FALSE;
      } else tol = PetscMin(pep->errest[pep->nconv],tol/2);
      PetscCall((*pep->stopping)(pep,pep->its,pep->max_it,(pep->errest[pep->nconv]<pep->tol)?pep->nconv+sz:pep->nconv,pep->nev,&pep->reason,pep->stoppingctx));
      if (pep->errest[pep->nconv]<pep->tol) {
        /* Ritz pair converged */
        ini = PETSC_TRUE;
        minv = PetscMin(nv,(PetscInt)(pjd->keep*pep->ncv));
        if (pjd->ld>1) {
          PetscCall(BVGetColumn(pjd->X,pep->nconv,&v[0]));
          PetscCall(PEPJDCopyToExtendedVec(pep,v[0],pjd->T+pep->ncv*pep->nconv,pjd->ld-1,0,u[0],PETSC_TRUE));
          PetscCall(BVRestoreColumn(pjd->X,pep->nconv,&v[0]));
          PetscCall(BVSetActiveColumns(pjd->X,0,pep->nconv+1));
          PetscCall(BVNormColumn(pjd->X,pep->nconv,NORM_2,&norm));
          PetscCall(BVScaleColumn(pjd->X,pep->nconv,1.0/norm));
          for (k=0;k<pep->nconv;k++) pjd->T[pep->ncv*pep->nconv+k] *= PetscSqrtReal(np)/norm;
          pjd->T[(pep->ncv+1)*pep->nconv] = ritz[0];
          eig[pep->nconv] = ritz[0];
          idx++;
#if !defined(PETSC_USE_COMPLEX)
          if (sz==2) {
            PetscCall(BVGetColumn(pjd->X,pep->nconv+1,&v[0]));
            PetscCall(PEPJDCopyToExtendedVec(pep,v[0],pjd->T+pep->ncv*(pep->nconv+1),pjd->ld-1,0,u[1],PETSC_TRUE));
            PetscCall(BVRestoreColumn(pjd->X,pep->nconv+1,&v[0]));
            PetscCall(BVSetActiveColumns(pjd->X,0,pep->nconv+2));
            PetscCall(BVNormColumn(pjd->X,pep->nconv+1,NORM_2,&norm1));
            PetscCall(BVScaleColumn(pjd->X,pep->nconv+1,1.0/norm1));
            for (k=0;k<pep->nconv;k++) pjd->T[pep->ncv*(pep->nconv+1)+k] *= PetscSqrtReal(np)/norm1;
            pjd->T[(pep->ncv+1)*(pep->nconv+1)] = ritz[0];
            pjd->T[(pep->ncv+1)*pep->nconv+1] = -ritz[1]*norm1/norm;
            pjd->T[(pep->ncv+1)*(pep->nconv+1)-1] = ritz[1]*norm/norm1;
            eig[pep->nconv+1] = ritz[0];
            eigi[pep->nconv] = ritz[1]; eigi[pep->nconv+1] = -ritz[1];
            idx++;
          }
#endif
        } else PetscCall(BVInsertVec(pep->V,pep->nconv,u[0]));
        pep->nconv += sz;
      }
    } while (pep->errest[pep->nconv]<pep->tol && pep->nconv<nv);

    if (pep->reason==PEP_CONVERGED_ITERATING) {
      nvc = nv;
      if (idx) {
        pjd->nlock +=idx;
        PetscCall(PEPJDLockConverged(pep,&nv,idx));
      }
      if (nv+sz>=pep->ncv-1) {
        /* Basis full, force restart */
        minv = PetscMin(nv,(PetscInt)(pjd->keep*pep->ncv));
        PetscCall(DSGetDimensions(pep->ds,&dim,NULL,NULL,NULL));
        PetscCall(DSGetArray(pep->ds,DS_MAT_X,&pX));
        PetscCall(PEPJDOrthogonalize(dim,minv,pX,ld,&minv,NULL,NULL,ld));
        PetscCall(DSRestoreArray(pep->ds,DS_MAT_X,&pX));
        PetscCall(DSGetArray(pep->ds,DS_MAT_Y,&pX));
        PetscCall(PEPJDOrthogonalize(dim,minv,pX,ld,&minv,NULL,NULL,ld));
        PetscCall(DSRestoreArray(pep->ds,DS_MAT_Y,&pX));
        PetscCall(DSGetMat(pep->ds,DS_MAT_X,&X));
        PetscCall(BVMultInPlace(pjd->V,X,0,minv));
        PetscCall(DSRestoreMat(pep->ds,DS_MAT_X,&X));
        if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
         PetscCall(DSGetMat(pep->ds,DS_MAT_Y,&Y));
         PetscCall(BVMultInPlace(pjd->W,Y,pep->nconv,minv));
         PetscCall(DSRestoreMat(pep->ds,DS_MAT_Y,&Y));
        }
        nv = minv;
        bupdated = 0;
      } else {
        if (!idx && pep->errest[pep->nconv]<pjd->fix) {theta[0] = ritz[0]; theta[1] = ritz[1];}
        else {theta[0] = pep->target; theta[1] = 0.0;}
        /* Update system mat */
        PetscCall(PEPJDSystemSetUp(pep,sz,theta,u,p,ww));
        /* Solve correction equation to expand basis */
        PetscCall(BVGetColumn(pjd->V,nv,&t[0]));
        rr[0] = r[0];
        if (sz==2) {
          PetscCall(BVGetColumn(pjd->V,nv+1,&t[1]));
          rr[1] = r[1];
        } else {
          t[1] = NULL;
          rr[1] = NULL;
        }
        PetscCall(VecCreateCompWithVecs(t,kspsf,pjd->vtempl,&tc));
        PetscCall(VecCreateCompWithVecs(rr,kspsf,pjd->vtempl,&rc));
        PetscCall(VecCompSetSubVecs(pjd->vtempl,sz,NULL));
        tol  = PetscMax(rtol,tol/2);
        PetscCall(KSPSetTolerances(ksp,tol,abstol,dtol,maxits));
        PetscCall(KSPSolve(ksp,rc,tc));
        PetscCall(VecDestroy(&tc));
        PetscCall(VecDestroy(&rc));
        PetscCall(VecGetArray(t[0],&array));
        PetscCall(PetscMPIIntCast(pep->nconv,&count));
        PetscCallMPI(MPI_Bcast(array+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
        PetscCall(VecRestoreArray(t[0],&array));
        PetscCall(BVRestoreColumn(pjd->V,nv,&t[0]));
        PetscCall(BVOrthogonalizeColumn(pjd->V,nv,NULL,&norm,&lindep));
        if (lindep || norm==0.0) {
          PetscCheck(sz!=1,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Linearly dependent continuation vector");
          off = 1;
        } else {
          off = 0;
          PetscCall(BVScaleColumn(pjd->V,nv,1.0/norm));
        }
#if !defined(PETSC_USE_COMPLEX)
        if (sz==2) {
          PetscCall(VecGetArray(t[1],&array));
          PetscCallMPI(MPI_Bcast(array+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep)));
          PetscCall(VecRestoreArray(t[1],&array));
          PetscCall(BVRestoreColumn(pjd->V,nv+1,&t[1]));
          if (off) PetscCall(BVCopyColumn(pjd->V,nv+1,nv));
          PetscCall(BVOrthogonalizeColumn(pjd->V,nv+1-off,NULL,&norm,&lindep));
          if (lindep || norm==0.0) {
            PetscCheck(off==0,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Linearly dependent continuation vector");
            off = 1;
          } else PetscCall(BVScaleColumn(pjd->V,nv+1-off,1.0/norm));
        }
#endif
        if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
          PetscCall(BVInsertVec(pjd->W,nv,r[0]));
          if (sz==2 && !off) PetscCall(BVInsertVec(pjd->W,nv+1,r[1]));
          PetscCall(BVOrthogonalizeColumn(pjd->W,nv,NULL,&norm,&lindep));
          PetscCheck(!lindep && norm>0.0,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Linearly dependent continuation vector");
          PetscCall(BVScaleColumn(pjd->W,nv,1.0/norm));
          if (sz==2 && !off) {
            PetscCall(BVOrthogonalizeColumn(pjd->W,nv+1,NULL,&norm,&lindep));
            PetscCheck(!lindep && norm>0.0,PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED,"Linearly dependent continuation vector");
            PetscCall(BVScaleColumn(pjd->W,nv+1,1.0/norm));
          }
        }
        bupdated = idx?0:nv;
        nv += sz-off;
      }
      for (k=0;k<nvc;k++) {
        eig[pep->nconv-idx+k] = pep->eigr[k];
#if !defined(PETSC_USE_COMPLEX)
        eigi[pep->nconv-idx+k] = pep->eigi[k];
#endif
      }
      PetscCall(PEPMonitor(pep,pep->its,pep->nconv,eig,eigi,pep->errest,pep->nconv+1));
    }
  }
  if (pjd->ld>1) {
    for (k=0;k<pep->nconv;k++) {
      pep->eigr[k] = eig[k];
      pep->eigi[k] = eigi[k];
    }
    if (pep->nconv>0) PetscCall(PEPJDEigenvectors(pep));
    PetscCall(PetscFree2(pcctx->M,pcctx->ps));
  }
  PetscCall(VecDestroy(&u[0]));
  PetscCall(VecDestroy(&r[0]));
  PetscCall(VecDestroy(&p[0]));
#if !defined (PETSC_USE_COMPLEX)
  PetscCall(VecDestroy(&u[1]));
  PetscCall(VecDestroy(&r[1]));
  PetscCall(VecDestroy(&p[1]));
#endif
  PetscCall(KSPSetTolerances(ksp,rtol,abstol,dtol,maxits));
  PetscCall(KSPSetPC(ksp,pcctx->pc));
  PetscCall(VecDestroy(&pcctx->Bp[0]));
  PetscCall(VecDestroy(&pcctx->Bp[1]));
  PetscCall(MatShellGetContext(pjd->Pshell,&matctx));
  PetscCall(MatDestroy(&matctx->Pr));
  PetscCall(MatDestroy(&matctx->Pi));
  PetscCall(MatDestroy(&pjd->Pshell));
  PetscCall(MatDestroy(&pcctx->PPr));
  PetscCall(PCDestroy(&pcctx->pc));
  PetscCall(PetscFree(pcctx));
  PetscCall(PetscFree(matctx));
  PetscCall(PCDestroy(&pjd->pcshell));
  PetscCall(PetscFree3(eig,eigi,res));
  PetscCall(VecDestroy(&pjd->vtempl));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSetRestart_JD(PEP pep,PetscReal keep)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (keep==(PetscReal)PETSC_DEFAULT || keep==(PetscReal)PETSC_DECIDE) pjd->keep = 0.5;
  else {
    PetscCheck(keep>=0.1 && keep<=0.9,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    pjd->keep = keep;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDSetRestart - Sets the restart parameter for the Jacobi-Davidson
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -pep_jd_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: PEPJDGetRestart()
@*/
PetscErrorCode PEPJDSetRestart(PEP pep,PetscReal keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,keep,2);
  PetscTryMethod(pep,"PEPJDSetRestart_C",(PEP,PetscReal),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDGetRestart_JD(PEP pep,PetscReal *keep)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *keep = pjd->keep;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDGetRestart - Gets the restart parameter used in the Jacobi-Davidson method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: PEPJDSetRestart()
@*/
PetscErrorCode PEPJDGetRestart(PEP pep,PetscReal *keep)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(keep,2);
  PetscUseMethod(pep,"PEPJDGetRestart_C",(PEP,PetscReal*),(pep,keep));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSetFix_JD(PEP pep,PetscReal fix)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (fix == (PetscReal)PETSC_DEFAULT || fix == (PetscReal)PETSC_DECIDE) pjd->fix = 0.01;
  else {
    PetscCheck(fix>=0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid fix value, must be >0");
    pjd->fix = fix;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDSetFix - Sets the threshold for changing the target in the correction
   equation.

   Logically Collective

   Input Parameters:
+  pep - the eigenproblem solver context
-  fix - threshold for changing the target

   Options Database Key:
.  -pep_jd_fix - the fix value

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: PEPJDGetFix()
@*/
PetscErrorCode PEPJDSetFix(PEP pep,PetscReal fix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,fix,2);
  PetscTryMethod(pep,"PEPJDSetFix_C",(PEP,PetscReal),(pep,fix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDGetFix_JD(PEP pep,PetscReal *fix)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *fix = pjd->fix;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDGetFix - Returns the threshold for changing the target in the correction
   equation.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  fix - threshold for changing the target

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: PEPJDSetFix()
@*/
PetscErrorCode PEPJDGetFix(PEP pep,PetscReal *fix)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(fix,2);
  PetscUseMethod(pep,"PEPJDGetFix_C",(PEP,PetscReal*),(pep,fix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSetReusePreconditioner_JD(PEP pep,PetscBool reusepc)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  pjd->reusepc = reusepc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDSetReusePreconditioner - Sets a flag indicating whether the preconditioner
   must be reused or not.

   Logically Collective

   Input Parameters:
+  pep     - the eigenproblem solver context
-  reusepc - the reuse flag

   Options Database Key:
.  -pep_jd_reuse_preconditioner - the reuse flag

   Note:
   The default value is False. If set to True, the preconditioner is built
   only at the beginning, using the target value. Otherwise, it may be rebuilt
   (depending on the fix parameter) at each iteration from the Ritz value.

   Level: advanced

.seealso: PEPJDGetReusePreconditioner(), PEPJDSetFix()
@*/
PetscErrorCode PEPJDSetReusePreconditioner(PEP pep,PetscBool reusepc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,reusepc,2);
  PetscTryMethod(pep,"PEPJDSetReusePreconditioner_C",(PEP,PetscBool),(pep,reusepc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDGetReusePreconditioner_JD(PEP pep,PetscBool *reusepc)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *reusepc = pjd->reusepc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDGetReusePreconditioner - Returns the flag for reusing the preconditioner.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  reusepc - the reuse flag

   Level: advanced

.seealso: PEPJDSetReusePreconditioner()
@*/
PetscErrorCode PEPJDGetReusePreconditioner(PEP pep,PetscBool *reusepc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(reusepc,2);
  PetscUseMethod(pep,"PEPJDGetReusePreconditioner_C",(PEP,PetscBool*),(pep,reusepc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSetMinimalityIndex_JD(PEP pep,PetscInt mmidx)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (mmidx == PETSC_DEFAULT || mmidx == PETSC_DECIDE) {
    if (pjd->mmidx != 1) pep->state = PEP_STATE_INITIAL;
    pjd->mmidx = 1;
  } else {
    PetscCheck(mmidx>0,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmidx value, should be >0");
    if (pjd->mmidx != mmidx) pep->state = PEP_STATE_INITIAL;
    pjd->mmidx = mmidx;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDSetMinimalityIndex - Sets the maximum allowed value for the minimality index.

   Logically Collective

   Input Parameters:
+  pep   - the eigenproblem solver context
-  mmidx - maximum minimality index

   Options Database Key:
.  -pep_jd_minimality_index - the minimality index value

   Note:
   The default value is equal to the degree of the polynomial. A smaller value
   can be used if the wanted eigenvectors are known to be linearly independent.

   Level: advanced

.seealso: PEPJDGetMinimalityIndex()
@*/
PetscErrorCode PEPJDSetMinimalityIndex(PEP pep,PetscInt mmidx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,mmidx,2);
  PetscTryMethod(pep,"PEPJDSetMinimalityIndex_C",(PEP,PetscInt),(pep,mmidx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDGetMinimalityIndex_JD(PEP pep,PetscInt *mmidx)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *mmidx = pjd->mmidx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDGetMinimalityIndex - Returns the maximum allowed value of the minimality
   index.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  mmidx - minimality index

   Level: advanced

.seealso: PEPJDSetMinimalityIndex()
@*/
PetscErrorCode PEPJDGetMinimalityIndex(PEP pep,PetscInt *mmidx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(mmidx,2);
  PetscUseMethod(pep,"PEPJDGetMinimalityIndex_C",(PEP,PetscInt*),(pep,mmidx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDSetProjection_JD(PEP pep,PEPJDProjection proj)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  switch (proj) {
    case PEP_JD_PROJECTION_HARMONIC:
    case PEP_JD_PROJECTION_ORTHOGONAL:
      if (pjd->proj != proj) {
        pep->state = PEP_STATE_INITIAL;
        pjd->proj = proj;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'proj' value");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDSetProjection - Sets the type of projection to be used in the Jacobi-Davidson solver.

   Logically Collective

   Input Parameters:
+  pep  - the eigenproblem solver context
-  proj - the type of projection

   Options Database Key:
.  -pep_jd_projection - the projection type, either orthogonal or harmonic

   Level: advanced

.seealso: PEPJDGetProjection()
@*/
PetscErrorCode PEPJDSetProjection(PEP pep,PEPJDProjection proj)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,proj,2);
  PetscTryMethod(pep,"PEPJDSetProjection_C",(PEP,PEPJDProjection),(pep,proj));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPJDGetProjection_JD(PEP pep,PEPJDProjection *proj)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *proj = pjd->proj;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPJDGetProjection - Returns the type of projection used by the Jacobi-Davidson solver.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  proj - the type of projection

   Level: advanced

.seealso: PEPJDSetProjection()
@*/
PetscErrorCode PEPJDGetProjection(PEP pep,PEPJDProjection *proj)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscAssertPointer(proj,2);
  PetscUseMethod(pep,"PEPJDGetProjection_C",(PEP,PEPJDProjection*),(pep,proj));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetFromOptions_JD(PEP pep,PetscOptionItems PetscOptionsObject)
{
  PetscBool       flg,b1;
  PetscReal       r1;
  PetscInt        i1;
  PEPJDProjection proj;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"PEP JD Options");

    PetscCall(PetscOptionsReal("-pep_jd_restart","Proportion of vectors kept after restart","PEPJDSetRestart",0.5,&r1,&flg));
    if (flg) PetscCall(PEPJDSetRestart(pep,r1));

    PetscCall(PetscOptionsReal("-pep_jd_fix","Tolerance for changing the target in the correction equation","PEPJDSetFix",0.01,&r1,&flg));
    if (flg) PetscCall(PEPJDSetFix(pep,r1));

    PetscCall(PetscOptionsBool("-pep_jd_reuse_preconditioner","Whether to reuse the preconditioner","PEPJDSetReusePreconditoiner",PETSC_FALSE,&b1,&flg));
    if (flg) PetscCall(PEPJDSetReusePreconditioner(pep,b1));

    PetscCall(PetscOptionsInt("-pep_jd_minimality_index","Maximum allowed minimality index","PEPJDSetMinimalityIndex",1,&i1,&flg));
    if (flg) PetscCall(PEPJDSetMinimalityIndex(pep,i1));

    PetscCall(PetscOptionsEnum("-pep_jd_projection","Type of projection","PEPJDSetProjection",PEPJDProjectionTypes,(PetscEnum)PEP_JD_PROJECTION_HARMONIC,(PetscEnum*)&proj,&flg));
    if (flg) PetscCall(PEPJDSetProjection(pep,proj));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPView_JD(PEP pep,PetscViewer viewer)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*pjd->keep)));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  threshold for changing the target in the correction equation (fix): %g\n",(double)pjd->fix));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  projection type: %s\n",PEPJDProjectionTypes[pjd->proj]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  maximum allowed minimality index: %" PetscInt_FMT "\n",pjd->mmidx));
    if (pjd->reusepc) PetscCall(PetscViewerASCIIPrintf(viewer,"  reusing the preconditioner\n"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSetDefaultST_JD(PEP pep)
{
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)pep->st)->type_name) {
    PetscCall(STSetType(pep->st,STPRECOND));
    PetscCall(STPrecondSetKSPHasMat(pep->st,PETSC_TRUE));
  }
  PetscCall(STSetTransform(pep->st,PETSC_FALSE));
  PetscCall(STGetKSP(pep->st,&ksp));
  if (!((PetscObject)ksp)->type_name) {
    PetscCall(KSPSetType(ksp,KSPBCGSL));
    PetscCall(KSPSetTolerances(ksp,1e-5,PETSC_CURRENT,PETSC_CURRENT,100));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPReset_JD(PEP pep)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<pep->nmat;i++) PetscCall(BVDestroy(pjd->TV+i));
  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) PetscCall(BVDestroy(&pjd->W));
  if (pjd->ld>1) {
    PetscCall(BVDestroy(&pjd->V));
    for (i=0;i<pep->nmat;i++) PetscCall(BVDestroy(pjd->AX+i));
    PetscCall(BVDestroy(&pjd->N[0]));
    PetscCall(BVDestroy(&pjd->N[1]));
    PetscCall(PetscFree3(pjd->XpX,pjd->T,pjd->Tj));
  }
  PetscCall(PetscFree2(pjd->TV,pjd->AX));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPDestroy_JD(PEP pep)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(pep->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetRestart_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetFix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetFix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetReusePreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetReusePreconditioner_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetMinimalityIndex_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetMinimalityIndex_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetProjection_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetProjection_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_JD(PEP pep)
{
  PEP_JD         *pjd;

  PetscFunctionBegin;
  PetscCall(PetscNew(&pjd));
  pep->data = (void*)pjd;

  pep->lineariz = PETSC_FALSE;
  pjd->fix      = 0.01;
  pjd->mmidx    = 0;

  pep->ops->solve          = PEPSolve_JD;
  pep->ops->setup          = PEPSetUp_JD;
  pep->ops->setfromoptions = PEPSetFromOptions_JD;
  pep->ops->destroy        = PEPDestroy_JD;
  pep->ops->reset          = PEPReset_JD;
  pep->ops->view           = PEPView_JD;
  pep->ops->setdefaultst   = PEPSetDefaultST_JD;

  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetRestart_C",PEPJDSetRestart_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetRestart_C",PEPJDGetRestart_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetFix_C",PEPJDSetFix_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetFix_C",PEPJDGetFix_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetReusePreconditioner_C",PEPJDSetReusePreconditioner_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetReusePreconditioner_C",PEPJDGetReusePreconditioner_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetMinimalityIndex_C",PEPJDSetMinimalityIndex_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetMinimalityIndex_C",PEPJDGetMinimalityIndex_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetProjection_C",PEPJDSetProjection_JD));
  PetscCall(PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetProjection_C",PEPJDGetProjection_JD));
  PetscFunctionReturn(PETSC_SUCCESS);
}
