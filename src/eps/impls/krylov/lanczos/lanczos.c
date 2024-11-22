/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lanczos"

   Method: Explicitly Restarted Symmetric/Hermitian Lanczos

   Algorithm:

       Lanczos method for symmetric (Hermitian) problems, with explicit
       restart and deflation. Several reorthogonalization strategies can
       be selected.

   References:

       [1] "Lanczos Methods in SLEPc", SLEPc Technical Report STR-5,
           available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  EPSLanczosReorthogType reorthog;      /* user-provided reorthogonalization parameter */
  PetscInt               allocsize;     /* number of columns of work BV's allocated at setup */
  BV                     AV;            /* work BV used in selective reorthogonalization */
} EPS_LANCZOS;

static PetscErrorCode EPSSetUp_Lanczos(EPS eps)
{
  EPS_LANCZOS        *lanczos = (EPS_LANCZOS*)eps->data;
  BVOrthogRefineType refine;
  BVOrthogBlockType  btype;
  PetscReal          eta;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(EPSSetDimensions_Default(eps,&eps->nev,&eps->ncv,&eps->mpd));
  PetscCheck(eps->ncv<=eps->nev+eps->mpd,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues");
  EPSCheckUnsupported(eps,EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION | EPS_FEATURE_EXTRACTION | EPS_FEATURE_THRESHOLD);
  EPSCheckIgnored(eps,EPS_FEATURE_BALANCE);

  PetscCheck(lanczos->reorthog!=(EPSLanczosReorthogType)-1,PetscObjectComm((PetscObject)eps),PETSC_ERR_USER_INPUT,"You should explicitly provide the reorthogonalization type, e.g., -eps_lanczos_reorthog local. Note that the EPSLANCZOS solver is *NOT RECOMMENDED* for general use, because it uses explicit restart which typically has slow convergence. The recommended solver is EPSKRYLOVSCHUR (the default), which implements Lanczos with thick restart in the case of symmetric/Hermitian problems");

  PetscCall(EPSAllocateSolution(eps,1));
  PetscCall(EPS_SetInnerProduct(eps));
  if (lanczos->reorthog != EPS_LANCZOS_REORTHOG_FULL) {
    PetscCall(BVGetOrthogonalization(eps->V,NULL,&refine,&eta,&btype));
    PetscCall(BVSetOrthogonalization(eps->V,BV_ORTHOG_MGS,refine,eta,btype));
    PetscCall(PetscInfo(eps,"Switching to MGS orthogonalization\n"));
  }
  if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_SELECTIVE) {
    if (!lanczos->allocsize) {
      PetscCall(BVDuplicate(eps->V,&lanczos->AV));
      PetscCall(BVGetSizes(lanczos->AV,NULL,NULL,&lanczos->allocsize));
    } else { /* make sure V and AV have the same size */
      PetscCall(BVGetSizes(eps->V,NULL,NULL,&lanczos->allocsize));
      PetscCall(BVResize(lanczos->AV,lanczos->allocsize,PETSC_FALSE));
    }
  }

  PetscCall(DSSetType(eps->ds,DSHEP));
  PetscCall(DSSetCompact(eps->ds,PETSC_TRUE));
  PetscCall(DSAllocate(eps->ds,eps->ncv+1));
  if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) PetscCall(EPSSetWorkVecs(eps,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSLocalLanczos - Local reorthogonalization.

   This is the simplest variant. At each Lanczos step, the corresponding Lanczos vector
   is orthogonalized with respect to the two previous Lanczos vectors, according to
   the three term Lanczos recurrence. WARNING: This variant does not track the loss of
   orthogonality that occurs in finite-precision arithmetic and, therefore, the
   generated vectors are not guaranteed to be (semi-)orthogonal.
*/
static PetscErrorCode EPSLocalLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscInt       i,j,m = *M;
  Mat            Op;
  PetscBool      *which,lwhich[100];
  PetscScalar    *hwork,lhwork[100];

  PetscFunctionBegin;
  if (m > 100) PetscCall(PetscMalloc2(m,&which,m,&hwork));
  else {
    which = lwhich;
    hwork = lhwork;
  }
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;

  PetscCall(BVSetActiveColumns(eps->V,0,m));
  PetscCall(STGetOperator(eps->st,&Op));
  for (j=k;j<m;j++) {
    PetscCall(BVMatMultColumn(eps->V,Op,j));
    which[j] = PETSC_TRUE;
    if (j-2>=k) which[j-2] = PETSC_FALSE;
    PetscCall(BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,beta+j,breakdown));
    alpha[j] = PetscRealPart(hwork[j]);
    if (PetscUnlikely(*breakdown)) {
      *M = j+1;
      break;
    } else PetscCall(BVScaleColumn(eps->V,j+1,1/beta[j]));
  }
  PetscCall(STRestoreOperator(eps->st,&Op));
  if (m > 100) PetscCall(PetscFree2(which,hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   DenseTridiagonal - Solves a real tridiagonal Hermitian Eigenvalue Problem.

   Input Parameters:
+  n   - dimension of the eigenproblem
.  D   - pointer to the array containing the diagonal elements
-  E   - pointer to the array containing the off-diagonal elements

   Output Parameters:
+  w  - pointer to the array to store the computed eigenvalues
-  V  - pointer to the array to store the eigenvectors

   Notes:
   If V is NULL then the eigenvectors are not computed.

   This routine use LAPACK routines xSTEVR.
*/
static PetscErrorCode DenseTridiagonal(PetscInt n_,PetscReal *D,PetscReal *E,PetscReal *w,PetscScalar *V)
{
  PetscReal      abstol = 0.0,vl,vu,*work;
  PetscBLASInt   il,iu,m,*isuppz,n,lwork,*iwork,liwork,info;
  const char     *jobz;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       i,j;
  PetscReal      *VV=NULL;
#endif

  PetscFunctionBegin;
  PetscCall(PetscBLASIntCast(n_,&n));
  PetscCall(PetscBLASIntCast(20*n_,&lwork));
  PetscCall(PetscBLASIntCast(10*n_,&liwork));
  if (V) {
    jobz = "V";
#if defined(PETSC_USE_COMPLEX)
    PetscCall(PetscMalloc1(n*n,&VV));
#endif
  } else jobz = "N";
  PetscCall(PetscMalloc3(2*n,&isuppz,lwork,&work,liwork,&iwork));
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
#if defined(PETSC_USE_COMPLEX)
  PetscCallBLAS("LAPACKstevr",LAPACKstevr_(jobz,"A",&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,VV,&n,isuppz,work,&lwork,iwork,&liwork,&info));
#else
  PetscCallBLAS("LAPACKstevr",LAPACKstevr_(jobz,"A",&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,V,&n,isuppz,work,&lwork,iwork,&liwork,&info));
#endif
  PetscCall(PetscFPTrapPop());
  SlepcCheckLapackInfo("stevr",info);
#if defined(PETSC_USE_COMPLEX)
  if (V) {
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        V[i*n+j] = VV[i*n+j];
    PetscCall(PetscFree(VV));
  }
#endif
  PetscCall(PetscFree3(isuppz,work,iwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSSelectiveLanczos - Selective reorthogonalization.
*/
static PetscErrorCode EPSSelectiveLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscReal anorm)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscInt       i,j,m = *M,n,nritz=0,nritzo;
  Vec            vj1,av;
  Mat            Op;
  PetscReal      *d,*e,*ritz,norm;
  PetscScalar    *Y,*hwork;
  PetscBool      *which;

  PetscFunctionBegin;
  PetscCall(PetscCalloc6(m+1,&d,m,&e,m,&ritz,m*m,&Y,m,&which,m,&hwork));
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;
  PetscCall(STGetOperator(eps->st,&Op));

  for (j=k;j<m;j++) {
    PetscCall(BVSetActiveColumns(eps->V,0,m));

    /* Lanczos step */
    PetscCall(BVMatMultColumn(eps->V,Op,j));
    which[j] = PETSC_TRUE;
    if (j-2>=k) which[j-2] = PETSC_FALSE;
    PetscCall(BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,&norm,breakdown));
    alpha[j] = PetscRealPart(hwork[j]);
    beta[j] = norm;
    if (PetscUnlikely(*breakdown)) {
      *M = j+1;
      break;
    }

    /* Compute eigenvalues and eigenvectors Y of the tridiagonal block */
    n = j-k+1;
    for (i=0;i<n;i++) {
      d[i] = alpha[i+k];
      e[i] = beta[i+k];
    }
    PetscCall(DenseTridiagonal(n,d,e,ritz,Y));

    /* Estimate ||A|| */
    for (i=0;i<n;i++)
      if (PetscAbsReal(ritz[i]) > anorm) anorm = PetscAbsReal(ritz[i]);

    /* Compute nearly converged Ritz vectors */
    nritzo = 0;
    for (i=0;i<n;i++) {
      if (norm*PetscAbsScalar(Y[i*n+n-1]) < PETSC_SQRT_MACHINE_EPSILON*anorm) nritzo++;
    }
    if (nritzo>nritz) {
      nritz = 0;
      for (i=0;i<n;i++) {
        if (norm*PetscAbsScalar(Y[i*n+n-1]) < PETSC_SQRT_MACHINE_EPSILON*anorm) {
          PetscCall(BVSetActiveColumns(eps->V,k,k+n));
          PetscCall(BVGetColumn(lanczos->AV,nritz,&av));
          PetscCall(BVMultVec(eps->V,1.0,0.0,av,Y+i*n));
          PetscCall(BVRestoreColumn(lanczos->AV,nritz,&av));
          nritz++;
        }
      }
    }
    if (nritz > 0) {
      PetscCall(BVGetColumn(eps->V,j+1,&vj1));
      PetscCall(BVSetActiveColumns(lanczos->AV,0,nritz));
      PetscCall(BVOrthogonalizeVec(lanczos->AV,vj1,hwork,&norm,breakdown));
      PetscCall(BVRestoreColumn(eps->V,j+1,&vj1));
      if (PetscUnlikely(*breakdown)) {
        *M = j+1;
        break;
      }
    }
    PetscCall(BVScaleColumn(eps->V,j+1,1.0/norm));
  }

  PetscCall(STRestoreOperator(eps->st,&Op));
  PetscCall(PetscFree6(d,e,ritz,Y,which,hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static void update_omega(PetscReal *omega,PetscReal *omega_old,PetscInt j,PetscReal *alpha,PetscReal *beta,PetscReal eps1,PetscReal anorm)
{
  PetscInt  k;
  PetscReal T,binv;

  PetscFunctionBegin;
  /* Estimate of contribution to roundoff errors from A*v
       fl(A*v) = A*v + f,
     where ||f|| \approx eps1*||A||.
     For a full matrix A, a rule-of-thumb estimate is eps1 = sqrt(n)*eps */
  T = eps1*anorm;
  binv = 1.0/beta[j+1];

  /* Update omega(1) using omega(0)==0 */
  omega_old[0]= beta[1]*omega[1] + (alpha[0]-alpha[j])*omega[0] - beta[j]*omega_old[0];
  if (omega_old[0] > 0) omega_old[0] = binv*(omega_old[0] + T);
  else omega_old[0] = binv*(omega_old[0] - T);

  /* Update remaining components */
  for (k=1;k<j-1;k++) {
    omega_old[k] = beta[k+1]*omega[k+1] + (alpha[k]-alpha[j])*omega[k] + beta[k]*omega[k-1] - beta[j]*omega_old[k];
    if (omega_old[k] > 0) omega_old[k] = binv*(omega_old[k] + T);
    else omega_old[k] = binv*(omega_old[k] - T);
  }
  omega_old[j-1] = binv*T;

  /* Swap omega and omega_old */
  for (k=0;k<j;k++) {
    omega[k] = omega_old[k];
    omega_old[k] = omega[k];
  }
  omega[j] = eps1;
  PetscFunctionReturnVoid();
}

static void compute_int(PetscBool *which,PetscReal *mu,PetscInt j,PetscReal delta,PetscReal eta)
{
  PetscInt  i,k,maxpos;
  PetscReal max;
  PetscBool found;

  PetscFunctionBegin;
  /* initialize which */
  found = PETSC_FALSE;
  maxpos = 0;
  max = 0.0;
  for (i=0;i<j;i++) {
    if (PetscAbsReal(mu[i]) >= delta) {
      which[i] = PETSC_TRUE;
      found = PETSC_TRUE;
    } else which[i] = PETSC_FALSE;
    if (PetscAbsReal(mu[i]) > max) {
      maxpos = i;
      max = PetscAbsReal(mu[i]);
    }
  }
  if (!found) which[maxpos] = PETSC_TRUE;

  for (i=0;i<j;i++) {
    if (which[i]) {
      /* find left interval */
      for (k=i;k>=0;k--) {
        if (PetscAbsReal(mu[k])<eta || which[k]) break;
        else which[k] = PETSC_TRUE;
      }
      /* find right interval */
      for (k=i+1;k<j;k++) {
        if (PetscAbsReal(mu[k])<eta || which[k]) break;
        else which[k] = PETSC_TRUE;
      }
    }
  }
  PetscFunctionReturnVoid();
}

/*
   EPSPartialLanczos - Partial reorthogonalization.
*/
static PetscErrorCode EPSPartialLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscReal anorm)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscInt       i,j,m = *M;
  Mat            Op;
  PetscReal      norm,*omega,lomega[100],*omega_old,lomega_old[100],eps1,delta,eta;
  PetscBool      *which,lwhich[100],*which2,lwhich2[100];
  PetscBool      reorth = PETSC_FALSE,force_reorth = PETSC_FALSE;
  PetscBool      fro = PETSC_FALSE,estimate_anorm = PETSC_FALSE;
  PetscScalar    *hwork,lhwork[100];

  PetscFunctionBegin;
  if (m>100) PetscCall(PetscMalloc5(m,&omega,m,&omega_old,m,&which,m,&which2,m,&hwork));
  else {
    omega     = lomega;
    omega_old = lomega_old;
    which     = lwhich;
    which2    = lwhich2;
    hwork     = lhwork;
  }

  eps1 = PetscSqrtReal((PetscReal)eps->n)*PETSC_MACHINE_EPSILON/2;
  delta = PETSC_SQRT_MACHINE_EPSILON/PetscSqrtReal((PetscReal)eps->ncv);
  eta = PetscPowReal(PETSC_MACHINE_EPSILON,3.0/4.0)/PetscSqrtReal((PetscReal)eps->ncv);
  if (anorm < 0.0) {
    anorm = 1.0;
    estimate_anorm = PETSC_TRUE;
  }
  for (i=0;i<PetscMax(100,m);i++) omega[i] = omega_old[i] = 0.0;
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;

  PetscCall(BVSetActiveColumns(eps->V,0,m));
  PetscCall(STGetOperator(eps->st,&Op));
  for (j=k;j<m;j++) {
    PetscCall(BVMatMultColumn(eps->V,Op,j));
    if (fro) {
      /* Lanczos step with full reorthogonalization */
      PetscCall(BVOrthogonalizeColumn(eps->V,j+1,hwork,&norm,breakdown));
      alpha[j] = PetscRealPart(hwork[j]);
    } else {
      /* Lanczos step */
      which[j] = PETSC_TRUE;
      if (j-2>=k) which[j-2] = PETSC_FALSE;
      PetscCall(BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,&norm,breakdown));
      alpha[j] = PetscRealPart(hwork[j]);
      beta[j] = norm;

      /* Estimate ||A|| if needed */
      if (estimate_anorm) {
        if (j>k) anorm = PetscMax(anorm,PetscAbsReal(alpha[j])+norm+beta[j-1]);
        else anorm = PetscMax(anorm,PetscAbsReal(alpha[j])+norm);
      }

      /* Check if reorthogonalization is needed */
      reorth = PETSC_FALSE;
      if (j>k) {
        update_omega(omega,omega_old,j,alpha,beta-1,eps1,anorm);
        for (i=0;i<j-k;i++) {
          if (PetscAbsReal(omega[i]) > delta) reorth = PETSC_TRUE;
        }
      }
      if (reorth || force_reorth) {
        for (i=0;i<k;i++) which2[i] = PETSC_FALSE;
        for (i=k;i<=j;i++) which2[i] = PETSC_TRUE;
        if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_PERIODIC) {
          /* Periodic reorthogonalization */
          if (force_reorth) force_reorth = PETSC_FALSE;
          else force_reorth = PETSC_TRUE;
          for (i=0;i<j-k;i++) omega[i] = eps1;
        } else {
          /* Partial reorthogonalization */
          if (force_reorth) force_reorth = PETSC_FALSE;
          else {
            force_reorth = PETSC_TRUE;
            compute_int(which2+k,omega,j-k,delta,eta);
            for (i=0;i<j-k;i++) {
              if (which2[i+k]) omega[i] = eps1;
            }
          }
        }
        PetscCall(BVOrthogonalizeSomeColumn(eps->V,j+1,which2,hwork,&norm,breakdown));
      }
    }

    if (PetscUnlikely(*breakdown || norm < eps->n*anorm*PETSC_MACHINE_EPSILON)) {
      *M = j+1;
      break;
    }
    if (!fro && norm*delta < anorm*eps1) {
      fro = PETSC_TRUE;
      PetscCall(PetscInfo(eps,"Switching to full reorthogonalization at iteration %" PetscInt_FMT "\n",eps->its));
    }
    beta[j] = norm;
    PetscCall(BVScaleColumn(eps->V,j+1,1.0/norm));
  }

  PetscCall(STRestoreOperator(eps->st,&Op));
  if (m>100) PetscCall(PetscFree5(omega,omega_old,which,which2,hwork));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   EPSBasicLanczos - Computes an m-step Lanczos factorization. The first k
   columns are assumed to be locked and therefore they are not modified. On
   exit, the following relation is satisfied:

                    OP * V - V * T = f * e_m^T

   where the columns of V are the Lanczos vectors, T is a tridiagonal matrix,
   f is the residual vector and e_m is the m-th vector of the canonical basis.
   The Lanczos vectors (together with vector f) are B-orthogonal (to working
   accuracy) if full reorthogonalization is being used, otherwise they are
   (B-)semi-orthogonal. On exit, beta contains the B-norm of f and the next
   Lanczos vector can be computed as v_{m+1} = f / beta.

   This function simply calls another function which depends on the selected
   reorthogonalization strategy.
*/
static PetscErrorCode EPSBasicLanczos(EPS eps,PetscInt k,PetscInt *m,PetscReal *betam,PetscBool *breakdown,PetscReal anorm)
{
  EPS_LANCZOS        *lanczos = (EPS_LANCZOS*)eps->data;
  PetscScalar        *T;
  PetscInt           i,n=*m,ld;
  PetscReal          *alpha,*beta;
  BVOrthogRefineType orthog_ref;
  Mat                Op,M;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  switch (lanczos->reorthog) {
    case EPS_LANCZOS_REORTHOG_LOCAL:
      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
      beta = alpha + ld;
      PetscCall(EPSLocalLanczos(eps,alpha,beta,k,m,breakdown));
      *betam = beta[*m-1];
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
      break;
    case EPS_LANCZOS_REORTHOG_FULL:
      PetscCall(STGetOperator(eps->st,&Op));
      PetscCall(DSGetMat(eps->ds,DS_MAT_T,&M));
      PetscCall(BVMatLanczos(eps->V,Op,M,k,m,betam,breakdown));
      PetscCall(DSRestoreMat(eps->ds,DS_MAT_T,&M));
      PetscCall(STRestoreOperator(eps->st,&Op));
      break;
    case EPS_LANCZOS_REORTHOG_SELECTIVE:
      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
      beta = alpha + ld;
      PetscCall(EPSSelectiveLanczos(eps,alpha,beta,k,m,breakdown,anorm));
      *betam = beta[*m-1];
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
      break;
    case EPS_LANCZOS_REORTHOG_PERIODIC:
    case EPS_LANCZOS_REORTHOG_PARTIAL:
      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
      beta = alpha + ld;
      PetscCall(EPSPartialLanczos(eps,alpha,beta,k,m,breakdown,anorm));
      *betam = beta[*m-1];
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
      break;
    case EPS_LANCZOS_REORTHOG_DELAYED:
      PetscCall(PetscMalloc1(n*n,&T));
      PetscCall(BVGetOrthogonalization(eps->V,NULL,&orthog_ref,NULL,NULL));
      if (orthog_ref == BV_ORTHOG_REFINE_NEVER) PetscCall(EPSDelayedArnoldi1(eps,T,n,k,m,betam,breakdown));
      else PetscCall(EPSDelayedArnoldi(eps,T,n,k,m,betam,breakdown));
      n = *m;
      PetscCall(DSGetArrayReal(eps->ds,DS_MAT_T,&alpha));
      beta = alpha + ld;
      for (i=k;i<n-1;i++) {
        alpha[i] = PetscRealPart(T[n*i+i]);
        beta[i] = PetscRealPart(T[n*i+i+1]);
      }
      alpha[n-1] = PetscRealPart(T[n*(n-1)+n-1]);
      beta[n-1] = *betam;
      PetscCall(DSRestoreArrayReal(eps->ds,DS_MAT_T,&alpha));
      PetscCall(PetscFree(T));
      break;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_Lanczos(EPS eps)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscInt       nconv,i,j,k,l,x,n,*perm,restart,ncv=eps->ncv,r,ld;
  Vec            vi,vj,w;
  Mat            U;
  PetscScalar    *Y,*ritz,stmp;
  PetscReal      *bnd,anorm,beta,norm,rtmp,resnorm;
  PetscBool      breakdown;
  char           *conv,ctmp;

  PetscFunctionBegin;
  PetscCall(DSGetLeadingDimension(eps->ds,&ld));
  PetscCall(PetscMalloc4(ncv,&ritz,ncv,&bnd,ncv,&perm,ncv,&conv));

  /* The first Lanczos vector is the normalized initial vector */
  PetscCall(EPSGetStartVector(eps,0,NULL));

  anorm = -1.0;
  nconv = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an ncv-step Lanczos factorization */
    n = PetscMin(nconv+eps->mpd,ncv);
    PetscCall(DSSetDimensions(eps->ds,n,nconv,PETSC_DETERMINE));
    PetscCall(EPSBasicLanczos(eps,nconv,&n,&beta,&breakdown,anorm));
    PetscCall(DSSetDimensions(eps->ds,n,nconv,0));
    PetscCall(DSSetState(eps->ds,DS_STATE_INTERMEDIATE));
    PetscCall(BVSetActiveColumns(eps->V,nconv,n));

    /* Solve projected problem */
    PetscCall(DSSolve(eps->ds,ritz,NULL));
    PetscCall(DSSort(eps->ds,ritz,NULL,NULL,NULL,NULL));
    PetscCall(DSSynchronize(eps->ds,ritz,NULL));

    /* Estimate ||A|| */
    for (i=nconv;i<n;i++)
      anorm = PetscMax(anorm,PetscAbsReal(PetscRealPart(ritz[i])));

    /* Compute residual norm estimates as beta*abs(Y(m,:)) + eps*||A|| */
    PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Y));
    for (i=nconv;i<n;i++) {
      resnorm = beta*PetscAbsScalar(Y[n-1+i*ld]) + PETSC_MACHINE_EPSILON*anorm;
      PetscCall((*eps->converged)(eps,ritz[i],eps->eigi[i],resnorm,&bnd[i],eps->convergedctx));
      if (bnd[i]<eps->tol) conv[i] = 'C';
      else conv[i] = 'N';
    }
    PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Y));

    /* purge repeated ritz values */
    if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) {
      for (i=nconv+1;i<n;i++) {
        if (conv[i] == 'C' && PetscAbsScalar((ritz[i]-ritz[i-1])/ritz[i]) < eps->tol) conv[i] = 'R';
      }
    }

    /* Compute restart vector */
    if (breakdown) PetscCall(PetscInfo(eps,"Breakdown in Lanczos method (it=%" PetscInt_FMT " norm=%g)\n",eps->its,(double)beta));
    else {
      restart = nconv;
      while (restart<n && conv[restart] != 'N') restart++;
      if (restart >= n) {
        breakdown = PETSC_TRUE;
      } else {
        for (i=restart+1;i<n;i++) {
          if (conv[i] == 'N') {
            PetscCall(SlepcSCCompare(eps->sc,ritz[restart],0.0,ritz[i],0.0,&r));
            if (r>0) restart = i;
          }
        }
        PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Y));
        PetscCall(BVMultColumn(eps->V,1.0,0.0,n,Y+restart*ld+nconv));
        PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Y));
      }
    }

    /* Count and put converged eigenvalues first */
    for (i=nconv;i<n;i++) perm[i] = i;
    for (k=nconv;k<n;k++) {
      if (conv[perm[k]] != 'C') {
        j = k + 1;
        while (j<n && conv[perm[j]] != 'C') j++;
        if (j>=n) break;
        l = perm[k]; perm[k] = perm[j]; perm[j] = l;
      }
    }

    /* Sort eigenvectors according to permutation */
    PetscCall(DSGetArray(eps->ds,DS_MAT_Q,&Y));
    for (i=nconv;i<k;i++) {
      x = perm[i];
      if (x != i) {
        j = i + 1;
        while (perm[j] != i) j++;
        /* swap eigenvalues i and j */
        stmp = ritz[x]; ritz[x] = ritz[i]; ritz[i] = stmp;
        rtmp = bnd[x]; bnd[x] = bnd[i]; bnd[i] = rtmp;
        ctmp = conv[x]; conv[x] = conv[i]; conv[i] = ctmp;
        perm[j] = x; perm[i] = i;
        /* swap eigenvectors i and j */
        for (l=0;l<n;l++) {
          stmp = Y[l+x*ld]; Y[l+x*ld] = Y[l+i*ld]; Y[l+i*ld] = stmp;
        }
      }
    }
    PetscCall(DSRestoreArray(eps->ds,DS_MAT_Q,&Y));

    /* compute converged eigenvectors */
    PetscCall(DSGetMat(eps->ds,DS_MAT_Q,&U));
    PetscCall(BVMultInPlace(eps->V,U,nconv,k));
    PetscCall(DSRestoreMat(eps->ds,DS_MAT_Q,&U));

    /* purge spurious ritz values */
    if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) {
      for (i=nconv;i<k;i++) {
        PetscCall(BVGetColumn(eps->V,i,&vi));
        PetscCall(VecNorm(vi,NORM_2,&norm));
        PetscCall(VecScale(vi,1.0/norm));
        w = eps->work[0];
        PetscCall(STApply(eps->st,vi,w));
        PetscCall(VecAXPY(w,-ritz[i],vi));
        PetscCall(BVRestoreColumn(eps->V,i,&vi));
        PetscCall(VecNorm(w,NORM_2,&norm));
        PetscCall((*eps->converged)(eps,ritz[i],eps->eigi[i],norm,&bnd[i],eps->convergedctx));
        if (bnd[i]>=eps->tol) conv[i] = 'S';
      }
      for (i=nconv;i<k;i++) {
        if (conv[i] != 'C') {
          j = i + 1;
          while (j<k && conv[j] != 'C') j++;
          if (j>=k) break;
          /* swap eigenvalues i and j */
          stmp = ritz[j]; ritz[j] = ritz[i]; ritz[i] = stmp;
          rtmp = bnd[j]; bnd[j] = bnd[i]; bnd[i] = rtmp;
          ctmp = conv[j]; conv[j] = conv[i]; conv[i] = ctmp;
          /* swap eigenvectors i and j */
          PetscCall(BVGetColumn(eps->V,i,&vi));
          PetscCall(BVGetColumn(eps->V,j,&vj));
          PetscCall(VecSwap(vi,vj));
          PetscCall(BVRestoreColumn(eps->V,i,&vi));
          PetscCall(BVRestoreColumn(eps->V,j,&vj));
        }
      }
      k = i;
    }

    /* store ritz values and estimated errors */
    for (i=nconv;i<n;i++) {
      eps->eigr[i] = ritz[i];
      eps->errest[i] = bnd[i];
    }
    nconv = k;
    PetscCall(EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,n));
    PetscCall((*eps->stopping)(eps,eps->its,eps->max_it,nconv,eps->nev,&eps->reason,eps->stoppingctx));

    if (eps->reason == EPS_CONVERGED_ITERATING) { /* copy restart vector */
      PetscCall(BVCopyColumn(eps->V,n,nconv));
      if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL && !breakdown) {
        /* Reorthonormalize restart vector */
        PetscCall(BVOrthonormalizeColumn(eps->V,nconv,PETSC_FALSE,NULL,&breakdown));
      }
      if (breakdown) {
        /* Use random vector for restarting */
        PetscCall(PetscInfo(eps,"Using random vector for restart\n"));
        PetscCall(EPSGetStartVector(eps,nconv,&breakdown));
      }
      if (PetscUnlikely(breakdown)) { /* give up */
        eps->reason = EPS_DIVERGED_BREAKDOWN;
        PetscCall(PetscInfo(eps,"Unable to generate more start vectors\n"));
      }
    }
  }
  eps->nconv = nconv;

  PetscCall(PetscFree4(ritz,bnd,perm,conv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSetFromOptions_Lanczos(EPS eps,PetscOptionItems *PetscOptionsObject)
{
  EPS_LANCZOS            *lanczos = (EPS_LANCZOS*)eps->data;
  PetscBool              flg;
  EPSLanczosReorthogType reorthog=EPS_LANCZOS_REORTHOG_LOCAL,curval;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"EPS Lanczos Options");

    curval = (lanczos->reorthog==(EPSLanczosReorthogType)-1)? EPS_LANCZOS_REORTHOG_LOCAL: lanczos->reorthog;
    PetscCall(PetscOptionsEnum("-eps_lanczos_reorthog","Lanczos reorthogonalization","EPSLanczosSetReorthog",EPSLanczosReorthogTypes,(PetscEnum)curval,(PetscEnum*)&reorthog,&flg));
    if (flg) PetscCall(EPSLanczosSetReorthog(eps,reorthog));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLanczosSetReorthog_Lanczos(EPS eps,EPSLanczosReorthogType reorthog)
{
  EPS_LANCZOS *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  switch (reorthog) {
    case EPS_LANCZOS_REORTHOG_LOCAL:
    case EPS_LANCZOS_REORTHOG_FULL:
    case EPS_LANCZOS_REORTHOG_DELAYED:
    case EPS_LANCZOS_REORTHOG_SELECTIVE:
    case EPS_LANCZOS_REORTHOG_PERIODIC:
    case EPS_LANCZOS_REORTHOG_PARTIAL:
      if (lanczos->reorthog != reorthog) {
        lanczos->reorthog = reorthog;
        eps->state = EPS_STATE_INITIAL;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid reorthogonalization type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLanczosSetReorthog - Sets the type of reorthogonalization used during the Lanczos
   iteration.

   Logically Collective

   Input Parameters:
+  eps - the eigenproblem solver context
-  reorthog - the type of reorthogonalization

   Options Database Key:
.  -eps_lanczos_reorthog - Sets the reorthogonalization type (either 'local', 'selective',
                         'periodic', 'partial', 'full' or 'delayed')

   Level: advanced

.seealso: EPSLanczosGetReorthog(), EPSLanczosReorthogType
@*/
PetscErrorCode EPSLanczosSetReorthog(EPS eps,EPSLanczosReorthogType reorthog)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,reorthog,2);
  PetscTryMethod(eps,"EPSLanczosSetReorthog_C",(EPS,EPSLanczosReorthogType),(eps,reorthog));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSLanczosGetReorthog_Lanczos(EPS eps,EPSLanczosReorthogType *reorthog)
{
  EPS_LANCZOS *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  *reorthog = lanczos->reorthog;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   EPSLanczosGetReorthog - Gets the type of reorthogonalization used during
   the Lanczos iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  reorthog - the type of reorthogonalization

   Level: advanced

.seealso: EPSLanczosSetReorthog(), EPSLanczosReorthogType
@*/
PetscErrorCode EPSLanczosGetReorthog(EPS eps,EPSLanczosReorthogType *reorthog)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscAssertPointer(reorthog,2);
  PetscUseMethod(eps,"EPSLanczosGetReorthog_C",(EPS,EPSLanczosReorthogType*),(eps,reorthog));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_Lanczos(EPS eps)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  PetscCall(BVDestroy(&lanczos->AV));
  lanczos->allocsize = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_Lanczos(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosSetReorthog_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosGetReorthog_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSView_Lanczos(EPS eps,PetscViewer viewer)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (lanczos->reorthog != (EPSLanczosReorthogType)-1) PetscCall(PetscViewerASCIIPrintf(viewer,"  %s reorthogonalization\n",EPSLanczosReorthogTypes[lanczos->reorthog]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Lanczos(EPS eps)
{
  EPS_LANCZOS    *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;
  ctx->reorthog = (EPSLanczosReorthogType)-1;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_Lanczos;
  eps->ops->setup          = EPSSetUp_Lanczos;
  eps->ops->setupsort      = EPSSetUpSort_Default;
  eps->ops->setfromoptions = EPSSetFromOptions_Lanczos;
  eps->ops->destroy        = EPSDestroy_Lanczos;
  eps->ops->reset          = EPSReset_Lanczos;
  eps->ops->view           = EPSView_Lanczos;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Hermitian;

  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosSetReorthog_C",EPSLanczosSetReorthog_Lanczos));
  PetscCall(PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosGetReorthog_C",EPSLanczosGetReorthog_Lanczos));
  PetscFunctionReturn(PETSC_SUCCESS);
}
