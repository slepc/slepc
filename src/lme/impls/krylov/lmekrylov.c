/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix equation solver: "krylov"

   Method: Arnoldi with Eiermann-Ernst restart

   Algorithm:

       Project the equation onto the Arnoldi basis and solve the compressed
       equation the Hessenberg matrix H, restart by discarding the Krylov
       basis but keeping H.

   References:

       [1] Y. Saad, "Numerical solution of large Lyapunov equations", in
           Signal processing, scattering and operator theory, and numerical
           methods, vol. 5 of Progr. Systems Control Theory, pages 503-511,
           1990.

       [2] D. Kressner, "Memory-efficient Krylov subspace techniques for
           solving large-scale Lyapunov equations", in 2008 IEEE Int. Conf.
           Computer-Aided Control Systems, pages 613-618, 2008.
*/

#include <slepc/private/lmeimpl.h>
#include <slepcblaslapack.h>

static PetscErrorCode LMESetUp_Krylov(LME lme)
{
  PetscInt       N;

  PetscFunctionBegin;
  PetscCall(MatGetSize(lme->A,&N,NULL));
  if (lme->ncv==PETSC_DETERMINE) lme->ncv = PetscMin(30,N);
  if (lme->max_it==PETSC_DETERMINE) lme->max_it = 100;
  PetscCall(LMEAllocateSolution(lme,1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode LMESolve_Krylov_Lyapunov_Vec(LME lme,Vec b,PetscBool fixed,PetscInt rrank,BV C1,BV *X1,PetscInt *col,PetscBool *fail,PetscInt *totalits)
{
  PetscInt       n=0,m,ldh,ldg=0,i,j,rank=0,lrank,pass,nouter=0,its;
  PetscReal      bnorm,beta,errest;
  PetscBool      breakdown;
  PetscScalar    *Harray,*G=NULL,*Gnew=NULL,*L=NULL,*U=NULL,*r,*Qarray,sone=1.0,zero=0.0;
  PetscBLASInt   n_,m_,rk_;
  Mat            Q,H;

  PetscFunctionBegin;
  *fail = PETSC_FALSE;
  its = 0;
  m  = lme->ncv;
  ldh = m+1;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ldh,m,NULL,&H));
  PetscCall(MatDenseGetArray(H,&Harray));

  PetscCall(VecNorm(b,NORM_2,&bnorm));
  PetscCheck(bnorm,PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Cannot process a zero vector in the right-hand side");

  for (pass=0;pass<2;pass++) {

    /* set initial vector to b/||b|| */
    PetscCall(BVInsertVec(lme->V,0,b));
    PetscCall(BVScaleColumn(lme->V,0,1.0/bnorm));

    /* Restart loop */
    while ((pass==0 && !*fail) || (pass==1 && its+1<nouter)) {
      its++;

      /* compute Arnoldi factorization */
      PetscCall(BVMatArnoldi(lme->V,lme->A,H,0,&m,&beta,&breakdown));
      PetscCall(BVSetActiveColumns(lme->V,0,m));

      if (pass==0) {
        /* glue together the previous H and the new H obtained with Arnoldi */
        ldg = n+m+1;
        PetscCall(PetscCalloc1(ldg*(n+m),&Gnew));
        for (j=0;j<m;j++) PetscCall(PetscArraycpy(Gnew+n+(j+n)*ldg,Harray+j*ldh,m));
        Gnew[n+m+(n+m-1)*ldg] = beta;
        if (G) {
          for (j=0;j<n;j++) PetscCall(PetscArraycpy(Gnew+j*ldg,G+j*(n+1),n+1));
          PetscCall(PetscFree(G));
        }
        G = Gnew;
        n += m;
      } else {
        /* update Z = Z + V(:,1:m)*Q    with   Q=U(blk,:)*P(1:nrk,:)'  */
        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,m,*col+rank,NULL,&Q));
        PetscCall(MatDenseGetArray(Q,&Qarray));
        PetscCall(PetscBLASIntCast(m,&m_));
        PetscCall(PetscBLASIntCast(n,&n_));
        PetscCall(PetscBLASIntCast(rank,&rk_));
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&m_,&rk_,&rk_,&sone,U+its*m,&n_,L,&n_,&zero,Qarray+(*col)*m,&m_));
        PetscCall(MatDenseRestoreArray(Q,&Qarray));
        PetscCall(BVSetActiveColumns(*X1,*col,*col+rank));
        PetscCall(BVMult(*X1,1.0,1.0,lme->V,Q));
        PetscCall(MatDestroy(&Q));
      }

      if (pass==0) {
        /* solve compressed Lyapunov equation */
        PetscCall(PetscCalloc1(n,&r));
        PetscCall(PetscCalloc1(n*n,&L));
        r[0] = bnorm;
        errest = PetscAbsScalar(G[n+(n-1)*ldg]);
        PetscCall(LMEDenseHessLyapunovChol(lme,n,G,ldg,1,r,n,L,n,&errest));
        PetscCall(LMEMonitor(lme,*totalits+its,errest));
        PetscCall(PetscFree(r));

        /* check convergence */
        if (errest<lme->tol) {
          lme->errest += errest;
          PetscCall(PetscMalloc1(n*n,&U));
          /* transpose L */
          for (j=0;j<n;j++) {
            for (i=j+1;i<n;i++) {
              L[i+j*n] = PetscConj(L[j+i*n]);
              L[j+i*n] = 0.0;
            }
          }
          PetscCall(LMEDenseRankSVD(lme,n,L,n,U,n,&lrank));
          PetscCall(PetscInfo(lme,"Rank of the Cholesky factor = %" PetscInt_FMT "\n",lrank));
          nouter = its;
          its = -1;
          if (!fixed) {  /* X1 was not set by user, allocate it with rank columns */
            rank = lrank;
            if (*col) PetscCall(BVResize(*X1,*col+rank,PETSC_TRUE));
            else PetscCall(BVDuplicateResize(C1,rank,X1));
          } else rank = PetscMin(lrank,rrank);
          PetscCall(PetscFree(G));
          break;
        } else {
          PetscCall(PetscFree(L));
          if (*totalits+its>=lme->max_it) *fail = PETSC_TRUE;
        }
      }

      /* restart with vector v_{m+1} */
      if (!*fail) PetscCall(BVCopyColumn(lme->V,m,0));
    }
  }

  *col += rank;
  *totalits += its+1;
  PetscCall(MatDenseRestoreArray(H,&Harray));
  PetscCall(MatDestroy(&H));
  if (L) PetscCall(PetscFree(L));
  if (U) PetscCall(PetscFree(U));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode LMESolve_Krylov_Lyapunov(LME lme)
{
  PetscBool      fail,fixed = lme->X? PETSC_TRUE: PETSC_FALSE;
  PetscInt       i,k,rank=0,col=0;
  Vec            b;
  BV             X1=NULL,C1;
  Mat            X1m,X1t,C1m;

  PetscFunctionBegin;
  PetscCall(MatLRCGetMats(lme->C,NULL,&C1m,NULL,NULL));
  PetscCall(BVCreateFromMat(C1m,&C1));
  PetscCall(BVSetFromOptions(C1));
  PetscCall(BVGetActiveColumns(C1,NULL,&k));
  if (fixed) {
    PetscCall(MatLRCGetMats(lme->X,NULL,&X1m,NULL,NULL));
    PetscCall(BVCreateFromMat(X1m,&X1));
    PetscCall(BVSetFromOptions(X1));
    PetscCall(BVGetActiveColumns(X1,NULL,&rank));
    rank = rank/k;
  }
  for (i=0;i<k;i++) {
    PetscCall(BVGetColumn(C1,i,&b));
    PetscCall(LMESolve_Krylov_Lyapunov_Vec(lme,b,fixed,rank,C1,&X1,&col,&fail,&lme->its));
    PetscCall(BVRestoreColumn(C1,i,&b));
    if (fail) {
      lme->reason = LME_DIVERGED_ITS;
      break;
    }
  }
  if (lme->reason==LME_CONVERGED_ITERATING) lme->reason = LME_CONVERGED_TOL;
  PetscCall(BVCreateMat(X1,&X1t));
  if (fixed) PetscCall(MatCopy(X1t,X1m,SAME_NONZERO_PATTERN));
  else PetscCall(MatCreateLRC(NULL,X1t,NULL,NULL,&lme->X));
  PetscCall(MatDestroy(&X1t));
  PetscCall(BVDestroy(&C1));
  PetscCall(BVDestroy(&X1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   LMEKRYLOV - LMEKRYLOV = "krylov" - A restarted Krylov method for the solution of
   linear matrix equations.

   Note:
   A description of the method is provided in chapter [](ch:lme)
   of the Users Manual.

   Level: beginner

.seealso: [](ch:lme), `LME`, `LMEType`, `LMESetType()`
M*/
SLEPC_EXTERN PetscErrorCode LMECreate_Krylov(LME lme)
{
  PetscFunctionBegin;
  lme->ops->solve[LME_LYAPUNOV]      = LMESolve_Krylov_Lyapunov;
  lme->ops->setup                    = LMESetUp_Krylov;
  PetscFunctionReturn(PETSC_SUCCESS);
}
