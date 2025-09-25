/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for NEP, simple version
*/

#include <slepc/private/nepimpl.h>
#include <slepcblaslapack.h>

#define NREF_MAXIT 10

typedef struct {
  VecScatter    *scatter_id,nst;
  Mat           *A;
  Vec           nv,vg,v,w;
  FN            *fn;
} NEPSimpNRefctx;

typedef struct {
  Mat          M1;
  Vec          M2,M3;
  PetscScalar  M4,m3;
} NEP_REFINE_MATSHELL;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  NEP_REFINE_MATSHELL *ctx;
  PetscScalar         t;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(VecDot(x,ctx->M3,&t));
  t *= ctx->m3/ctx->M4;
  PetscCall(MatMult(ctx->M1,x,y));
  PetscCall(VecAXPY(y,-t,ctx->M2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSimpleNRefSetUp(NEP nep,NEPSimpNRefctx **ctx_)
{
  PetscInt       i,si,j,n0,m0,nloc,*idx1,*idx2,ne;
  IS             is1,is2;
  NEPSimpNRefctx *ctx;
  Vec            v;
  PetscMPIInt    rank,size;
  MPI_Comm       child;

  PetscFunctionBegin;
  PetscCall(PetscCalloc1(1,ctx_));
  ctx = *ctx_;
  if (nep->npart==1) {
    ctx->A  = nep->A;
    ctx->fn = nep->f;
    nep->refinesubc = NULL;
    ctx->scatter_id = NULL;
  } else {
    PetscCall(PetscSubcommGetChild(nep->refinesubc,&child));
    PetscCall(PetscMalloc2(nep->nt,&ctx->A,nep->npart,&ctx->scatter_id));

    /* Duplicate matrices */
    for (i=0;i<nep->nt;i++) PetscCall(MatCreateRedundantMatrix(nep->A[i],0,child,MAT_INITIAL_MATRIX,&ctx->A[i]));
    PetscCall(MatCreateVecs(ctx->A[0],&ctx->v,NULL));

    /* Duplicate FNs */
    PetscCall(PetscMalloc1(nep->nt,&ctx->fn));
    for (i=0;i<nep->nt;i++) PetscCall(FNDuplicate(nep->f[i],child,&ctx->fn[i]));

    /* Create scatters for sending vectors to each subcommucator */
    PetscCall(BVGetColumn(nep->V,0,&v));
    PetscCall(VecGetOwnershipRange(v,&n0,&m0));
    PetscCall(BVRestoreColumn(nep->V,0,&v));
    PetscCall(VecGetLocalSize(ctx->v,&nloc));
    PetscCall(PetscMalloc2(m0-n0,&idx1,m0-n0,&idx2));
    PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)nep),nloc,PETSC_DECIDE,&ctx->vg));
    for (si=0;si<nep->npart;si++) {
      j = 0;
      for (i=n0;i<m0;i++) {
        idx1[j]   = i;
        idx2[j++] = i+nep->n*si;
      }
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)nep),(m0-n0),idx1,PETSC_COPY_VALUES,&is1));
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)nep),(m0-n0),idx2,PETSC_COPY_VALUES,&is2));
      PetscCall(BVGetColumn(nep->V,0,&v));
      PetscCall(VecScatterCreate(v,is1,ctx->vg,is2,&ctx->scatter_id[si]));
      PetscCall(BVRestoreColumn(nep->V,0,&v));
      PetscCall(ISDestroy(&is1));
      PetscCall(ISDestroy(&is2));
    }
    PetscCall(PetscFree2(idx1,idx2));
  }
  if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ctx->A[0]),&rank));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ctx->A[0]),&size));
    if (size>1) {
      if (nep->npart==1) PetscCall(BVGetColumn(nep->V,0,&v));
      else v = ctx->v;
      PetscCall(VecGetOwnershipRange(v,&n0,&m0));
      ne = (rank == size-1)?nep->n:0;
      PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)ctx->A[0]),ne,PETSC_DECIDE,&ctx->nv));
      PetscCall(PetscMalloc1(m0-n0,&idx1));
      for (i=n0;i<m0;i++) {
        idx1[i-n0] = i;
      }
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)ctx->A[0]),(m0-n0),idx1,PETSC_COPY_VALUES,&is1));
      PetscCall(VecScatterCreate(v,is1,ctx->nv,is1,&ctx->nst));
      if (nep->npart==1) PetscCall(BVRestoreColumn(nep->V,0,&v));
      PetscCall(PetscFree(idx1));
      PetscCall(ISDestroy(&is1));
    }
  }  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Gather Eigenpair idx from subcommunicator with color sc
*/
static PetscErrorCode NEPSimpleNRefGatherEigenpair(NEP nep,NEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx,PetscInt *fail)
{
  PetscMPIInt    nproc,p;
  MPI_Comm       comm=((PetscObject)nep)->comm;
  Vec            v;
  PetscScalar    *array;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm,&nproc));
  PetscCall(PetscMPIIntCast((nproc/nep->npart)*(sc+1)+PetscMin(sc+1,nproc%nep->npart)-1,&p));
  if (nep->npart>1) {
    /* Communicate convergence successful */
    PetscCallMPI(MPI_Bcast(fail,1,MPIU_INT,p,comm));
    if (!(*fail)) {
      /* Process 0 of subcommunicator sc broadcasts the eigenvalue */
      PetscCallMPI(MPI_Bcast(&nep->eigr[idx],1,MPIU_SCALAR,p,comm));
      /* Gather nep->V[idx] from the subcommuniator sc */
      PetscCall(BVGetColumn(nep->V,idx,&v));
      if (nep->refinesubc->color==sc) {
        PetscCall(VecGetArray(ctx->v,&array));
        PetscCall(VecPlaceArray(ctx->vg,array));
      }
      PetscCall(VecScatterBegin(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE));
      PetscCall(VecScatterEnd(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE));
      if (nep->refinesubc->color==sc) {
        PetscCall(VecResetArray(ctx->vg));
        PetscCall(VecRestoreArray(ctx->v,&array));
      }
      PetscCall(BVRestoreColumn(nep->V,idx,&v));
    }
  } else {
    if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT && !(*fail)) PetscCallMPI(MPI_Bcast(&nep->eigr[idx],1,MPIU_SCALAR,p,comm));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSimpleNRefScatterEigenvector(NEP nep,NEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx)
{
  Vec            v;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (nep->npart>1) {
    PetscCall(BVGetColumn(nep->V,idx,&v));
    if (nep->refinesubc->color==sc) {
      PetscCall(VecGetArray(ctx->v,&array));
      PetscCall(VecPlaceArray(ctx->vg,array));
    }
    PetscCall(VecScatterBegin(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD));
    if (nep->refinesubc->color==sc) {
      PetscCall(VecResetArray(ctx->vg));
      PetscCall(VecRestoreArray(ctx->v,&array));
    }
    PetscCall(BVRestoreColumn(nep->V,idx,&v));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPSimpleNRefSetUpSystem(NEP nep,NEPSimpNRefctx *ctx,Mat *A,PetscInt idx,Mat *Mt,Mat *T,Mat *P,PetscBool ini,Vec t,Vec v)
{
  PetscInt            i,st,ml,m0,n0,m1,mg;
  PetscInt            ncols,*cols2=NULL,nt=nep->nt;
  PetscScalar         zero=0.0,*coeffs,*coeffs2;
  PetscMPIInt         rank,size;
  MPI_Comm            comm;
  const PetscInt      *cols;
  const PetscScalar   *vals,*array;
  NEP_REFINE_MATSHELL *fctx;
  Vec                 w=ctx->w;
  Mat                 M;

  PetscFunctionBegin;
  PetscCall(PetscMalloc2(nt,&coeffs,nt,&coeffs2));
  switch (nep->scheme) {
  case NEP_REFINE_SCHEME_SCHUR:
    if (ini) {
      PetscCall(PetscCalloc1(1,&fctx));
      PetscCall(MatGetSize(A[0],&m0,&n0));
      PetscCall(MatCreateShell(PetscObjectComm((PetscObject)A[0]),PETSC_DECIDE,PETSC_DECIDE,m0,n0,fctx,T));
      PetscCall(MatShellSetOperation(*T,MATOP_MULT,(PetscErrorCodeFn*)MatMult_FS));
    } else PetscCall(MatShellGetContext(*T,&fctx));
    M=fctx->M1;
    break;
  case NEP_REFINE_SCHEME_MBE:
    M=*T;
    break;
  case NEP_REFINE_SCHEME_EXPLICIT:
    M=*Mt;
    break;
  }
  if (ini) PetscCall(MatDuplicate(A[0],MAT_COPY_VALUES,&M));
  else PetscCall(MatCopy(A[0],M,DIFFERENT_NONZERO_PATTERN));
  for (i=0;i<nt;i++) PetscCall(FNEvaluateFunction(ctx->fn[i],nep->eigr[idx],coeffs+i));
  if (coeffs[0]!=1.0) PetscCall(MatScale(M,coeffs[0]));
  for (i=1;i<nt;i++) PetscCall(MatAXPY(M,coeffs[i],A[i],(ini)?nep->mstr:SUBSET_NONZERO_PATTERN));
  for (i=0;i<nt;i++) PetscCall(FNEvaluateDerivative(ctx->fn[i],nep->eigr[idx],coeffs2+i));
  st = 0;
  for (i=0;i<nt && PetscAbsScalar(coeffs2[i])==0.0;i++) st++;
  PetscCall(MatMult(A[st],v,w));
  if (coeffs2[st]!=1.0) PetscCall(VecScale(w,coeffs2[st]));
  for (i=st+1;i<nt;i++) {
    PetscCall(MatMult(A[i],v,t));
    PetscCall(VecAXPY(w,coeffs2[i],t));
  }

  switch (nep->scheme) {
  case NEP_REFINE_SCHEME_EXPLICIT:
    comm = PetscObjectComm((PetscObject)A[0]);
    PetscCallMPI(MPI_Comm_rank(comm,&rank));
    PetscCallMPI(MPI_Comm_size(comm,&size));
    PetscCall(MatGetSize(M,&mg,NULL));
    PetscCall(MatGetOwnershipRange(M,&m0,&m1));
    if (ini) {
      PetscCall(MatCreate(comm,T));
      PetscCall(MatGetLocalSize(M,&ml,NULL));
      if (rank==size-1) ml++;
      PetscCall(MatSetSizes(*T,ml,ml,mg+1,mg+1));
      PetscCall(MatSetFromOptions(*T));
      *Mt = M;
      *P  = *T;
    }

    /* Set values */
    PetscCall(VecGetArrayRead(w,&array));
    for (i=m0;i<m1;i++) {
      PetscCall(MatGetRow(M,i,&ncols,&cols,&vals));
      PetscCall(MatSetValues(*T,1,&i,ncols,cols,vals,INSERT_VALUES));
      PetscCall(MatRestoreRow(M,i,&ncols,&cols,&vals));
      PetscCall(MatSetValues(*T,1,&i,1,&mg,array+i-m0,INSERT_VALUES));
    }
    PetscCall(VecRestoreArrayRead(w,&array));
    PetscCall(VecConjugate(v));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)A[0]),&size));
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)A[0]),&rank));
    if (size>1) {
      if (rank==size-1) {
        PetscCall(PetscMalloc1(nep->n,&cols2));
        for (i=0;i<nep->n;i++) cols2[i]=i;
      }
      PetscCall(VecScatterBegin(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecGetArrayRead(ctx->nv,&array));
      if (rank==size-1) {
        PetscCall(MatSetValues(*T,1,&mg,nep->n,cols2,array,INSERT_VALUES));
        PetscCall(MatSetValues(*T,1,&mg,1,&mg,&zero,INSERT_VALUES));
      }
      PetscCall(VecRestoreArrayRead(ctx->nv,&array));
    } else {
      PetscCall(PetscMalloc1(m1-m0,&cols2));
      for (i=0;i<m1-m0;i++) cols2[i]=m0+i;
      PetscCall(VecGetArrayRead(v,&array));
      PetscCall(MatSetValues(*T,1,&mg,m1-m0,cols2,array,INSERT_VALUES));
      PetscCall(MatSetValues(*T,1,&mg,1,&mg,&zero,INSERT_VALUES));
      PetscCall(VecRestoreArrayRead(v,&array));
    }
    PetscCall(VecConjugate(v));
    PetscCall(MatAssemblyBegin(*T,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*T,MAT_FINAL_ASSEMBLY));
    PetscCall(PetscFree(cols2));
    break;
  case NEP_REFINE_SCHEME_SCHUR:
    fctx->M2 = ctx->w;
    fctx->M3 = v;
    fctx->m3 = 1.0+PetscConj(nep->eigr[idx])*nep->eigr[idx];
    fctx->M4 = PetscConj(nep->eigr[idx]);
    fctx->M1 = M;
    if (ini) PetscCall(MatDuplicate(M,MAT_COPY_VALUES,P));
    else PetscCall(MatCopy(M,*P,SAME_NONZERO_PATTERN));
    if (fctx->M4!=0.0) {
      PetscCall(VecConjugate(v));
      PetscCall(VecPointwiseMult(t,v,w));
      PetscCall(VecConjugate(v));
      PetscCall(VecScale(t,-fctx->m3/fctx->M4));
      PetscCall(MatDiagonalSet(*P,t,ADD_VALUES));
    }
    break;
  case NEP_REFINE_SCHEME_MBE:
    *T = M;
    *P = M;
    break;
  }
  PetscCall(PetscFree2(coeffs,coeffs2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPNewtonRefinementSimple(NEP nep,PetscInt *maxits,PetscReal tol,PetscInt k)
{
  PetscInt            i,n,its,idx=0,*idx_sc,*its_sc,color,*fail_sc;
  PetscMPIInt         rank,size;
  Mat                 Mt=NULL,T=NULL,P=NULL;
  MPI_Comm            comm;
  Vec                 r,v,dv,rr=NULL,dvv=NULL,t[2];
  const PetscScalar   *array;
  PetscScalar         *array2,deig=0.0,tt[2],ttt;
  PetscReal           norm,error;
  PetscBool           ini=PETSC_TRUE,sc_pend,solved=PETSC_FALSE;
  NEPSimpNRefctx      *ctx;
  NEP_REFINE_MATSHELL *fctx=NULL;
  KSPConvergedReason  reason;

  PetscFunctionBegin;
  PetscCall(PetscLogEventBegin(NEP_Refine,nep,0,0,0));
  PetscCall(NEPSimpleNRefSetUp(nep,&ctx));
  its = (maxits)?*maxits:NREF_MAXIT;
  if (!nep->refineksp) PetscCall(NEPRefineGetKSP(nep,&nep->refineksp));
  if (nep->npart==1) PetscCall(BVGetColumn(nep->V,0,&v));
  else v = ctx->v;
  PetscCall(VecDuplicate(v,&ctx->w));
  PetscCall(VecDuplicate(v,&r));
  PetscCall(VecDuplicate(v,&dv));
  PetscCall(VecDuplicate(v,&t[0]));
  PetscCall(VecDuplicate(v,&t[1]));
  if (nep->npart==1) {
    PetscCall(BVRestoreColumn(nep->V,0,&v));
    PetscCall(PetscObjectGetComm((PetscObject)nep,&comm));
  } else PetscCall(PetscSubcommGetChild(nep->refinesubc,&comm));
  PetscCallMPI(MPI_Comm_size(comm,&size));
  PetscCallMPI(MPI_Comm_rank(comm,&rank));
  PetscCall(VecGetLocalSize(r,&n));
  PetscCall(PetscMalloc3(nep->npart,&idx_sc,nep->npart,&its_sc,nep->npart,&fail_sc));
  for (i=0;i<nep->npart;i++) fail_sc[i] = 0;
  for (i=0;i<nep->npart;i++) its_sc[i] = 0;
  color = (nep->npart==1)?0:nep->refinesubc->color;

  /* Loop performing iterative refinements */
  while (!solved) {
    for (i=0;i<nep->npart;i++) {
      sc_pend = PETSC_TRUE;
      if (its_sc[i]==0) {
        idx_sc[i] = idx++;
        if (idx_sc[i]>=k) {
          sc_pend = PETSC_FALSE;
        } else PetscCall(NEPSimpleNRefScatterEigenvector(nep,ctx,i,idx_sc[i]));
      }  else { /* Gather Eigenpair from subcommunicator i */
        PetscCall(NEPSimpleNRefGatherEigenpair(nep,ctx,i,idx_sc[i],&fail_sc[i]));
      }
      while (sc_pend) {
        if (!fail_sc[i]) PetscCall(NEPComputeError(nep,idx_sc[i],NEP_ERROR_RELATIVE,&error));
        if (error<=tol || its_sc[i]>=its || fail_sc[i]) {
          idx_sc[i] = idx++;
          its_sc[i] = 0;
          fail_sc[i] = 0;
          if (idx_sc[i]<k) PetscCall(NEPSimpleNRefScatterEigenvector(nep,ctx,i,idx_sc[i]));
        } else {
          sc_pend = PETSC_FALSE;
          its_sc[i]++;
        }
        if (idx_sc[i]>=k) sc_pend = PETSC_FALSE;
      }
    }
    solved = PETSC_TRUE;
    for (i=0;i<nep->npart&&solved;i++) solved = PetscNot(idx_sc[i]<k);
    if (idx_sc[color]<k) {
#if !defined(PETSC_USE_COMPLEX)
      PetscCheck(nep->eigi[idx_sc[color]]==0.0,PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Simple Refinement not implemented in real scalar for complex eigenvalues");
#endif
      if (nep->npart==1) PetscCall(BVGetColumn(nep->V,idx_sc[color],&v));
      else v = ctx->v;
      PetscCall(NEPSimpleNRefSetUpSystem(nep,ctx,ctx->A,idx_sc[color],&Mt,&T,&P,ini,t[0],v));
      PetscCall(NEP_KSPSetOperators(nep->refineksp,T,P));
      if (ini) {
        PetscCall(KSPSetFromOptions(nep->refineksp));
        if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
          PetscCall(MatCreateVecs(T,&dvv,NULL));
          PetscCall(VecDuplicate(dvv,&rr));
        }
        ini = PETSC_FALSE;
      }
      switch (nep->scheme) {
      case NEP_REFINE_SCHEME_EXPLICIT:
        PetscCall(MatMult(Mt,v,r));
        PetscCall(VecGetArrayRead(r,&array));
        if (rank==size-1) {
          PetscCall(VecGetArray(rr,&array2));
          PetscCall(PetscArraycpy(array2,array,n));
          array2[n] = 0.0;
          PetscCall(VecRestoreArray(rr,&array2));
        } else PetscCall(VecPlaceArray(rr,array));
        PetscCall(KSPSolve(nep->refineksp,rr,dvv));
        PetscCall(KSPGetConvergedReason(nep->refineksp,&reason));
        if (reason>0) {
          if (rank != size-1) PetscCall(VecResetArray(rr));
          PetscCall(VecRestoreArrayRead(r,&array));
          PetscCall(VecGetArrayRead(dvv,&array));
          PetscCall(VecPlaceArray(dv,array));
          PetscCall(VecAXPY(v,-1.0,dv));
          PetscCall(VecNorm(v,NORM_2,&norm));
          PetscCall(VecScale(v,1.0/norm));
          PetscCall(VecResetArray(dv));
          if (rank==size-1) nep->eigr[idx_sc[color]] -= array[n];
          PetscCall(VecRestoreArrayRead(dvv,&array));
        } else fail_sc[color] = 1;
        break;
      case NEP_REFINE_SCHEME_MBE:
        PetscCall(MatMult(T,v,r));
        /* Mixed block elimination */
        PetscCall(VecConjugate(v));
        PetscCall(KSPSolveTranspose(nep->refineksp,v,t[0]));
        PetscCall(KSPGetConvergedReason(nep->refineksp,&reason));
        if (reason>0) {
          PetscCall(VecConjugate(t[0]));
          PetscCall(VecDot(ctx->w,t[0],&tt[0]));
          PetscCall(KSPSolve(nep->refineksp,ctx->w,t[1]));
          PetscCall(KSPGetConvergedReason(nep->refineksp,&reason));
          if (reason>0) {
            PetscCall(VecDot(t[1],v,&tt[1]));
            PetscCall(VecDot(r,t[0],&ttt));
            tt[0] = ttt/tt[0];
            PetscCall(VecAXPY(r,-tt[0],ctx->w));
            PetscCall(KSPSolve(nep->refineksp,r,dv));
            PetscCall(KSPGetConvergedReason(nep->refineksp,&reason));
            if (reason>0) {
              PetscCall(VecDot(dv,v,&ttt));
              tt[1] = ttt/tt[1];
              PetscCall(VecAXPY(dv,-tt[1],t[1]));
              deig = tt[0]+tt[1];
            }
          }
          PetscCall(VecConjugate(v));
          PetscCall(VecAXPY(v,-1.0,dv));
          PetscCall(VecNorm(v,NORM_2,&norm));
          PetscCall(VecScale(v,1.0/norm));
          nep->eigr[idx_sc[color]] -= deig;
          fail_sc[color] = 0;
        } else {
          PetscCall(VecConjugate(v));
          fail_sc[color] = 1;
        }
        break;
      case NEP_REFINE_SCHEME_SCHUR:
        fail_sc[color] = 1;
        PetscCall(MatShellGetContext(T,&fctx));
        if (fctx->M4!=0.0) {
          PetscCall(MatMult(fctx->M1,v,r));
          PetscCall(KSPSolve(nep->refineksp,r,dv));
          PetscCall(KSPGetConvergedReason(nep->refineksp,&reason));
          if (reason>0) {
            PetscCall(VecDot(dv,v,&deig));
            deig *= -fctx->m3/fctx->M4;
            PetscCall(VecAXPY(v,-1.0,dv));
            PetscCall(VecNorm(v,NORM_2,&norm));
            PetscCall(VecScale(v,1.0/norm));
            nep->eigr[idx_sc[color]] -= deig;
            fail_sc[color] = 0;
          }
        }
        break;
      }
      if (nep->npart==1) PetscCall(BVRestoreColumn(nep->V,idx_sc[color],&v));
    }
  }
  PetscCall(VecDestroy(&t[0]));
  PetscCall(VecDestroy(&t[1]));
  PetscCall(VecDestroy(&dv));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(VecDestroy(&r));
  PetscCall(PetscFree3(idx_sc,its_sc,fail_sc));
  PetscCall(VecScatterDestroy(&ctx->nst));
  if (nep->npart>1) {
    PetscCall(VecDestroy(&ctx->vg));
    PetscCall(VecDestroy(&ctx->v));
    for (i=0;i<nep->nt;i++) PetscCall(MatDestroy(&ctx->A[i]));
    for (i=0;i<nep->npart;i++) PetscCall(VecScatterDestroy(&ctx->scatter_id[i]));
    PetscCall(PetscFree2(ctx->A,ctx->scatter_id));
  }
  if (fctx && nep->scheme==NEP_REFINE_SCHEME_SCHUR) {
    PetscCall(MatDestroy(&P));
    PetscCall(MatDestroy(&fctx->M1));
    PetscCall(PetscFree(fctx));
  }
  if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
    PetscCall(MatDestroy(&Mt));
    PetscCall(VecDestroy(&dvv));
    PetscCall(VecDestroy(&rr));
    PetscCall(VecDestroy(&ctx->nv));
    if (nep->npart>1) {
      for (i=0;i<nep->nt;i++) PetscCall(FNDestroy(&ctx->fn[i]));
      PetscCall(PetscFree(ctx->fn));
    }
  }
  if (nep->scheme==NEP_REFINE_SCHEME_MBE) {
    if (nep->npart>1) {
      for (i=0;i<nep->nt;i++) PetscCall(FNDestroy(&ctx->fn[i]));
      PetscCall(PetscFree(ctx->fn));
    }
  }
  PetscCall(MatDestroy(&T));
  PetscCall(PetscFree(ctx));
  PetscCall(PetscLogEventEnd(NEP_Refine,nep,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
