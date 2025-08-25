/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for PEP, simple version
*/

#include <slepc/private/pepimpl.h>
#include <slepcblaslapack.h>

#define NREF_MAXIT 10

typedef struct {
  VecScatter *scatter_id,nst;
  Mat        *A;
  Vec        nv,vg,v,w;
} PEPSimpNRefctx;

typedef struct {
  Mat          M1;
  Vec          M2,M3;
  PetscScalar  M4,m3;
} PEP_REFINES_MATSHELL;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  PEP_REFINES_MATSHELL *ctx;
  PetscScalar          t;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&ctx));
  PetscCall(VecDot(x,ctx->M3,&t));
  t *= ctx->m3/ctx->M4;
  PetscCall(MatMult(ctx->M1,x,y));
  PetscCall(VecAXPY(y,-t,ctx->M2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSimpleNRefSetUp(PEP pep,PEPSimpNRefctx **ctx_)
{
  PetscInt       i,si,j,n0,m0,nloc,*idx1,*idx2,ne;
  IS             is1,is2;
  PEPSimpNRefctx *ctx;
  Vec            v;
  PetscMPIInt    rank,size;
  MPI_Comm       child;

  PetscFunctionBegin;
  PetscCall(PetscCalloc1(1,ctx_));
  ctx = *ctx_;
  if (pep->npart==1) {
    pep->refinesubc = NULL;
    ctx->scatter_id = NULL;
    ctx->A = pep->A;
  } else {
    PetscCall(PetscSubcommGetChild(pep->refinesubc,&child));
    PetscCall(PetscMalloc2(pep->nmat,&ctx->A,pep->npart,&ctx->scatter_id));

    /* Duplicate matrices */
    for (i=0;i<pep->nmat;i++) PetscCall(MatCreateRedundantMatrix(pep->A[i],0,child,MAT_INITIAL_MATRIX,&ctx->A[i]));
    PetscCall(MatCreateVecs(ctx->A[0],&ctx->v,NULL));

    /* Create scatters for sending vectors to each subcommucator */
    PetscCall(BVGetColumn(pep->V,0,&v));
    PetscCall(VecGetOwnershipRange(v,&n0,&m0));
    PetscCall(BVRestoreColumn(pep->V,0,&v));
    PetscCall(VecGetLocalSize(ctx->v,&nloc));
    PetscCall(PetscMalloc2(m0-n0,&idx1,m0-n0,&idx2));
    PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)pep),nloc,PETSC_DECIDE,&ctx->vg));
    for (si=0;si<pep->npart;si++) {
      j = 0;
      for (i=n0;i<m0;i++) {
        idx1[j]   = i;
        idx2[j++] = i+pep->n*si;
      }
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),(m0-n0),idx1,PETSC_COPY_VALUES,&is1));
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pep),(m0-n0),idx2,PETSC_COPY_VALUES,&is2));
      PetscCall(BVGetColumn(pep->V,0,&v));
      PetscCall(VecScatterCreate(v,is1,ctx->vg,is2,&ctx->scatter_id[si]));
      PetscCall(BVRestoreColumn(pep->V,0,&v));
      PetscCall(ISDestroy(&is1));
      PetscCall(ISDestroy(&is2));
    }
    PetscCall(PetscFree2(idx1,idx2));
  }
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)ctx->A[0]),&rank));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ctx->A[0]),&size));
    if (size>1) {
      if (pep->npart==1) PetscCall(BVGetColumn(pep->V,0,&v));
      else v = ctx->v;
      PetscCall(VecGetOwnershipRange(v,&n0,&m0));
      ne = (rank == size-1)?pep->n:0;
      PetscCall(VecCreateMPI(PetscObjectComm((PetscObject)ctx->A[0]),ne,PETSC_DECIDE,&ctx->nv));
      PetscCall(PetscMalloc1(m0-n0,&idx1));
      for (i=n0;i<m0;i++) idx1[i-n0] = i;
      PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)ctx->A[0]),(m0-n0),idx1,PETSC_COPY_VALUES,&is1));
      PetscCall(VecScatterCreate(v,is1,ctx->nv,is1,&ctx->nst));
      if (pep->npart==1) PetscCall(BVRestoreColumn(pep->V,0,&v));
      PetscCall(PetscFree(idx1));
      PetscCall(ISDestroy(&is1));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Gather Eigenpair idx from subcommunicator with color sc
*/
static PetscErrorCode PEPSimpleNRefGatherEigenpair(PEP pep,PEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx,PetscInt *fail)
{
  PetscMPIInt       nproc,p;
  MPI_Comm          comm=((PetscObject)pep)->comm;
  Vec               v;
  const PetscScalar *array;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm,&nproc));
  PetscCall(PetscMPIIntCast((nproc/pep->npart)*(sc+1)+PetscMin(sc+1,nproc%pep->npart)-1,&p));
  if (pep->npart>1) {
    /* Communicate convergence successful */
    PetscCallMPI(MPI_Bcast(fail,1,MPIU_INT,p,comm));
    if (!(*fail)) {
      /* Process 0 of subcommunicator sc broadcasts the eigenvalue */
      PetscCallMPI(MPI_Bcast(&pep->eigr[idx],1,MPIU_SCALAR,p,comm));
      /* Gather pep->V[idx] from the subcommuniator sc */
      PetscCall(BVGetColumn(pep->V,idx,&v));
      if (pep->refinesubc->color==sc) {
        PetscCall(VecGetArrayRead(ctx->v,&array));
        PetscCall(VecPlaceArray(ctx->vg,array));
      }
      PetscCall(VecScatterBegin(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE));
      PetscCall(VecScatterEnd(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE));
      if (pep->refinesubc->color==sc) {
        PetscCall(VecResetArray(ctx->vg));
        PetscCall(VecRestoreArrayRead(ctx->v,&array));
      }
      PetscCall(BVRestoreColumn(pep->V,idx,&v));
    }
  } else {
    if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT && !(*fail)) PetscCallMPI(MPI_Bcast(&pep->eigr[idx],1,MPIU_SCALAR,p,comm));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSimpleNRefScatterEigenvector(PEP pep,PEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx)
{
  Vec               v;
  const PetscScalar *array;

  PetscFunctionBegin;
  if (pep->npart>1) {
    PetscCall(BVGetColumn(pep->V,idx,&v));
    if (pep->refinesubc->color==sc) {
      PetscCall(VecGetArrayRead(ctx->v,&array));
      PetscCall(VecPlaceArray(ctx->vg,array));
    }
    PetscCall(VecScatterBegin(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD));
    if (pep->refinesubc->color==sc) {
      PetscCall(VecResetArray(ctx->vg));
      PetscCall(VecRestoreArrayRead(ctx->v,&array));
    }
    PetscCall(BVRestoreColumn(pep->V,idx,&v));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPEvaluateFunctionDerivatives(PEP pep,PetscScalar alpha,PetscScalar *vals)
{
  PetscInt    i,nmat=pep->nmat;
  PetscScalar a0,a1,a2;
  PetscReal   *a=pep->pbc,*b=a+nmat,*g=b+nmat;

  PetscFunctionBegin;
  a0 = 0.0;
  a1 = 1.0;
  vals[0] = 0.0;
  if (nmat>1) vals[1] = 1/a[0];
  for (i=2;i<nmat;i++) {
    a2 = ((alpha-b[i-2])*a1-g[i-2]*a0)/a[i-2];
    vals[i] = (a2+(alpha-b[i-1])*vals[i-1]-g[i-1]*vals[i-2])/a[i-1];
    a0 = a1; a1 = a2;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSimpleNRefSetUpSystem(PEP pep,Mat *A,PEPSimpNRefctx *ctx,PetscInt idx,Mat *Mt,Mat *T,Mat *P,PetscBool ini,Vec t,Vec v)
{
  PetscInt             i,nmat=pep->nmat,ml,m0,n0,m1,mg;
  PetscInt             ncols,*cols2=NULL;
  PetscScalar          zero=0.0,*coeffs,*coeffs2;
  PetscMPIInt          rank,size;
  MPI_Comm             comm;
  const PetscInt       *cols;
  const PetscScalar    *vals,*array;
  MatStructure         str;
  PEP_REFINES_MATSHELL *fctx;
  PEPRefineScheme      scheme=pep->scheme;
  Vec                  w=ctx->w;
  Mat                  M;

  PetscFunctionBegin;
  PetscCall(STGetMatStructure(pep->st,&str));
  PetscCall(PetscMalloc2(nmat,&coeffs,nmat,&coeffs2));
  switch (scheme) {
  case PEP_REFINE_SCHEME_SCHUR:
    if (ini) {
      PetscCall(PetscCalloc1(1,&fctx));
      PetscCall(MatGetSize(A[0],&m0,&n0));
      PetscCall(MatCreateShell(PetscObjectComm((PetscObject)A[0]),PETSC_DECIDE,PETSC_DECIDE,m0,n0,fctx,T));
      PetscCall(MatShellSetOperation(*T,MATOP_MULT,(PetscErrorCodeFn*)MatMult_FS));
    } else PetscCall(MatShellGetContext(*T,&fctx));
    M=fctx->M1;
    break;
  case PEP_REFINE_SCHEME_MBE:
    M=*T;
    break;
  case PEP_REFINE_SCHEME_EXPLICIT:
    M=*Mt;
    break;
  }
  if (ini) PetscCall(MatDuplicate(A[0],MAT_COPY_VALUES,&M));
  else PetscCall(MatCopy(A[0],M,DIFFERENT_NONZERO_PATTERN));
  PetscCall(PEPEvaluateBasis(pep,pep->eigr[idx],0,coeffs,NULL));
  PetscCall(MatScale(M,coeffs[0]));
  for (i=1;i<nmat;i++) PetscCall(MatAXPY(M,coeffs[i],A[i],(ini)?str:SUBSET_NONZERO_PATTERN));
  PetscCall(PEPEvaluateFunctionDerivatives(pep,pep->eigr[idx],coeffs2));
  for (i=0;i<nmat && PetscAbsScalar(coeffs2[i])==0.0;i++);
  PetscCall(MatMult(A[i],v,w));
  if (coeffs2[i]!=1.0) PetscCall(VecScale(w,coeffs2[i]));
  for (i++;i<nmat;i++) {
    PetscCall(MatMult(A[i],v,t));
    PetscCall(VecAXPY(w,coeffs2[i],t));
  }
  switch (scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
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
        PetscCall(PetscMalloc1(pep->n,&cols2));
        for (i=0;i<pep->n;i++) cols2[i]=i;
      }
      PetscCall(VecScatterBegin(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecGetArrayRead(ctx->nv,&array));
      if (rank==size-1) {
        PetscCall(MatSetValues(*T,1,&mg,pep->n,cols2,array,INSERT_VALUES));
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
  case PEP_REFINE_SCHEME_SCHUR:
    fctx->M2 = w;
    fctx->M3 = v;
    fctx->m3 = 0.0;
    for (i=1;i<nmat-1;i++) fctx->m3 += PetscConj(coeffs[i])*coeffs[i];
    fctx->M4 = 0.0;
    for (i=1;i<nmat-1;i++) fctx->M4 += PetscConj(coeffs[i])*coeffs2[i];
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
  case PEP_REFINE_SCHEME_MBE:
    *T = M;
    *P = M;
    break;
  }
  PetscCall(PetscFree2(coeffs,coeffs2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPNewtonRefinementSimple(PEP pep,PetscInt *maxits,PetscReal tol,PetscInt k)
{
  PetscInt             i,n,its,idx=0,*idx_sc,*its_sc,color,*fail_sc;
  PetscMPIInt          rank,size;
  Mat                  Mt=NULL,T=NULL,P=NULL;
  MPI_Comm             comm;
  Vec                  r,v,dv,rr=NULL,dvv=NULL,t[2];
  PetscScalar          *array2,deig=0.0,tt[2],ttt;
  const PetscScalar    *array;
  PetscReal            norm,error;
  PetscBool            ini=PETSC_TRUE,sc_pend,solved=PETSC_FALSE;
  PEPSimpNRefctx       *ctx;
  PEP_REFINES_MATSHELL *fctx=NULL;
  KSPConvergedReason   reason;

  PetscFunctionBegin;
  PetscCall(PetscLogEventBegin(PEP_Refine,pep,0,0,0));
  PetscCall(PEPSimpleNRefSetUp(pep,&ctx));
  its = (maxits)?*maxits:NREF_MAXIT;
  if (!pep->refineksp) PetscCall(PEPRefineGetKSP(pep,&pep->refineksp));
  if (pep->npart==1) PetscCall(BVGetColumn(pep->V,0,&v));
  else v = ctx->v;
  PetscCall(VecDuplicate(v,&ctx->w));
  PetscCall(VecDuplicate(v,&r));
  PetscCall(VecDuplicate(v,&dv));
  PetscCall(VecDuplicate(v,&t[0]));
  PetscCall(VecDuplicate(v,&t[1]));
  if (pep->npart==1) {
    PetscCall(BVRestoreColumn(pep->V,0,&v));
    PetscCall(PetscObjectGetComm((PetscObject)pep,&comm));
  } else PetscCall(PetscSubcommGetChild(pep->refinesubc,&comm));
  PetscCallMPI(MPI_Comm_size(comm,&size));
  PetscCallMPI(MPI_Comm_rank(comm,&rank));
  PetscCall(VecGetLocalSize(r,&n));
  PetscCall(PetscMalloc3(pep->npart,&idx_sc,pep->npart,&its_sc,pep->npart,&fail_sc));
  for (i=0;i<pep->npart;i++) fail_sc[i] = 0;
  for (i=0;i<pep->npart;i++) its_sc[i] = 0;
  color = (pep->npart==1)?0:pep->refinesubc->color;

  /* Loop performing iterative refinements */
  while (!solved) {
    for (i=0;i<pep->npart;i++) {
      sc_pend = PETSC_TRUE;
      if (its_sc[i]==0) {
        idx_sc[i] = idx++;
        if (idx_sc[i]>=k) {
          sc_pend = PETSC_FALSE;
        } else PetscCall(PEPSimpleNRefScatterEigenvector(pep,ctx,i,idx_sc[i]));
      }  else { /* Gather Eigenpair from subcommunicator i */
        PetscCall(PEPSimpleNRefGatherEigenpair(pep,ctx,i,idx_sc[i],&fail_sc[i]));
      }
      while (sc_pend) {
        if (!fail_sc[i]) PetscCall(PEPComputeError(pep,idx_sc[i],PEP_ERROR_BACKWARD,&error));
        if (error<=tol || its_sc[i]>=its || fail_sc[i]) {
          idx_sc[i] = idx++;
          its_sc[i] = 0;
          fail_sc[i] = 0;
          if (idx_sc[i]<k) PetscCall(PEPSimpleNRefScatterEigenvector(pep,ctx,i,idx_sc[i]));
        } else {
          sc_pend = PETSC_FALSE;
          its_sc[i]++;
        }
        if (idx_sc[i]>=k) sc_pend = PETSC_FALSE;
      }
    }
    solved = PETSC_TRUE;
    for (i=0;i<pep->npart&&solved;i++) solved = PetscNot(idx_sc[i]<k);
    if (idx_sc[color]<k) {
#if !defined(PETSC_USE_COMPLEX)
      PetscCheck(pep->eigi[idx_sc[color]]==0.0,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Simple Refinement not implemented in real scalars for complex eigenvalues");
#endif
      if (pep->npart==1) PetscCall(BVGetColumn(pep->V,idx_sc[color],&v));
      else v = ctx->v;
      PetscCall(PEPSimpleNRefSetUpSystem(pep,ctx->A,ctx,idx_sc[color],&Mt,&T,&P,ini,t[0],v));
      PetscCall(PEP_KSPSetOperators(pep->refineksp,T,P));
      if (ini) {
        PetscCall(KSPSetFromOptions(pep->refineksp));
        if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
          PetscCall(MatCreateVecs(T,&dvv,NULL));
          PetscCall(VecDuplicate(dvv,&rr));
        }
        ini = PETSC_FALSE;
      }

      switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
        PetscCall(MatMult(Mt,v,r));
        PetscCall(VecGetArrayRead(r,&array));
        if (rank==size-1) {
          PetscCall(VecGetArray(rr,&array2));
          PetscCall(PetscArraycpy(array2,array,n));
          array2[n] = 0.0;
          PetscCall(VecRestoreArray(rr,&array2));
        } else PetscCall(VecPlaceArray(rr,array));
        PetscCall(KSPSolve(pep->refineksp,rr,dvv));
        PetscCall(KSPGetConvergedReason(pep->refineksp,&reason));
        if (reason>0) {
          if (rank != size-1) PetscCall(VecResetArray(rr));
          PetscCall(VecRestoreArrayRead(r,&array));
          PetscCall(VecGetArrayRead(dvv,&array));
          PetscCall(VecPlaceArray(dv,array));
          PetscCall(VecAXPY(v,-1.0,dv));
          PetscCall(VecNorm(v,NORM_2,&norm));
          PetscCall(VecScale(v,1.0/norm));
          PetscCall(VecResetArray(dv));
          if (rank==size-1) pep->eigr[idx_sc[color]] -= array[n];
          PetscCall(VecRestoreArrayRead(dvv,&array));
        } else fail_sc[color] = 1;
        break;
      case PEP_REFINE_SCHEME_MBE:
        PetscCall(MatMult(T,v,r));
        /* Mixed block elimination */
        PetscCall(VecConjugate(v));
        PetscCall(KSPSolveTranspose(pep->refineksp,v,t[0]));
        PetscCall(KSPGetConvergedReason(pep->refineksp,&reason));
        if (reason>0) {
          PetscCall(VecConjugate(t[0]));
          PetscCall(VecDot(ctx->w,t[0],&tt[0]));
          PetscCall(KSPSolve(pep->refineksp,ctx->w,t[1]));
          PetscCall(KSPGetConvergedReason(pep->refineksp,&reason));
          if (reason>0) {
            PetscCall(VecDot(t[1],v,&tt[1]));
            PetscCall(VecDot(r,t[0],&ttt));
            tt[0] = ttt/tt[0];
            PetscCall(VecAXPY(r,-tt[0],ctx->w));
            PetscCall(KSPSolve(pep->refineksp,r,dv));
            PetscCall(KSPGetConvergedReason(pep->refineksp,&reason));
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
          pep->eigr[idx_sc[color]] -= deig;
          fail_sc[color] = 0;
        } else {
          PetscCall(VecConjugate(v));
          fail_sc[color] = 1;
        }
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        fail_sc[color] = 1;
        PetscCall(MatShellGetContext(T,&fctx));
        if (fctx->M4!=0.0) {
          PetscCall(MatMult(fctx->M1,v,r));
          PetscCall(KSPSolve(pep->refineksp,r,dv));
          PetscCall(KSPGetConvergedReason(pep->refineksp,&reason));
          if (reason>0) {
            PetscCall(VecDot(dv,v,&deig));
            deig *= -fctx->m3/fctx->M4;
            PetscCall(VecAXPY(v,-1.0,dv));
            PetscCall(VecNorm(v,NORM_2,&norm));
            PetscCall(VecScale(v,1.0/norm));
            pep->eigr[idx_sc[color]] -= deig;
            fail_sc[color] = 0;
          }
        }
        break;
      }
      if (pep->npart==1) PetscCall(BVRestoreColumn(pep->V,idx_sc[color],&v));
    }
  }
  PetscCall(VecDestroy(&t[0]));
  PetscCall(VecDestroy(&t[1]));
  PetscCall(VecDestroy(&dv));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(VecDestroy(&r));
  PetscCall(PetscFree3(idx_sc,its_sc,fail_sc));
  PetscCall(VecScatterDestroy(&ctx->nst));
  if (pep->npart>1) {
    PetscCall(VecDestroy(&ctx->vg));
    PetscCall(VecDestroy(&ctx->v));
    for (i=0;i<pep->nmat;i++) PetscCall(MatDestroy(&ctx->A[i]));
    for (i=0;i<pep->npart;i++) PetscCall(VecScatterDestroy(&ctx->scatter_id[i]));
    PetscCall(PetscFree2(ctx->A,ctx->scatter_id));
  }
  if (fctx && pep->scheme==PEP_REFINE_SCHEME_SCHUR) {
    PetscCall(MatDestroy(&P));
    PetscCall(MatDestroy(&fctx->M1));
    PetscCall(PetscFree(fctx));
  }
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    PetscCall(MatDestroy(&Mt));
    PetscCall(VecDestroy(&dvv));
    PetscCall(VecDestroy(&rr));
    PetscCall(VecDestroy(&ctx->nv));
  }
  PetscCall(MatDestroy(&T));
  PetscCall(PetscFree(ctx));
  PetscCall(PetscLogEventEnd(PEP_Refine,pep,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
