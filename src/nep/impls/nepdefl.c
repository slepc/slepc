/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepcblaslapack.h>
#include "nepdefl.h"

PetscErrorCode NEPDeflationGetInvariantPair(NEP_EXT_OP extop,BV *X,Mat *H)
{
  PetscFunctionBegin;
  if (X) *X = extop->X;
  if (H) PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,extop->szd+1,extop->szd+1,extop->H,H));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationExtendInvariantPair(NEP_EXT_OP extop,Vec u,PetscScalar lambda,PetscInt k)
{
  Vec            uu;
  PetscInt       ld,i;
  PetscMPIInt    np;
  PetscReal      norm;

  PetscFunctionBegin;
  PetscCall(BVGetColumn(extop->X,k,&uu));
  ld = extop->szd+1;
  PetscCall(NEPDeflationCopyToExtendedVec(extop,uu,extop->H+k*ld,u,PETSC_TRUE));
  PetscCall(BVRestoreColumn(extop->X,k,&uu));
  PetscCall(BVNormColumn(extop->X,k,NORM_2,&norm));
  PetscCall(BVScaleColumn(extop->X,k,1.0/norm));
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)u),&np));
  for (i=0;i<k;i++) extop->H[k*ld+i] *= PetscSqrtReal(np)/norm;
  extop->H[k*(ld+1)] = lambda;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationExtractEigenpair(NEP_EXT_OP extop,PetscInt k,Vec u,PetscScalar lambda,DS ds)
{
  Mat            A,H;
  PetscInt       ldh=extop->szd+1,ldds,k1=k+1;
  PetscScalar    *eigr,*eigi,*t,*Z;
  Vec            x;

  PetscFunctionBegin;
  PetscCall(NEPDeflationExtendInvariantPair(extop,u,lambda,k));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,k1,k1,extop->H,&H));
  PetscCall(MatDenseSetLDA(H,ldh));
  PetscCall(PetscCalloc3(k1,&eigr,k1,&eigi,extop->szd,&t));
  PetscCall(DSReset(ds));
  PetscCall(DSSetType(ds,DSNHEP));
  PetscCall(DSAllocate(ds,ldh));
  PetscCall(DSGetLeadingDimension(ds,&ldds));
  PetscCall(DSSetDimensions(ds,k1,0,k1));
  PetscCall(DSGetMat(ds,DS_MAT_A,&A));
  PetscCall(MatCopy(H,A,SAME_NONZERO_PATTERN));
  PetscCall(DSRestoreMat(ds,DS_MAT_A,&A));
  PetscCall(MatDestroy(&H));
  PetscCall(DSSolve(ds,eigr,eigi));
  PetscCall(DSVectors(ds,DS_MAT_X,&k,NULL));
  PetscCall(DSGetArray(ds,DS_MAT_X,&Z));
  PetscCall(BVMultColumn(extop->X,1.0,Z[k*ldds+k],k,Z+k*ldds));
  PetscCall(DSRestoreArray(ds,DS_MAT_X,&Z));
  PetscCall(BVGetColumn(extop->X,k,&x));
  PetscCall(NEPDeflationCopyToExtendedVec(extop,x,t,u,PETSC_FALSE));
  PetscCall(BVRestoreColumn(extop->X,k,&x));
  PetscCall(PetscFree3(eigr,eigi,t));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationCopyToExtendedVec(NEP_EXT_OP extop,Vec v,PetscScalar *a,Vec vex,PetscBool back)
{
  PetscMPIInt    np,rk,count;
  PetscScalar    *array1,*array2;
  PetscInt       nloc;

  PetscFunctionBegin;
  if (extop->szd) {
    PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)vex),&rk));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)vex),&np));
    PetscCall(BVGetSizes(extop->nep->V,&nloc,NULL,NULL));
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
        PetscCall(PetscArraycpy(a,array2+nloc,extop->szd));
        PetscCall(PetscMPIIntCast(extop->szd,&count));
        PetscCallMPI(MPI_Bcast(a,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)vex)));
      } else {
        PetscCall(PetscArraycpy(array2+nloc,a,extop->szd));
        PetscCall(PetscMPIIntCast(extop->szd,&count));
        PetscCallMPI(MPI_Bcast(array2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)vex)));
      }
      PetscCall(VecRestoreArray(vex,&array2));
    }
  } else {
    if (back) PetscCall(VecCopy(vex,v));
    else PetscCall(VecCopy(v,vex));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationCreateVec(NEP_EXT_OP extop,Vec *v)
{
  PetscInt       nloc;
  Vec            u;
  VecType        type;

  PetscFunctionBegin;
  if (extop->szd) {
    PetscCall(BVGetColumn(extop->nep->V,0,&u));
    PetscCall(VecGetType(u,&type));
    PetscCall(BVRestoreColumn(extop->nep->V,0,&u));
    PetscCall(VecCreate(PetscObjectComm((PetscObject)extop->nep),v));
    PetscCall(VecSetType(*v,type));
    PetscCall(BVGetSizes(extop->nep->V,&nloc,NULL,NULL));
    nloc += extop->szd;
    PetscCall(VecSetSizes(*v,nloc,PETSC_DECIDE));
  } else PetscCall(BVCreateVec(extop->nep->V,v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationCreateBV(NEP_EXT_OP extop,PetscInt sz,BV *V)
{
  PetscInt           nloc;
  BVType             type;
  BVOrthogType       otype;
  BVOrthogRefineType oref;
  PetscReal          oeta;
  BVOrthogBlockType  oblock;
  NEP                nep=extop->nep;
  VecType            vtype;

  PetscFunctionBegin;
  if (extop->szd) {
    PetscCall(BVGetSizes(nep->V,&nloc,NULL,NULL));
    PetscCall(BVCreate(PetscObjectComm((PetscObject)nep),V));
    PetscCall(BVSetSizes(*V,nloc+extop->szd,PETSC_DECIDE,sz));
    PetscCall(BVGetType(nep->V,&type));
    PetscCall(BVSetType(*V,type));
    PetscCall(BVGetVecType(nep->V,&vtype));
    PetscCall(BVSetVecType(*V,vtype));
    PetscCall(BVGetOrthogonalization(nep->V,&otype,&oref,&oeta,&oblock));
    PetscCall(BVSetOrthogonalization(*V,otype,oref,oeta,oblock));
    PetscCall(PetscObjectStateIncrease((PetscObject)*V));
  } else PetscCall(BVDuplicateResize(nep->V,sz,V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationSetRandomVec(NEP_EXT_OP extop,Vec v)
{
  PetscInt       n,next,i;
  PetscRandom    rand;
  PetscScalar    *array;
  PetscMPIInt    nn,np;

  PetscFunctionBegin;
  PetscCall(BVGetRandomContext(extop->nep->V,&rand));
  PetscCall(VecSetRandom(v,rand));
  if (extop->szd) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)v),&np));
    PetscCall(BVGetSizes(extop->nep->V,&n,NULL,NULL));
    PetscCall(VecGetLocalSize(v,&next));
    PetscCall(VecGetArray(v,&array));
    for (i=n+extop->n;i<next;i++) array[i] = 0.0;
    for (i=n;i<n+extop->n;i++) array[i] /= PetscSqrtReal(np);
    PetscCall(PetscMPIIntCast(extop->n,&nn));
    PetscCallMPI(MPI_Bcast(array+n,nn,MPIU_SCALAR,0,PetscObjectComm((PetscObject)v)));
    PetscCall(VecRestoreArray(v,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationEvaluateBasisMat(NEP_EXT_OP extop,PetscInt idx,PetscBool hat,PetscScalar *bval,PetscScalar *Hj,PetscScalar *Hjprev)
{
  PetscInt       i,k,n=extop->n,ldhj=extop->szd,ldh=extop->szd+1;
  PetscScalar    sone=1.0,zero=0.0;
  PetscBLASInt   ldh_,ldhj_,n_;

  PetscFunctionBegin;
  i = (idx<0)?extop->szd*extop->szd*(-idx):extop->szd*extop->szd;
  PetscCall(PetscArrayzero(Hj,i));
  PetscCall(PetscBLASIntCast(ldhj+1,&ldh_));
  PetscCall(PetscBLASIntCast(ldhj,&ldhj_));
  PetscCall(PetscBLASIntCast(n,&n_));
  if (idx<1) {
    if (!hat) for (i=0;i<extop->n;i++) Hj[i+i*ldhj] = 1.0;
    else for (i=0;i<extop->n;i++) Hj[i+i*ldhj] = 0.0;
  } else {
      for (i=0;i<n;i++) extop->H[i*ldh+i] -= extop->bc[idx-1];
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->H,&ldh_,Hjprev,&ldhj_,&zero,Hj,&ldhj_));
      for (i=0;i<n;i++) extop->H[i*ldh+i] += extop->bc[idx-1];
      if (hat) for (i=0;i<n;i++) Hj[i*(ldhj+1)] += bval[idx-1];
  }
  if (idx<0) {
    idx = -idx;
    for (k=1;k<idx;k++) {
      for (i=0;i<n;i++) extop->H[i*ldh+i] -= extop->bc[k-1];
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->H,&ldh_,Hj+(k-1)*ldhj*ldhj,&ldhj_,&zero,Hj+k*ldhj*ldhj,&ldhj_));
      for (i=0;i<n;i++) extop->H[i*ldh+i] += extop->bc[k-1];
      if (hat) for (i=0;i<n;i++) Hj[i*(ldhj+1)] += bval[k-1];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationLocking(NEP_EXT_OP extop,Vec u,PetscScalar lambda)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(NEPDeflationExtendInvariantPair(extop,u,lambda,extop->n));
  extop->n++;
  PetscCall(BVSetActiveColumns(extop->X,0,extop->n));
  if (extop->n <= extop->szd) {
    /* update XpX */
    PetscCall(BVDotColumn(extop->X,extop->n-1,extop->XpX+(extop->n-1)*extop->szd));
    extop->XpX[(extop->n-1)*(1+extop->szd)] = 1.0;
    for (i=0;i<extop->n-1;i++) extop->XpX[i*extop->szd+extop->n-1] = PetscConj(extop->XpX[(extop->n-1)*extop->szd+i]);
    /* determine minimality index */
    extop->midx = PetscMin(extop->max_midx,extop->n);
    /* polynominal basis coefficients */
    for (i=0;i<extop->midx;i++) extop->bc[i] = extop->nep->target;
    /* evaluate the polynomial basis in H */
    PetscCall(NEPDeflationEvaluateBasisMat(extop,-extop->midx,PETSC_FALSE,NULL,extop->Hj,NULL));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationEvaluateHatFunction(NEP_EXT_OP extop, PetscInt idx,PetscScalar lambda,PetscScalar *y,PetscScalar *hfj,PetscScalar *hfjp,PetscInt ld)
{
  PetscInt          i,j,k,off,ini,fin,sz,ldh,n=extop->n;
  Mat               A,B;
  PetscScalar       *array;
  const PetscScalar *barray;

  PetscFunctionBegin;
  if (idx<0) {ini = 0; fin = extop->nep->nt;}
  else {ini = idx; fin = idx+1;}
  if (y) sz = hfjp?n+2:n+1;
  else sz = hfjp?3*n:2*n;
  ldh = extop->szd+1;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,sz,sz,NULL,&A));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,sz,sz,NULL,&B));
  PetscCall(MatDenseGetArray(A,&array));
  for (j=0;j<n;j++)
    for (i=0;i<n;i++) array[j*sz+i] = extop->H[j*ldh+i];
  PetscCall(MatDenseRestoreArrayWrite(A,&array));
  if (y) {
    PetscCall(MatDenseGetArray(A,&array));
    array[extop->n*(sz+1)] = lambda;
    if (hfjp) { array[(n+1)*sz+n] = 1.0; array[(n+1)*sz+n+1] = lambda;}
    for (i=0;i<n;i++) array[n*sz+i] = y[i];
    PetscCall(MatDenseRestoreArrayWrite(A,&array));
    for (j=ini;j<fin;j++) {
      PetscCall(FNEvaluateFunctionMat(extop->nep->f[j],A,B));
      PetscCall(MatDenseGetArrayRead(B,&barray));
      for (i=0;i<n;i++) hfj[j*ld+i] = barray[n*sz+i];
      if (hfjp) for (i=0;i<n;i++) hfjp[j*ld+i] = barray[(n+1)*sz+i];
      PetscCall(MatDenseRestoreArrayRead(B,&barray));
    }
  } else {
    off = idx<0?ld*n:0;
    PetscCall(MatDenseGetArray(A,&array));
    for (k=0;k<n;k++) {
      array[(n+k)*sz+k] = 1.0;
      array[(n+k)*sz+n+k] = lambda;
    }
    if (hfjp) for (k=0;k<n;k++) {
      array[(2*n+k)*sz+n+k] = 1.0;
      array[(2*n+k)*sz+2*n+k] = lambda;
    }
    PetscCall(MatDenseRestoreArray(A,&array));
    for (j=ini;j<fin;j++) {
      PetscCall(FNEvaluateFunctionMat(extop->nep->f[j],A,B));
      PetscCall(MatDenseGetArrayRead(B,&barray));
      for (i=0;i<n;i++) for (k=0;k<n;k++) hfj[j*off+i*ld+k] = barray[n*sz+i*sz+k];
      if (hfjp) for (k=0;k<n;k++) for (i=0;i<n;i++) hfjp[j*off+i*ld+k] = barray[2*n*sz+i*sz+k];
      PetscCall(MatDenseRestoreArrayRead(B,&barray));
    }
  }
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatMult_NEPDeflation(Mat M,Vec x,Vec y)
{
  NEP_DEF_MATSHELL  *matctx;
  NEP_EXT_OP        extop;
  Vec               x1,y1;
  PetscScalar       *yy,sone=1.0,zero=0.0;
  const PetscScalar *xx;
  PetscInt          nloc,i;
  PetscMPIInt       np;
  PetscBLASInt      n_,one=1,szd_;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)M),&np));
  PetscCall(MatShellGetContext(M,&matctx));
  extop = matctx->extop;
  if (extop->ref) PetscCall(VecZeroEntries(y));
  if (extop->szd) {
    x1 = matctx->w[0]; y1 = matctx->w[1];
    PetscCall(VecGetArrayRead(x,&xx));
    PetscCall(VecPlaceArray(x1,xx));
    PetscCall(VecGetArray(y,&yy));
    PetscCall(VecPlaceArray(y1,yy));
    PetscCall(MatMult(matctx->T,x1,y1));
    if (!extop->ref && extop->n) {
      PetscCall(VecGetLocalSize(x1,&nloc));
      /* copy for avoiding warning of constant array xx */
      for (i=0;i<extop->n;i++) matctx->work[i] = xx[nloc+i]*PetscSqrtReal(np);
      PetscCall(BVMultVec(matctx->U,1.0,1.0,y1,matctx->work));
      PetscCall(BVDotVec(extop->X,x1,matctx->work));
      PetscCall(PetscBLASIntCast(extop->n,&n_));
      PetscCall(PetscBLASIntCast(extop->szd,&szd_));
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,matctx->A,&szd_,matctx->work,&one,&zero,yy+nloc,&one));
      PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,matctx->B,&szd_,xx+nloc,&one,&sone,yy+nloc,&one));
      for (i=0;i<extop->n;i++) yy[nloc+i] /= PetscSqrtReal(np);
    }
    PetscCall(VecResetArray(x1));
    PetscCall(VecRestoreArrayRead(x,&xx));
    PetscCall(VecResetArray(y1));
    PetscCall(VecRestoreArray(y,&yy));
  } else PetscCall(MatMult(matctx->T,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatCreateVecs_NEPDeflation(Mat M,Vec *right,Vec *left)
{
  NEP_DEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (right) PetscCall(VecDuplicate(matctx->w[0],right));
  if (left) PetscCall(VecDuplicate(matctx->w[0],left));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatDestroy_NEPDeflation(Mat M)
{
  NEP_DEF_MATSHELL *matctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(M,&matctx));
  if (matctx->extop->szd) {
    PetscCall(BVDestroy(&matctx->U));
    PetscCall(PetscFree4(matctx->hfj,matctx->work,matctx->A,matctx->B));
    PetscCall(VecDestroy(&matctx->w[0]));
    PetscCall(VecDestroy(&matctx->w[1]));
  }
  if (matctx->P != matctx->T) PetscCall(MatDestroy(&matctx->P));
  PetscCall(MatDestroy(&matctx->T));
  PetscCall(PetscFree(matctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationEvaluateBasis(NEP_EXT_OP extop,PetscScalar lambda,PetscInt n,PetscScalar *val,PetscBool jacobian)
{
  PetscScalar p;
  PetscInt    i;

  PetscFunctionBegin;
  if (!jacobian) {
    val[0] = 1.0;
    for (i=1;i<extop->n;i++) val[i] = val[i-1]*(lambda-extop->bc[i-1]);
  } else {
    val[0] = 0.0;
    p = 1.0;
    for (i=1;i<extop->n;i++) {
      val[i] = val[i-1]*(lambda-extop->bc[i-1])+p;
      p *= (lambda-extop->bc[i-1]);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationComputeShellMat(NEP_EXT_OP extop,PetscScalar lambda,PetscBool jacobian,Mat *M)
{
  NEP_DEF_MATSHELL *matctx,*matctxC;
  PetscInt         nloc,mloc,n=extop->n,j,i,szd=extop->szd,ldh=szd+1,k;
  Mat              F,Mshell,Mcomp;
  PetscBool        ini=PETSC_FALSE;
  PetscScalar      *hf,*hfj,*hfjp,sone=1.0,*hH,*hHprev,*pts,*B,*A,*Hj=extop->Hj,*basisv,zero=0.0;
  PetscBLASInt     n_,info,szd_;
  size_t           len;

  PetscFunctionBegin;
  if (!M) Mshell = jacobian?extop->MJ:extop->MF;
  else Mshell = *M;
  Mcomp = jacobian?extop->MF:extop->MJ;
  if (!Mshell) {
    ini = PETSC_TRUE;
    PetscCall(PetscNew(&matctx));
    PetscCall(MatGetLocalSize(extop->nep->function,&mloc,&nloc));
    nloc += szd; mloc += szd;
    PetscCall(MatCreateShell(PetscObjectComm((PetscObject)extop->nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&Mshell));
    PetscCall(MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_NEPDeflation));
    PetscCall(MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_NEPDeflation));
    PetscCall(MatShellSetOperation(Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_NEPDeflation));
    matctx->nep = extop->nep;
    matctx->extop = extop;
    if (!M) {
      if (jacobian) { matctx->jacob = PETSC_TRUE; matctx->T = extop->nep->jacobian; extop->MJ = Mshell; }
      else { matctx->jacob = PETSC_FALSE; matctx->T = extop->nep->function; extop->MF = Mshell; }
      PetscCall(PetscObjectReference((PetscObject)matctx->T));
      if (!jacobian) {
        if (extop->nep->function_pre && extop->nep->function_pre != extop->nep->function) {
          matctx->P = extop->nep->function_pre;
          PetscCall(PetscObjectReference((PetscObject)matctx->P));
        } else matctx->P = matctx->T;
      }
    } else {
      matctx->jacob = jacobian;
      PetscCall(MatDuplicate(jacobian?extop->nep->jacobian:extop->nep->function,MAT_DO_NOT_COPY_VALUES,&matctx->T));
      *M = Mshell;
      if (!jacobian) {
        if (extop->nep->function_pre && extop->nep->function_pre != extop->nep->function) PetscCall(MatDuplicate(extop->nep->function_pre,MAT_DO_NOT_COPY_VALUES,&matctx->P));
        else matctx->P = matctx->T;
      }
    }
    if (szd) {
      PetscCall(BVCreateVec(extop->nep->V,matctx->w));
      PetscCall(VecDuplicate(matctx->w[0],matctx->w+1));
      PetscCall(BVDuplicateResize(extop->nep->V,szd,&matctx->U));
      len = extop->simpU? (size_t)2*szd*szd: (size_t)2*szd*szd*extop->nep->nt;
      PetscCall(PetscMalloc4(len,&matctx->hfj,szd,&matctx->work,szd*szd,&matctx->A,szd*szd,&matctx->B));
    }
  } else PetscCall(MatShellGetContext(Mshell,&matctx));
  if (ini || matctx->theta != lambda || matctx->n != extop->n) {
    if (ini || matctx->theta != lambda) {
      if (jacobian) PetscCall(NEPComputeJacobian(extop->nep,lambda,matctx->T));
      else PetscCall(NEPComputeFunction(extop->nep,lambda,matctx->T,matctx->P));
    }
    if (n) {
      matctx->hfjset = PETSC_FALSE;
      if (!extop->simpU) {
        /* likely hfjp has been already computed */
        if (Mcomp) {
          PetscCall(MatShellGetContext(Mcomp,&matctxC));
          if (matctxC->hfjset && matctxC->theta == lambda && matctxC->n == extop->n) {
            PetscCall(PetscArraycpy(matctx->hfj,matctxC->hfj,2*extop->szd*extop->szd*extop->nep->nt));
            matctx->hfjset = PETSC_TRUE;
          }
        }
        hfj = matctx->hfj; hfjp = matctx->hfj+extop->szd*extop->szd*extop->nep->nt;
        if (!matctx->hfjset) {
          PetscCall(NEPDeflationEvaluateHatFunction(extop,-1,lambda,NULL,hfj,hfjp,n));
          matctx->hfjset = PETSC_TRUE;
        }
        PetscCall(BVSetActiveColumns(matctx->U,0,n));
        hf = jacobian?hfjp:hfj;
        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,hf,&F));
        PetscCall(BVMatMult(extop->X,extop->nep->A[0],matctx->U));
        PetscCall(BVMultInPlace(matctx->U,F,0,n));
        PetscCall(BVSetActiveColumns(extop->W,0,extop->n));
        for (j=1;j<extop->nep->nt;j++) {
          PetscCall(BVMatMult(extop->X,extop->nep->A[j],extop->W));
          PetscCall(MatDensePlaceArray(F,hf+j*n*n));
          PetscCall(BVMult(matctx->U,1.0,1.0,extop->W,F));
          PetscCall(MatDenseResetArray(F));
        }
        PetscCall(MatDestroy(&F));
      } else {
        hfj = matctx->hfj;
        PetscCall(BVSetActiveColumns(matctx->U,0,n));
        PetscCall(BVMatMult(extop->X,matctx->T,matctx->U));
        for (j=0;j<n;j++) {
          for (i=0;i<n;i++) hfj[j*n+i] = -extop->H[j*ldh+i];
          hfj[j*(n+1)] += lambda;
        }
        PetscCall(PetscBLASIntCast(n,&n_));
        PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
        PetscCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,hfj,&n_,&info));
        PetscCall(PetscFPTrapPop());
        SlepcCheckLapackInfo("trtri",info);
        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,hfj,&F));
        PetscCall(BVMultInPlace(matctx->U,F,0,n));
        if (jacobian) {
          PetscCall(NEPDeflationComputeFunction(extop,lambda,NULL));
          PetscCall(MatShellGetContext(extop->MF,&matctxC));
          PetscCall(BVMult(matctx->U,-1.0,1.0,matctxC->U,F));
        }
        PetscCall(MatDestroy(&F));
      }
      PetscCall(PetscCalloc3(n,&basisv,szd*szd,&hH,szd*szd,&hHprev));
      PetscCall(NEPDeflationEvaluateBasis(extop,lambda,n,basisv,jacobian));
      A = matctx->A;
      PetscCall(PetscArrayzero(A,szd*szd));
      if (!jacobian) for (i=0;i<n;i++) A[i*(szd+1)] = 1.0;
      for (j=0;j<n;j++)
        for (i=0;i<n;i++)
          for (k=1;k<extop->midx;k++) A[j*szd+i] += basisv[k]*PetscConj(Hj[k*szd*szd+i*szd+j]);
      PetscCall(PetscBLASIntCast(n,&n_));
      PetscCall(PetscBLASIntCast(szd,&szd_));
      B = matctx->B;
      PetscCall(PetscArrayzero(B,szd*szd));
      for (i=1;i<extop->midx;i++) {
        PetscCall(NEPDeflationEvaluateBasisMat(extop,i,PETSC_TRUE,basisv,hH,hHprev));
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->XpX,&szd_,hH,&szd_,&zero,hHprev,&szd_));
        PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,extop->Hj+szd*szd*i,&szd_,hHprev,&szd_,&sone,B,&szd_));
        pts = hHprev; hHprev = hH; hH = pts;
      }
      PetscCall(PetscFree3(basisv,hH,hHprev));
    }
    matctx->theta = lambda;
    matctx->n = extop->n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationComputeFunction(NEP_EXT_OP extop,PetscScalar lambda,Mat *F)
{
  PetscFunctionBegin;
  PetscCall(NEPDeflationComputeShellMat(extop,lambda,PETSC_FALSE,NULL));
  if (F) *F = extop->MF;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationComputeJacobian(NEP_EXT_OP extop,PetscScalar lambda,Mat *J)
{
  PetscFunctionBegin;
  PetscCall(NEPDeflationComputeShellMat(extop,lambda,PETSC_TRUE,NULL));
  if (J) *J = extop->MJ;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationSolveSetUp(NEP_EXT_OP extop,PetscScalar lambda)
{
  NEP_DEF_MATSHELL  *matctx;
  NEP_DEF_FUN_SOLVE solve;
  PetscInt          i,j,n=extop->n;
  Vec               u,tu;
  Mat               F;
  PetscScalar       snone=-1.0,sone=1.0;
  PetscBLASInt      n_,szd_,ldh_,*p,info;
  Mat               Mshell;

  PetscFunctionBegin;
  solve = extop->solve;
  if (lambda!=solve->theta || n!=solve->n) {
    PetscCall(NEPDeflationComputeShellMat(extop,lambda,PETSC_FALSE,solve->sincf?NULL:&solve->T));
    Mshell = (solve->sincf)?extop->MF:solve->T;
    PetscCall(MatShellGetContext(Mshell,&matctx));
    PetscCall(NEP_KSPSetOperators(solve->ksp,matctx->T,matctx->P));
    if (!extop->ref && n) {
      PetscCall(PetscBLASIntCast(n,&n_));
      PetscCall(PetscBLASIntCast(extop->szd,&szd_));
      PetscCall(PetscBLASIntCast(extop->szd+1,&ldh_));
      if (!extop->simpU) {
        PetscCall(BVSetActiveColumns(solve->T_1U,0,n));
        for (i=0;i<n;i++) {
          PetscCall(BVGetColumn(matctx->U,i,&u));
          PetscCall(BVGetColumn(solve->T_1U,i,&tu));
          PetscCall(KSPSolve(solve->ksp,u,tu));
          PetscCall(BVRestoreColumn(solve->T_1U,i,&tu));
          PetscCall(BVRestoreColumn(matctx->U,i,&u));
        }
        PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,n,n,solve->work,&F));
        PetscCall(BVDot(solve->T_1U,extop->X,F));
        PetscCall(MatDestroy(&F));
      } else {
        for (j=0;j<n;j++)
          for (i=0;i<n;i++) solve->work[j*n+i] = extop->XpX[j*extop->szd+i];
        for (i=0;i<n;i++) extop->H[i*ldh_+i] -= lambda;
        PetscCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&n_,&n_,&snone,extop->H,&ldh_,solve->work,&n_));
        for (i=0;i<n;i++) extop->H[i*ldh_+i] += lambda;
      }
      PetscCall(PetscArraycpy(solve->M,matctx->B,extop->szd*extop->szd));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&snone,matctx->A,&szd_,solve->work,&n_,&sone,solve->M,&szd_));
      PetscCall(PetscMalloc1(n,&p));
      PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
      PetscCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,solve->M,&szd_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      PetscCallBLAS("LAPACKgetri",LAPACKgetri_(&n_,solve->M,&szd_,p,solve->work,&n_,&info));
      SlepcCheckLapackInfo("getri",info);
      PetscCall(PetscFPTrapPop());
      PetscCall(PetscFree(p));
    }
    solve->theta = lambda;
    solve->n = n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationFunctionSolve(NEP_EXT_OP extop,Vec b,Vec x)
{
  Vec               b1,x1;
  PetscScalar       *xx,*bb,*x2,*b2,*w,*w2,snone=-1.0,sone=1.0,zero=0.0;
  NEP_DEF_MATSHELL  *matctx;
  NEP_DEF_FUN_SOLVE solve=extop->solve;
  PetscBLASInt      one=1,szd_,n_,ldh_;
  PetscInt          nloc,i;
  PetscMPIInt       np,count;

  PetscFunctionBegin;
  if (extop->ref) PetscCall(VecZeroEntries(x));
  if (extop->szd) {
    x1 = solve->w[0]; b1 = solve->w[1];
    PetscCall(VecGetArray(x,&xx));
    PetscCall(VecPlaceArray(x1,xx));
    PetscCall(VecGetArray(b,&bb));
    PetscCall(VecPlaceArray(b1,bb));
  } else {
    b1 = b; x1 = x;
  }
  PetscCall(KSPSolve(extop->solve->ksp,b1,x1));
  if (!extop->ref && extop->n && extop->szd) {
    PetscCall(PetscBLASIntCast(extop->szd,&szd_));
    PetscCall(PetscBLASIntCast(extop->n,&n_));
    PetscCall(PetscBLASIntCast(extop->szd+1,&ldh_));
    PetscCall(BVGetSizes(extop->nep->V,&nloc,NULL,NULL));
    PetscCall(PetscMalloc2(extop->n,&b2,extop->n,&x2));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)b),&np));
    for (i=0;i<extop->n;i++) b2[i] = bb[nloc+i]*PetscSqrtReal(np);
    w = solve->work; w2 = solve->work+extop->n;
    PetscCall(MatShellGetContext(solve->sincf?extop->MF:solve->T,&matctx));
    PetscCall(PetscArraycpy(w2,b2,extop->n));
    PetscCall(BVDotVec(extop->X,x1,w));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&snone,matctx->A,&szd_,w,&one,&sone,w2,&one));
    PetscCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,solve->M,&szd_,w2,&one,&zero,x2,&one));
    if (extop->simpU) {
      for (i=0;i<extop->n;i++) extop->H[i+i*(extop->szd+1)] -= solve->theta;
      for (i=0;i<extop->n;i++) w[i] = x2[i];
      PetscCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&one,&snone,extop->H,&ldh_,w,&n_));
      for (i=0;i<extop->n;i++) extop->H[i+i*(extop->szd+1)] += solve->theta;
      PetscCall(BVMultVec(extop->X,-1.0,1.0,x1,w));
    } else PetscCall(BVMultVec(solve->T_1U,-1.0,1.0,x1,x2));
    for (i=0;i<extop->n;i++) xx[i+nloc] = x2[i]/PetscSqrtReal(np);
    PetscCall(PetscMPIIntCast(extop->n,&count));
    PetscCallMPI(MPI_Bcast(xx+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)b)));
  }
  if (extop->szd) {
    PetscCall(VecResetArray(x1));
    PetscCall(VecRestoreArray(x,&xx));
    PetscCall(VecResetArray(b1));
    PetscCall(VecRestoreArray(b,&bb));
    if (!extop->ref && extop->n) PetscCall(PetscFree2(b2,x2));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationSetRefine(NEP_EXT_OP extop,PetscBool ref)
{
  PetscFunctionBegin;
  extop->ref = ref;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationReset(NEP_EXT_OP extop)
{
  PetscInt          j;
  NEP_DEF_FUN_SOLVE solve;

  PetscFunctionBegin;
  if (!extop) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscFree(extop->H));
  PetscCall(BVDestroy(&extop->X));
  if (extop->szd) {
    PetscCall(VecDestroy(&extop->w));
    PetscCall(PetscFree3(extop->Hj,extop->XpX,extop->bc));
    PetscCall(BVDestroy(&extop->W));
  }
  PetscCall(MatDestroy(&extop->MF));
  PetscCall(MatDestroy(&extop->MJ));
  if (extop->solve) {
    solve = extop->solve;
    if (extop->szd) {
      if (!extop->simpU) PetscCall(BVDestroy(&solve->T_1U));
      PetscCall(PetscFree2(solve->M,solve->work));
      PetscCall(VecDestroy(&solve->w[0]));
      PetscCall(VecDestroy(&solve->w[1]));
    }
    if (!solve->sincf) PetscCall(MatDestroy(&solve->T));
    PetscCall(PetscFree(extop->solve));
  }
  if (extop->proj) {
    if (extop->szd) {
      for (j=0;j<extop->nep->nt;j++) PetscCall(MatDestroy(&extop->proj->V1pApX[j]));
      PetscCall(MatDestroy(&extop->proj->XpV1));
      PetscCall(PetscFree3(extop->proj->V2,extop->proj->V1pApX,extop->proj->work));
      PetscCall(VecDestroy(&extop->proj->w));
      PetscCall(BVDestroy(&extop->proj->V1));
    }
    PetscCall(PetscFree(extop->proj));
  }
  PetscCall(PetscFree(extop));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationInitialize(NEP nep,BV X,KSP ksp,PetscBool sincfun,PetscInt sz,NEP_EXT_OP *extop)
{
  NEP_EXT_OP        op;
  NEP_DEF_FUN_SOLVE solve;
  PetscInt          szd;
  Vec               x;

  PetscFunctionBegin;
  PetscCall(NEPDeflationReset(*extop));
  PetscCall(PetscNew(&op));
  *extop  = op;
  op->nep = nep;
  op->n   = 0;
  op->szd = szd = sz-1;
  op->max_midx = PetscMin(MAX_MINIDX,szd);
  op->X = X;
  if (!X) PetscCall(BVDuplicateResize(nep->V,sz,&op->X));
  else PetscCall(PetscObjectReference((PetscObject)X));
  PetscCall(PetscCalloc1(sz*sz,&(op)->H));
  if (op->szd) {
    PetscCall(BVGetColumn(op->X,0,&x));
    PetscCall(VecDuplicate(x,&op->w));
    PetscCall(BVRestoreColumn(op->X,0,&x));
    op->simpU = PETSC_FALSE;
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
      /* undocumented option to use the simple expression for U = T*X*inv(lambda-H) */
      PetscCall(PetscOptionsGetBool(NULL,NULL,"-nep_deflation_simpleu",&op->simpU,NULL));
    } else {
      op->simpU = PETSC_TRUE;
    }
    PetscCall(PetscCalloc3(szd*szd*op->max_midx,&(op)->Hj,szd*szd,&(op)->XpX,szd,&op->bc));
    PetscCall(BVDuplicateResize(op->X,op->szd,&op->W));
  }
  if (ksp) {
    PetscCall(PetscNew(&solve));
    op->solve    = solve;
    solve->ksp   = ksp;
    solve->sincf = sincfun;
    solve->n     = -1;
    if (op->szd) {
      if (!op->simpU) PetscCall(BVDuplicateResize(nep->V,szd,&solve->T_1U));
      PetscCall(PetscMalloc2(szd*szd,&solve->M,2*szd*szd,&solve->work));
      PetscCall(BVCreateVec(nep->V,&solve->w[0]));
      PetscCall(VecDuplicate(solve->w[0],&solve->w[1]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode NEPDeflationDSNEPComputeMatrix(DS ds,PetscScalar lambda,PetscBool deriv,DSMatType mat,void *ctx)
{
  Mat               A,Ei;
  PetscScalar       *T,*w1,*w2,*w=NULL,*ww,*hH,*hHprev,*pts;
  PetscScalar       alpha,alpha2,*AB,sone=1.0,zero=0.0,*basisv,s;
  const PetscScalar *E;
  PetscInt          i,ldds,nwork=0,szd,nv,j,k,n;
  PetscBLASInt      inc=1,nv_,ldds_,dim_,szdk,szd_,n_,ldh_;
  PetscMPIInt       np;
  NEP_DEF_PROJECT   proj=(NEP_DEF_PROJECT)ctx;
  NEP_EXT_OP        extop=proj->extop;
  NEP               nep=extop->nep;

  PetscFunctionBegin;
  PetscCall(DSGetDimensions(ds,&nv,NULL,NULL,NULL));
  PetscCall(DSGetLeadingDimension(ds,&ldds));
  PetscCall(DSGetMat(ds,mat,&A));
  PetscCall(MatZeroEntries(A));
  /* mat = V1^*T(lambda)V1 */
  for (i=0;i<nep->nt;i++) {
    if (deriv) PetscCall(FNEvaluateDerivative(nep->f[i],lambda,&alpha));
    else PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
    PetscCall(DSGetMat(ds,DSMatExtra[i],&Ei));
    PetscCall(MatAXPY(A,alpha,Ei,SAME_NONZERO_PATTERN));
    PetscCall(DSRestoreMat(ds,DSMatExtra[i],&Ei));
  }
  PetscCall(DSRestoreMat(ds,mat,&A));
  if (!extop->ref && extop->n) {
    PetscCall(DSGetArray(ds,mat,&T));
    n = extop->n;
    szd = extop->szd;
    PetscCall(PetscArrayzero(proj->work,proj->lwork));
    PetscCall(PetscBLASIntCast(nv,&nv_));
    PetscCall(PetscBLASIntCast(n,&n_));
    PetscCall(PetscBLASIntCast(ldds,&ldds_));
    PetscCall(PetscBLASIntCast(szd,&szd_));
    PetscCall(PetscBLASIntCast(proj->dim,&dim_));
    PetscCall(PetscBLASIntCast(extop->szd+1,&ldh_));
    w1 = proj->work; w2 = proj->work+proj->dim*proj->dim;
    nwork += 2*proj->dim*proj->dim;

    /* mat = mat + V1^*U(lambda)V2 */
    for (i=0;i<nep->nt;i++) {
      if (extop->simpU) {
        if (deriv) PetscCall(FNEvaluateDerivative(nep->f[i],lambda,&alpha));
        else PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha));
        ww = w1; w = w2;
        PetscCall(PetscArraycpy(ww,proj->V2,szd*nv));
        PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ds),&np));
        for (j=0;j<szd*nv;j++) ww[j] *= PetscSqrtReal(np);
        for (j=0;j<n;j++) extop->H[j*ldh_+j] -= lambda;
        alpha = -alpha;
        PetscCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha,extop->H,&ldh_,ww,&szd_));
        if (deriv) {
          PetscCall(PetscBLASIntCast(szd*nv,&szdk));
          PetscCall(FNEvaluateFunction(nep->f[i],lambda,&alpha2));
          PetscCall(PetscArraycpy(w,proj->V2,szd*nv));
          for (j=0;j<szd*nv;j++) w[j] *= PetscSqrtReal(np);
          alpha2 = -alpha2;
          PetscCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha2,extop->H,&ldh_,w,&szd_));
          alpha2 = 1.0;
          PetscCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha2,extop->H,&ldh_,w,&szd_));
          PetscCallBLAS("BLASaxpy",BLASaxpy_(&szdk,&sone,w,&inc,ww,&inc));
        }
        for (j=0;j<n;j++) extop->H[j*ldh_+j] += lambda;
      } else {
        PetscCall(NEPDeflationEvaluateHatFunction(extop,i,lambda,NULL,w1,w2,szd));
        w = deriv?w2:w1; ww = deriv?w1:w2;
        PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)ds),&np));
        s = PetscSqrtReal(np);
        PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&s,w,&szd_,proj->V2,&szd_,&zero,ww,&szd_));
      }
      PetscCall(MatDenseGetArrayRead(proj->V1pApX[i],&E));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&nv_,&nv_,&n_,&sone,E,&dim_,ww,&szd_,&sone,T,&ldds_));
      PetscCall(MatDenseRestoreArrayRead(proj->V1pApX[i],&E));
    }

    /* mat = mat + V2^*A(lambda)V1 */
    basisv = proj->work+nwork; nwork += szd;
    hH     = proj->work+nwork; nwork += szd*szd;
    hHprev = proj->work+nwork; nwork += szd*szd;
    AB     = proj->work+nwork;
    PetscCall(NEPDeflationEvaluateBasis(extop,lambda,n,basisv,deriv));
    if (!deriv) for (i=0;i<n;i++) AB[i*(szd+1)] = 1.0;
    for (j=0;j<n;j++)
      for (i=0;i<n;i++)
        for (k=1;k<extop->midx;k++) AB[j*szd+i] += basisv[k]*PetscConj(extop->Hj[k*szd*szd+i*szd+j]);
    PetscCall(MatDenseGetArrayRead(proj->XpV1,&E));
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&sone,AB,&szd_,E,&szd_,&zero,w,&szd_));
    PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&nv_,&nv_,&n_,&sone,proj->V2,&szd_,w,&szd_,&sone,T,&ldds_));
    PetscCall(MatDenseRestoreArrayRead(proj->XpV1,&E));

    /* mat = mat + V2^*B(lambda)V2 */
    PetscCall(PetscArrayzero(AB,szd*szd));
    for (i=1;i<extop->midx;i++) {
      PetscCall(NEPDeflationEvaluateBasisMat(extop,i,PETSC_TRUE,basisv,hH,hHprev));
      PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->XpX,&szd_,hH,&szd_,&zero,hHprev,&szd_));
      PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,extop->Hj+szd*szd*i,&szd_,hHprev,&szd_,&sone,AB,&szd_));
      pts = hHprev; hHprev = hH; hH = pts;
    }
    PetscCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&sone,AB,&szd_,proj->V2,&szd_,&zero,w,&szd_));
    PetscCallBLAS("BLASgemm",BLASgemm_("C","N",&nv_,&nv_,&n_,&sone,proj->V2,&szd_,w,&szd_,&sone,T,&ldds_));
    PetscCall(DSRestoreArray(ds,mat,&T));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode NEPDeflationProjectOperator(NEP_EXT_OP extop,BV Vext,DS ds,PetscInt j0,PetscInt j1)
{
  PetscInt        k,j,n=extop->n,dim;
  Vec             v,ve;
  BV              V1;
  Mat             G;
  NEP             nep=extop->nep;
  NEP_DEF_PROJECT proj;

  PetscFunctionBegin;
  NEPCheckSplit(extop->nep,1);
  proj = extop->proj;
  if (!proj) {
    /* Initialize the projection data structure */
    PetscCall(PetscNew(&proj));
    extop->proj = proj;
    proj->extop = extop;
    PetscCall(BVGetSizes(Vext,NULL,NULL,&dim));
    proj->dim = dim;
    if (extop->szd) {
      proj->lwork = 3*dim*dim+2*extop->szd*extop->szd+extop->szd;
      PetscCall(PetscMalloc3(dim*extop->szd,&proj->V2,nep->nt,&proj->V1pApX,proj->lwork,&proj->work));
      for (j=0;j<nep->nt;j++) PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,proj->dim,extop->szd,NULL,&proj->V1pApX[j]));
      PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,extop->szd,proj->dim,NULL,&proj->XpV1));
      PetscCall(BVCreateVec(extop->X,&proj->w));
      PetscCall(BVDuplicateResize(extop->X,proj->dim,&proj->V1));
    }
    PetscCall(DSNEPSetComputeMatrixFunction(ds,NEPDeflationDSNEPComputeMatrix,(void*)proj));
  }

  /* Split Vext in V1 and V2 */
  if (extop->szd) {
    for (j=j0;j<j1;j++) {
      PetscCall(BVGetColumn(Vext,j,&ve));
      PetscCall(BVGetColumn(proj->V1,j,&v));
      PetscCall(NEPDeflationCopyToExtendedVec(extop,v,proj->V2+j*extop->szd,ve,PETSC_TRUE));
      PetscCall(BVRestoreColumn(proj->V1,j,&v));
      PetscCall(BVRestoreColumn(Vext,j,&ve));
    }
    V1 = proj->V1;
  } else V1 = Vext;

  /* Compute matrices V1^* A_i V1 */
  PetscCall(BVSetActiveColumns(V1,j0,j1));
  for (k=0;k<nep->nt;k++) {
    PetscCall(DSGetMat(ds,DSMatExtra[k],&G));
    PetscCall(BVMatProject(V1,nep->A[k],V1,G));
    PetscCall(DSRestoreMat(ds,DSMatExtra[k],&G));
  }

  if (extop->n) {
    if (extop->szd) {
      /* Compute matrices V1^* A_i X  and V1^* X */
      PetscCall(BVSetActiveColumns(extop->W,0,n));
      for (k=0;k<nep->nt;k++) {
        PetscCall(BVMatMult(extop->X,nep->A[k],extop->W));
        PetscCall(BVDot(extop->W,V1,proj->V1pApX[k]));
      }
      PetscCall(BVDot(V1,extop->X,proj->XpV1));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
