/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: test for restarting, updateV, restartV
*/

#include "davidson.h"

typedef struct {
  PetscInt          min_size_V;        /* restart with this number of eigenvectors */
  PetscInt          plusk;             /* at restart, save plusk vectors from last iteration */
  PetscInt          mpd;               /* max size of the searching subspace */
  void              *old_updateV_data; /* old updateV data */
  PetscErrorCode    (*old_isRestarting)(dvdDashboard*,PetscBool*);  /* old isRestarting */
  Mat               oldU;              /* previous projected right igenvectors */
  Mat               oldV;              /* previous projected left eigenvectors */
  PetscInt          size_oldU;         /* size of oldU */
  PetscBool         allResiduals;      /* if computing all the residuals */
} dvdManagV_basic;

static PetscErrorCode dvd_updateV_start(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        i;

  PetscFunctionBegin;
  for (i=0;i<d->eps->ncv;i++) d->eigi[i] = 0.0;
  d->nR = d->real_nR;
  for (i=0;i<d->eps->ncv;i++) d->nR[i] = 1.0;
  d->nX = d->real_nX;
  for (i=0;i<d->eps->ncv;i++) d->errest[i] = 1.0;
  data->size_oldU = 0;
  d->nconv = 0;
  d->npreconv = 0;
  d->V_tra_s = d->V_tra_e = d->V_new_s = d->V_new_e = 0;
  d->size_D = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_isrestarting_fullV(dvdDashboard *d,PetscBool *r)
{
  PetscInt        l,k;
  PetscBool       restart;
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  restart = (k+2 > d->eps->ncv)? PETSC_TRUE: PETSC_FALSE;

  /* Check old isRestarting function */
  if (PetscUnlikely(!restart && data->old_isRestarting)) PetscCall(data->old_isRestarting(d,&restart));
  *r = restart;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_managementV_basic_d(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  /* Restore changes in dvdDashboard */
  d->updateV_data = data->old_updateV_data;

  /* Free local data */
  PetscCall(MatDestroy(&data->oldU));
  PetscCall(MatDestroy(&data->oldV));
  PetscCall(PetscFree(d->real_nR));
  PetscCall(PetscFree(d->real_nX));
  PetscCall(PetscFree(data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_updateV_conv_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        npreconv,cMT,cMTX,lV,kV,nV;
  Mat             Z,Z0,Q,Q0;
  PetscBool       t;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt        i;
#endif

  PetscFunctionBegin;
  npreconv = d->npreconv;
  /* Constrains the converged pairs to nev */
#if !defined(PETSC_USE_COMPLEX)
  /* Tries to maintain together conjugate eigenpairs */
  for (i=0; (i + (d->eigi[i]!=0.0?1:0) < npreconv) && (d->nconv + i < d->nev); i+= (d->eigi[i]!=0.0?2:1));
  npreconv = i;
#else
  npreconv = PetscMax(PetscMin(d->nev-d->nconv,npreconv),0);
#endif
  /* For GHEP without B-ortho, converge all of the requested pairs at once */
  PetscCall(PetscObjectTypeCompare((PetscObject)d->eps->ds,DSGHEP,&t));
  if (t && d->nconv+npreconv<d->nev) npreconv = 0;
  /* Quick exit */
  if (npreconv == 0) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  nV  = kV - lV;
  cMT = nV - npreconv;
  /* Harmonics restarts with right eigenvectors, and other with the left ones.
     If the problem is standard or hermitian, left and right vectors are the same */
  if (!(d->W||DVD_IS(d->sEP,DVD_EP_STD)||DVD_IS(d->sEP,DVD_EP_HERMITIAN))) {
    /* ps.Q <- [ps.Q(0:npreconv-1) ps.Z(npreconv:size_H-1)] */
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Z,&Z));
    PetscCall(MatDenseGetSubMatrix(Q,0,npreconv,nV,npreconv+cMT,&Q0));
    PetscCall(MatDenseGetSubMatrix(Z,0,npreconv,nV,npreconv+cMT,&Z0));
    PetscCall(MatCopy(Z0,Q0,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(Q,&Q0));
    PetscCall(MatDenseRestoreSubMatrix(Z,&Z0));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Z,&Z));
  }
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) PetscCall(DSPseudoOrthogonalize(d->eps->ds,DS_MAT_Q,nV,d->nBds,&cMTX,d->nBds));
  else PetscCall(DSOrthogonalize(d->eps->ds,DS_MAT_Q,nV,&cMTX));
  cMT = cMTX - npreconv;

  if (d->W) {
    PetscCall(DSOrthogonalize(d->eps->ds,DS_MAT_Z,nV,&cMTX));
    cMT = PetscMin(cMT,cMTX - npreconv);
  }

  /* Lock the converged pairs */
  d->eigr+= npreconv;
#if !defined(PETSC_USE_COMPLEX)
  if (d->eigi) d->eigi+= npreconv;
#endif
  d->nconv+= npreconv;
  d->errest+= npreconv;
  /* Notify the changes in V and update the other subspaces */
  d->V_tra_s = npreconv;          d->V_tra_e = nV;
  d->V_new_s = cMT;               d->V_new_e = d->V_new_s;
  /* Remove oldU */
  data->size_oldU = 0;

  d->npreconv-= npreconv;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_updateV_restart_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        lV,kV,nV,size_plusk,size_X,cMTX,cMTY,max_restart_size;
  Mat             Q,Q0,Z,Z0,U,V;

  PetscFunctionBegin;
  /* Select size_X desired pairs from V */
  /* The restarted basis should:
     - have at least one spot to add a new direction;
     - keep converged vectors, npreconv;
     - keep at least 1 oldU direction if possible.
  */
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  nV = kV - lV;
  max_restart_size = PetscMax(0,PetscMin(d->eps->mpd - 1,d->eps->ncv - lV - 2));
  size_X = PetscMin(PetscMin(data->min_size_V+d->npreconv,max_restart_size - (max_restart_size - d->npreconv > 1 && data->plusk > 0 && data->size_oldU > 0 ? 1 : 0)), nV);

  /* Add plusk eigenvectors from the previous iteration */
  size_plusk = PetscMax(0,PetscMin(PetscMin(PetscMin(data->plusk,data->size_oldU),max_restart_size - size_X),nV - size_X));

  d->size_MT = nV;
  /* ps.Q <- orth([pX(0:size_X-1) [oldU(0:size_plusk-1); 0] ]) */
  /* Harmonics restarts with right eigenvectors, and other with the left ones.
     If the problem is standard or hermitian, left and right vectors are the same */
  if (!(d->W||DVD_IS(d->sEP,DVD_EP_STD)||DVD_IS(d->sEP,DVD_EP_HERMITIAN))) {
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Z,&Z));
    PetscCall(MatDenseGetSubMatrix(Q,0,nV,0,size_X,&Q0));
    PetscCall(MatDenseGetSubMatrix(Z,0,nV,0,size_X,&Z0));
    PetscCall(MatCopy(Z0,Q0,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(Q,&Q0));
    PetscCall(MatDenseRestoreSubMatrix(Z,&Z0));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Z,&Z));
  }
  PetscCheck(size_plusk<=0 || !DVD_IS(d->sEP,DVD_EP_INDEFINITE),PETSC_COMM_SELF,PETSC_ERR_SUP,"Unsupported plusk>0 in indefinite eigenvalue problems");
  if (size_plusk > 0) {
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(MatDenseGetSubMatrix(Q,0,nV,size_X,size_X+size_plusk,&Q0));
    PetscCall(MatDenseGetSubMatrix(data->oldU,0,nV,0,size_plusk,&U));
    PetscCall(MatCopy(U,Q0,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(Q,&Q0));
    PetscCall(MatDenseRestoreSubMatrix(data->oldU,&U));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Q,&Q));
  }
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) PetscCall(DSPseudoOrthogonalize(d->eps->ds,DS_MAT_Q,size_X,d->nBds,&cMTX,d->nBds));
  else PetscCall(DSOrthogonalize(d->eps->ds,DS_MAT_Q,size_X+size_plusk,&cMTX));

  if (d->W && size_plusk > 0) {
    /* ps.Z <- orth([ps.Z(0:size_X-1) [oldV(0:size_plusk-1); 0] ]) */
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Z,&Z));
    PetscCall(MatDenseGetSubMatrix(Z,0,nV,size_X,size_X+size_plusk,&Z0));
    PetscCall(MatDenseGetSubMatrix(data->oldV,0,nV,0,size_plusk,&V));
    PetscCall(MatCopy(V,Z0,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(Z,&Z0));
    PetscCall(MatDenseRestoreSubMatrix(data->oldV,&V));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Z,&Z));
    PetscCall(DSOrthogonalize(d->eps->ds,DS_MAT_Z,size_X+size_plusk,&cMTY));
    cMTX = PetscMin(cMTX, cMTY);
  }
  PetscAssert(cMTX<=size_X+size_plusk,PETSC_COMM_SELF,PETSC_ERR_SUP,"Invalid number of columns to restart");

  /* Notify the changes in V and update the other subspaces */
  d->V_tra_s = 0;                     d->V_tra_e = cMTX;
  d->V_new_s = d->V_tra_e;            d->V_new_e = d->V_new_s;

  /* Remove oldU */
  data->size_oldU = 0;

  /* Remove npreconv */
  d->npreconv = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_updateV_testConv(dvdDashboard *d,PetscInt s,PetscInt pre,PetscInt e,PetscInt *nConv)
{
  PetscInt        i,j,b;
  PetscReal       norm;
  PetscBool       conv, c;
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  if (nConv) *nConv = s;
  for (i=s,conv=PETSC_TRUE;(conv || data->allResiduals) && (i < e);i+=b) {
#if !defined(PETSC_USE_COMPLEX)
    b = d->eigi[i]!=0.0?2:1;
#else
    b = 1;
#endif
    if (i+b-1 >= pre) PetscCall(d->calcpairs_residual(d,i,i+b));
    /* Test the Schur vector */
    for (j=0,c=PETSC_TRUE;j<b && c;j++) {
      norm = d->nR[i+j]/d->nX[i+j];
      c = d->testConv(d,d->eigr[i+j],d->eigi[i+j],norm,&d->errest[i+j]);
    }
    if (conv && c) {
      if (nConv) *nConv = i+b;
    } else conv = PETSC_FALSE;
  }
  pre = PetscMax(pre,i);

#if !defined(PETSC_USE_COMPLEX)
  /* Enforce converged conjugate complex eigenpairs */
  if (nConv) {
    for (j=0;j<*nConv;j++) if (d->eigi[j] != 0.0) j++;
    if (j>*nConv) (*nConv)--;
  }
#endif
  for (i=pre;i<e;i++) d->errest[i] = d->nR[i] = 1.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_updateV_update_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        size_D,s,lV,kV,nV;
  Mat             Q,Q0,Z,Z0,U,V;

  PetscFunctionBegin;
  /* Select the desired pairs */
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  nV = kV - lV;
  size_D = PetscMin(PetscMin(PetscMin(d->bs,nV),d->eps->ncv-nV),nV);
  if (size_D == 0) PetscFunctionReturn(PETSC_SUCCESS);

  /* Fill V with D */
  PetscCall(d->improveX(d,d->npreconv,d->npreconv+size_D,&size_D));

  /* If D is empty, exit */
  d->size_D = size_D;
  if (size_D == 0) PetscFunctionReturn(PETSC_SUCCESS);

  /* Get the residual of all pairs */
#if !defined(PETSC_USE_COMPLEX)
  s = (d->eigi[0]!=0.0)? 2: 1;
#else
  s = 1;
#endif
  PetscCall(BVGetActiveColumns(d->eps->V,&lV,&kV));
  nV = kV - lV;
  PetscCall(dvd_updateV_testConv(d,s,s,data->allResiduals?nV:size_D,NULL));

  /* Notify the changes in V */
  d->V_tra_s = 0;                 d->V_tra_e = 0;
  d->V_new_s = nV;                d->V_new_e = nV+size_D;

  /* Save the projected eigenvectors */
  if (data->plusk > 0) {
    PetscCall(MatZeroEntries(data->oldU));
    data->size_oldU = nV;
    PetscCall(DSGetMat(d->eps->ds,DS_MAT_Q,&Q));
    PetscCall(MatDenseGetSubMatrix(Q,0,nV,0,nV,&Q0));
    PetscCall(MatDenseGetSubMatrix(data->oldU,0,nV,0,nV,&U));
    PetscCall(MatCopy(Q0,U,SAME_NONZERO_PATTERN));
    PetscCall(MatDenseRestoreSubMatrix(Q,&Q0));
    PetscCall(MatDenseRestoreSubMatrix(data->oldU,&U));
    PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Q,&Q));
    if (d->W) {
      PetscCall(MatZeroEntries(data->oldV));
      PetscCall(DSGetMat(d->eps->ds,DS_MAT_Z,&Z));
      PetscCall(MatDenseGetSubMatrix(Z,0,nV,0,nV,&Z0));
      PetscCall(MatDenseGetSubMatrix(data->oldV,0,nV,0,nV,&V));
      PetscCall(MatCopy(Z0,V,SAME_NONZERO_PATTERN));
      PetscCall(MatDenseRestoreSubMatrix(Z,&Z0));
      PetscCall(MatDenseRestoreSubMatrix(data->oldV,&V));
      PetscCall(DSRestoreMat(d->eps->ds,DS_MAT_Z,&Z));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_updateV_extrapol(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        i;
  PetscBool       restart,t;

  PetscFunctionBegin;
  /* TODO: restrict select pairs to each case */
  PetscCall(d->calcpairs_selectPairs(d, data->min_size_V+d->npreconv));

  /* If the subspaces doesn't need restart, add new vector */
  PetscCall(d->isRestarting(d,&restart));
  if (!restart) {
    d->size_D = 0;
    PetscCall(dvd_updateV_update_gen(d));

    /* If no vector were converged, exit */
    /* For GHEP without B-ortho, converge all of the requested pairs at once */
    PetscCall(PetscObjectTypeCompare((PetscObject)d->eps->ds,DSGHEP,&t));
    if (d->nconv+d->npreconv < d->nev && (t || d->npreconv == 0)) PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* If some eigenpairs were converged, lock them  */
  if (d->npreconv > 0) {
    i = d->npreconv;
    PetscCall(dvd_updateV_conv_gen(d));

    /* If some eigenpair was locked, exit */
    if (i > d->npreconv) PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* Else, a restarting is performed */
  PetscCall(dvd_updateV_restart_gen(d));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_managementV_basic(dvdDashboard *d,dvdBlackboard *b,PetscInt bs,PetscInt mpd,PetscInt min_size_V,PetscInt plusk,PetscBool harm,PetscBool allResiduals)
{
  dvdManagV_basic *data;
#if !defined(PETSC_USE_COMPLEX)
  PetscBool       her_probl,std_probl;
#endif

  PetscFunctionBegin;
  /* Setting configuration constrains */
#if !defined(PETSC_USE_COMPLEX)
  /* if the last converged eigenvalue is complex its conjugate pair is also
     converged */
  her_probl = DVD_IS(d->sEP,DVD_EP_HERMITIAN)? PETSC_TRUE: PETSC_FALSE;
  std_probl = DVD_IS(d->sEP,DVD_EP_STD)? PETSC_TRUE: PETSC_FALSE;
  b->max_size_X = PetscMax(b->max_size_X,bs+((her_probl && std_probl)?0:1));
#else
  b->max_size_X = PetscMax(b->max_size_X,bs);
#endif

  b->max_size_V = PetscMax(b->max_size_V,mpd);
  min_size_V = PetscMin(min_size_V,mpd-bs);
  b->size_V = PetscMax(b->size_V,b->max_size_V+b->max_size_P+b->max_nev);
  b->max_size_oldX = plusk;

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    PetscCall(PetscNew(&data));
    data->mpd = b->max_size_V;
    data->min_size_V = min_size_V;
    d->bs = bs;
    data->plusk = plusk;
    data->allResiduals = allResiduals;

    d->eigr = d->eps->eigr;
    d->eigi = d->eps->eigi;
    d->errest = d->eps->errest;
    PetscCall(PetscMalloc1(d->eps->ncv,&d->real_nR));
    PetscCall(PetscMalloc1(d->eps->ncv,&d->real_nX));
    if (plusk > 0) PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&data->oldU));
    else data->oldU = NULL;
    if (harm && plusk>0) PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&data->oldV));
    else data->oldV = NULL;

    data->old_updateV_data = d->updateV_data;
    d->updateV_data = data;
    data->old_isRestarting = d->isRestarting;
    d->isRestarting = dvd_isrestarting_fullV;
    d->updateV = dvd_updateV_extrapol;
    d->preTestConv = dvd_updateV_testConv;
    PetscCall(EPSDavidsonFLAdd(&d->startList,dvd_updateV_start));
    PetscCall(EPSDavidsonFLAdd(&d->destroyList,dvd_managementV_basic_d));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
