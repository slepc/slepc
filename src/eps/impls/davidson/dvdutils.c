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

   Some utils
*/

#include "davidson.h"

typedef struct {
  PC pc;
} dvdPCWrapper;

/*
  Configure the harmonics.
  switch (mode) {
  DVD_HARM_RR:    harmonic RR
  DVD_HARM_RRR:   relative harmonic RR
  DVD_HARM_REIGS: rightmost eigenvalues
  DVD_HARM_LEIGS: largest eigenvalues
  }
  fixedTarged, if true use the target instead of the best eigenvalue
  target, the fixed target to be used
*/
typedef struct {
  PetscScalar Wa,Wb;       /* span{W} = span{Wa*AV - Wb*BV} */
  PetscScalar Pa,Pb;       /* H=W'*(Pa*AV - Pb*BV), G=W'*(Wa*AV - Wb*BV) */
  PetscBool   withTarget;
  HarmType_t  mode;
} dvdHarmonic;

typedef struct {
  Vec diagA, diagB;
} dvdJacobiPrecond;

static PetscErrorCode dvd_improvex_precond_d(dvdDashboard *d)
{
  dvdPCWrapper   *dvdpc = (dvdPCWrapper*)d->improvex_precond_data;

  PetscFunctionBegin;
  /* Free local data */
  PetscCall(PCDestroy(&dvdpc->pc));
  PetscCall(PetscFree(d->improvex_precond_data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_static_precond_PC_0(dvdDashboard *d,PetscInt i,Vec x,Vec Px)
{
  dvdPCWrapper   *dvdpc = (dvdPCWrapper*)d->improvex_precond_data;

  PetscFunctionBegin;
  PetscCall(PCApply(dvdpc->pc,x,Px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Create a trivial preconditioner
*/
static PetscErrorCode dvd_precond_none(dvdDashboard *d,PetscInt i,Vec x,Vec Px)
{
  PetscFunctionBegin;
  PetscCall(VecCopy(x,Px));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Create a static preconditioner from a PC
*/
PetscErrorCode dvd_static_precond_PC(dvdDashboard *d,dvdBlackboard *b,PC pc)
{
  dvdPCWrapper   *dvdpc;
  Mat            P;
  PetscBool      t0,t1,t2;

  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    /* If the preconditioner is valid */
    if (pc) {
      PetscCall(PetscNew(&dvdpc));
      dvdpc->pc = pc;
      PetscCall(PetscObjectReference((PetscObject)pc));
      d->improvex_precond_data = dvdpc;
      d->improvex_precond = dvd_static_precond_PC_0;

      /* PC saves the matrix associated with the linear system, and it has to
         be initialize to a valid matrix */
      PetscCall(PCGetOperatorsSet(pc,NULL,&t0));
      PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCNONE,&t1));
      PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCSHELL,&t2));
      if (t0 && !t1) {
        PetscCall(PCGetOperators(pc,NULL,&P));
        PetscCall(PetscObjectReference((PetscObject)P));
        PetscCall(PCSetOperators(pc,P,P));
        PetscCall(PCSetReusePreconditioner(pc,PETSC_TRUE));
        PetscCall(MatDestroy(&P));
      } else if (t2) {
        PetscCall(PCSetOperators(pc,d->A,d->A));
        PetscCall(PCSetReusePreconditioner(pc,PETSC_TRUE));
      } else {
        d->improvex_precond = dvd_precond_none;
      }

      PetscCall(EPSDavidsonFLAdd(&d->destroyList,dvd_improvex_precond_d));

    /* Else, use no preconditioner */
    } else d->improvex_precond = dvd_precond_none;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_d(dvdDashboard *d)
{
  PetscFunctionBegin;
  /* Free local data */
  PetscCall(PetscFree(d->calcpairs_W_data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_transf(dvdHarmonic *dvdh,PetscScalar t)
{
  PetscFunctionBegin;
  switch (dvdh->mode) {
  case DVD_HARM_RR:    /* harmonic RR */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 0.0; dvdh->Pb = -1.0;
    break;
  case DVD_HARM_RRR:   /* relative harmonic RR */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 1.0; dvdh->Pb = 0.0;
    break;
  case DVD_HARM_REIGS: /* rightmost eigenvalues */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 1.0; dvdh->Pb = -PetscConj(t);
    break;
  case DVD_HARM_LEIGS: /* largest eigenvalues */
    dvdh->Wa = 0.0; dvdh->Wb = 1.0; dvdh->Pa = 1.0; dvdh->Pb = 0.0;
    break;
  case DVD_HARM_NONE:
  default:
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Harmonic type not supported");
  }

  /* Check the transformation does not change the sign of the imaginary part */
#if !PetscDefined(USE_COMPLEX)
  if (dvdh->Pb*dvdh->Wa - dvdh->Wb*dvdh->Pa < 0.0) {
    dvdh->Pa *= -1.0;
    dvdh->Pb *= -1.0;
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_updateW(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       l,k;
  BV             BX = d->BX?d->BX:d->eps->V;

  PetscFunctionBegin;
  /* Update the target if it is necessary */
  if (!data->withTarget) PetscCall(dvd_harm_transf(data,d->eigr[0]));

  /* W(i) <- Wa*AV(i) - Wb*BV(i) */
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  PetscAssert(k==l+d->V_new_s,PETSC_COMM_SELF,PETSC_ERR_PLIB,"Consistency broken");
  PetscCall(BVSetActiveColumns(d->W,l+d->V_new_s,l+d->V_new_e));
  PetscCall(BVSetActiveColumns(d->AX,l+d->V_new_s,l+d->V_new_e));
  PetscCall(BVSetActiveColumns(BX,l+d->V_new_s,l+d->V_new_e));
  PetscCall(BVCopy(d->AX,d->W));
  PetscCall(BVScale(d->W,data->Wa));
  PetscCall(BVMult(d->W,-data->Wb,1.0,BX,NULL));
  PetscCall(BVSetActiveColumns(d->W,l,k));
  PetscCall(BVSetActiveColumns(d->AX,l,k));
  PetscCall(BVSetActiveColumns(BX,l,k));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_proj(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,j,l0,l,k,ld;
  PetscScalar    h,g,*H,*G;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l0,&k));
  l = l0 + d->V_new_s;
  k = l0 + d->V_new_e;
  PetscCall(MatGetSize(d->H,&ld,NULL));
  PetscCall(MatDenseGetArray(d->H,&H));
  PetscCall(MatDenseGetArray(d->G,&G));
  /* [H G] <- [Pa*H - Pb*G, Wa*H - Wb*G] */
  /* Right part */
  for (i=l;i<k;i++) {
    for (j=l0;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  /* Left part */
  for (i=l0;i<l;i++) {
    for (j=l;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  PetscCall(MatDenseRestoreArray(d->H,&H));
  PetscCall(MatDenseRestoreArray(d->G,&G));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_harm_updateproj(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,j,l,k,ld;
  PetscScalar    h,g,*H,*G;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  k = l + d->V_tra_s;
  PetscCall(MatGetSize(d->H,&ld,NULL));
  PetscCall(MatDenseGetArray(d->H,&H));
  PetscCall(MatDenseGetArray(d->G,&G));
  /* [H G] <- [Pa*H - Pb*G, Wa*H - Wb*G] */
  /* Right part */
  for (i=l;i<k;i++) {
    for (j=0;j<l;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  /* Lower triangular part */
  for (i=0;i<l;i++) {
    for (j=l;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  PetscCall(MatDenseRestoreArray(d->H,&H));
  PetscCall(MatDenseRestoreArray(d->G,&G));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_backtrans(dvdHarmonic *data,PetscScalar *ar,PetscScalar *ai)
{
  PetscScalar xr;
#if !PetscDefined(USE_COMPLEX)
  PetscScalar xi, k;
#endif

  PetscFunctionBegin;
  xr = *ar;
#if !PetscDefined(USE_COMPLEX)
  xi = *ai;
  if (PetscUnlikely(xi != 0.0)) {
    k = (data->Pa - data->Wa*xr)*(data->Pa - data->Wa*xr) + data->Wa*data->Wa*xi*xi;
    *ar = (data->Pb*data->Pa - (data->Pb*data->Wa + data->Wb*data->Pa)*xr + data->Wb*data->Wa*(xr*xr + xi*xi))/k;
    *ai = (data->Pb*data->Wa - data->Wb*data->Pa)*xi/k;
  } else
#endif
    *ar = (data->Pb - data->Wb*xr) / (data->Pa - data->Wa*xr);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_eig_backtrans(dvdDashboard *d,PetscScalar ar,PetscScalar ai,PetscScalar *br,PetscScalar *bi)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;

  PetscFunctionBegin;
  PetscCall(dvd_harm_backtrans(data,&ar,&ai));
  *br = ar;
  *bi = ai;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode dvd_harm_eigs_trans(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,l,k;

  PetscFunctionBegin;
  PetscCall(BVGetActiveColumns(d->eps->V,&l,&k));
  for (i=0;i<k-l;i++) PetscCall(dvd_harm_backtrans(data,&d->eigr[i],&d->eigi[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvd_harm_conf(dvdDashboard *d,dvdBlackboard *b,HarmType_t mode,PetscBool fixedTarget,PetscScalar t)
{
  dvdHarmonic    *dvdh;

  PetscFunctionBegin;
  /* Set the problem to GNHEP:
     d->G maybe is upper triangular due to biorthogonality of V and W */
  d->sEP = d->sA = d->sB = 0;

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    PetscCall(PetscNew(&dvdh));
    dvdh->withTarget = fixedTarget;
    dvdh->mode = mode;
    if (fixedTarget) PetscCall(dvd_harm_transf(dvdh, t));
    d->calcpairs_W_data = dvdh;
    d->calcpairs_W = dvd_harm_updateW;
    d->calcpairs_proj_trans = dvd_harm_proj;
    d->calcpairs_eigs_trans = dvd_harm_eigs_trans;
    d->calcpairs_eig_backtrans = dvd_harm_eig_backtrans;

    PetscCall(EPSDavidsonFLAdd(&d->destroyList,dvd_harm_d));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
