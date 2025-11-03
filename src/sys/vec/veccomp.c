/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/vecimplslepc.h>     /*I "slepcvec.h" I*/

/* Private MPI datatypes and operators */
static MPI_Datatype MPIU_NORM2=0, MPIU_NORM1_AND_2=0;
static PetscBool VecCompInitialized = PETSC_FALSE;
MPI_Op MPIU_NORM2_SUM=0;

/* Private functions */
static inline void SumNorm2(PetscReal*,PetscReal*,PetscReal*,PetscReal*);
static inline PetscReal GetNorm2(PetscReal,PetscReal);
static inline void AddNorm2(PetscReal*,PetscReal*,PetscReal);
static PetscErrorCode VecCompSetSubVecs_Comp(Vec,PetscInt,Vec*);
static PetscErrorCode VecCompGetSubVecs_Comp(Vec,PetscInt*,const Vec**);

#include "veccomp0.h"

#define __WITH_MPI__
#include "veccomp0.h"

static inline void SumNorm2(PetscReal *ssq0,PetscReal *scale0,PetscReal *ssq1,PetscReal *scale1)
{
  PetscReal q;
  if (*scale0 > *scale1) {
    q = *scale1/(*scale0);
    *ssq1 = *ssq0 + q*q*(*ssq1);
    *scale1 = *scale0;
  } else {
    q = *scale0/(*scale1);
    *ssq1 += q*q*(*ssq0);
  }
}

static inline PetscReal GetNorm2(PetscReal ssq,PetscReal scale)
{
  return scale*PetscSqrtReal(ssq);
}

static inline void AddNorm2(PetscReal *ssq,PetscReal *scale,PetscReal x)
{
  PetscReal absx,q;
  if (x != 0.0) {
    absx = PetscAbs(x);
    if (*scale < absx) {
      q = *scale/absx;
      *ssq = 1.0 + *ssq*q*q;
      *scale = absx;
    } else {
      q = absx/(*scale);
      *ssq += q*q;
    }
  }
}

SLEPC_EXTERN void MPIAPI SlepcSumNorm2_Local(void *in,void *out,PetscMPIInt *cnt,MPI_Datatype *datatype)
{
  PetscInt  i,count = *cnt;
  PetscReal *xin = (PetscReal*)in,*xout = (PetscReal*)out;

  PetscFunctionBegin;
  if (*datatype == MPIU_NORM2) {
    for (i=0;i<count;i++) {
      SumNorm2(&xin[i*2],&xin[i*2+1],&xout[i*2],&xout[i*2+1]);
    }
  } else if (*datatype == MPIU_NORM1_AND_2) {
    for (i=0;i<count;i++) {
      xout[i*3] += xin[i*3];
      SumNorm2(&xin[i*3+1],&xin[i*3+2],&xout[i*3+1],&xout[i*3+2]);
    }
  } else {
    (void)(*PetscErrorPrintf)("Can only handle MPIU_NORM* data types");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  PetscFunctionReturnVoid();
}

static PetscErrorCode VecCompNormEnd(void)
{
  PetscFunctionBegin;
  PetscCallMPI(MPI_Type_free(&MPIU_NORM2));
  PetscCallMPI(MPI_Type_free(&MPIU_NORM1_AND_2));
  PetscCallMPI(MPI_Op_free(&MPIU_NORM2_SUM));
  VecCompInitialized = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode VecCompNormInit(void)
{
  PetscFunctionBegin;
  PetscCallMPI(MPI_Type_contiguous(2,MPIU_REAL,&MPIU_NORM2));
  PetscCallMPI(MPI_Type_commit(&MPIU_NORM2));
  PetscCallMPI(MPI_Type_contiguous(3,MPIU_REAL,&MPIU_NORM1_AND_2));
  PetscCallMPI(MPI_Type_commit(&MPIU_NORM1_AND_2));
  PetscCallMPI(MPI_Op_create(SlepcSumNorm2_Local,PETSC_TRUE,&MPIU_NORM2_SUM));
  PetscCall(PetscRegisterFinalize(VecCompNormEnd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecDestroy_Comp(Vec v)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       i;

  PetscFunctionBegin;
#if defined(PETSC_USE_LOG)
  PetscCall(PetscLogObjectState((PetscObject)v,"Length=%" PetscInt_FMT,v->map->n));
#endif
  for (i=0;i<vs->nx;i++) PetscCall(VecDestroy(&vs->x[i]));
  if (--vs->n->friends <= 0) PetscCall(PetscFree(vs->n));
  PetscCall(PetscFree(vs->x));
  PetscCall(PetscFree(vs));
  PetscCall(PetscObjectComposeFunction((PetscObject)v,"VecCompSetSubVecs_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)v,"VecCompGetSubVecs_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static struct _VecOps DvOps = {
  PetscDesignatedInitializer(duplicate,VecDuplicate_Comp),
  PetscDesignatedInitializer(duplicatevecs,VecDuplicateVecs_Comp),
  PetscDesignatedInitializer(destroyvecs,VecDestroyVecs_Comp),
  PetscDesignatedInitializer(dot,VecDot_Comp_MPI),
  PetscDesignatedInitializer(mdot,VecMDot_Comp_MPI),
  PetscDesignatedInitializer(norm,VecNorm_Comp_MPI),
  PetscDesignatedInitializer(tdot,VecTDot_Comp_MPI),
  PetscDesignatedInitializer(mtdot,VecMTDot_Comp_MPI),
  PetscDesignatedInitializer(scale,VecScale_Comp),
  PetscDesignatedInitializer(copy,VecCopy_Comp),
  PetscDesignatedInitializer(set,VecSet_Comp),
  PetscDesignatedInitializer(swap,VecSwap_Comp),
  PetscDesignatedInitializer(axpy,VecAXPY_Comp),
  PetscDesignatedInitializer(axpby,VecAXPBY_Comp),
  PetscDesignatedInitializer(maxpy,VecMAXPY_Comp),
  PetscDesignatedInitializer(aypx,VecAYPX_Comp),
  PetscDesignatedInitializer(waxpy,VecWAXPY_Comp),
  PetscDesignatedInitializer(axpbypcz,VecAXPBYPCZ_Comp),
  PetscDesignatedInitializer(pointwisemult,VecPointwiseMult_Comp),
  PetscDesignatedInitializer(pointwisedivide,VecPointwiseDivide_Comp),
  PetscDesignatedInitializer(setvalues,NULL),
  PetscDesignatedInitializer(assemblybegin,NULL),
  PetscDesignatedInitializer(assemblyend,NULL),
  PetscDesignatedInitializer(getarray,NULL),
  PetscDesignatedInitializer(getsize,VecGetSize_Comp),
  PetscDesignatedInitializer(getlocalsize,VecGetLocalSize_Comp),
  PetscDesignatedInitializer(restorearray,NULL),
  PetscDesignatedInitializer(max,VecMax_Comp),
  PetscDesignatedInitializer(min,VecMin_Comp),
  PetscDesignatedInitializer(setrandom,VecSetRandom_Comp),
  PetscDesignatedInitializer(setoption,NULL),
  PetscDesignatedInitializer(setvaluesblocked,NULL),
  PetscDesignatedInitializer(destroy,VecDestroy_Comp),
  PetscDesignatedInitializer(view,VecView_Comp),
  PetscDesignatedInitializer(placearray,NULL),
  PetscDesignatedInitializer(replacearray,NULL),
  PetscDesignatedInitializer(dot_local,VecDot_Comp_Seq),
  PetscDesignatedInitializer(tdot_local,VecTDot_Comp_Seq),
  PetscDesignatedInitializer(norm_local,VecNorm_Comp_Seq),
  PetscDesignatedInitializer(mdot_local,VecMDot_Comp_Seq),
  PetscDesignatedInitializer(mtdot_local,VecMTDot_Comp_Seq),
  PetscDesignatedInitializer(load,NULL),
  PetscDesignatedInitializer(reciprocal,VecReciprocal_Comp),
  PetscDesignatedInitializer(conjugate,VecConjugate_Comp),
  PetscDesignatedInitializer(setlocaltoglobalmapping,NULL),
  PetscDesignatedInitializer(getlocaltoglobalmapping,NULL),
  PetscDesignatedInitializer(resetarray,NULL),
  PetscDesignatedInitializer(setfromoptions,NULL),
  PetscDesignatedInitializer(maxpointwisedivide,VecMaxPointwiseDivide_Comp),
  PetscDesignatedInitializer(pointwisemax,VecPointwiseMax_Comp),
  PetscDesignatedInitializer(pointwisemaxabs,VecPointwiseMaxAbs_Comp),
  PetscDesignatedInitializer(pointwisemin,VecPointwiseMin_Comp),
  PetscDesignatedInitializer(getvalues,NULL),
  PetscDesignatedInitializer(sqrt,VecSqrtAbs_Comp),
  PetscDesignatedInitializer(abs,VecAbs_Comp),
  PetscDesignatedInitializer(exp,VecExp_Comp),
  PetscDesignatedInitializer(log,VecLog_Comp),
  PetscDesignatedInitializer(shift,VecShift_Comp),
  PetscDesignatedInitializer(create,NULL),
  PetscDesignatedInitializer(stridegather,NULL),
  PetscDesignatedInitializer(stridescatter,NULL),
  PetscDesignatedInitializer(dotnorm2,VecDotNorm2_Comp_MPI),
  PetscDesignatedInitializer(getsubvector,NULL),
  PetscDesignatedInitializer(restoresubvector,NULL),
  PetscDesignatedInitializer(getarrayread,NULL),
  PetscDesignatedInitializer(restorearrayread,NULL),
  PetscDesignatedInitializer(stridesubsetgather,NULL),
  PetscDesignatedInitializer(stridesubsetscatter,NULL),
  PetscDesignatedInitializer(viewnative,NULL),
  PetscDesignatedInitializer(loadnative,NULL),
  PetscDesignatedInitializer(getlocalvector,NULL)
};

PetscErrorCode VecDuplicateVecs_Comp(Vec w,PetscInt m,Vec *V[])
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(w,VEC_CLASSID,1);
  PetscAssertPointer(V,3);
  PetscCheck(m>0,PetscObjectComm((PetscObject)w),PETSC_ERR_ARG_OUTOFRANGE,"m must be > 0: m = %" PetscInt_FMT,m);
  PetscCall(PetscMalloc1(m,V));
  for (i=0;i<m;i++) PetscCall(VecDuplicate(w,*V+i));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecDestroyVecs_Comp(PetscInt m,Vec v[])
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscAssertPointer(v,2);
  PetscCheck(m>0,PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"m must be > 0: m = %" PetscInt_FMT,m);
  for (i=0;i<m;i++) PetscCall(VecDestroy(&v[i]));
  PetscCall(PetscFree(v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode VecCreate_Comp_Private(Vec v,Vec *x,PetscInt nx,PetscBool x_to_me,Vec_Comp_N *n)
{
  Vec_Comp       *s;
  PetscInt       N=0,lN=0,i,k;

  PetscFunctionBegin;
  if (!VecCompInitialized) {
    VecCompInitialized = PETSC_TRUE;
    PetscCall(VecRegister(VECCOMP,VecCreate_Comp));
    PetscCall(VecCompNormInit());
  }

  /* Allocate a new Vec_Comp */
  if (v->data) PetscCall(PetscFree(v->data));
  PetscCall(PetscNew(&s));
  PetscCall(PetscMemcpy(v->ops,&DvOps,sizeof(DvOps)));
  v->data        = (void*)s;
  v->petscnative = PETSC_FALSE;

  /* Allocate the array of Vec, if it is needed to be done */
  if (!x_to_me) {
    if (nx) PetscCall(PetscMalloc1(nx,&s->x));
    if (x) PetscCall(PetscArraycpy(s->x,x,nx));
  } else s->x = x;

  s->nx = nx;

  if (nx && x) {
    /* Allocate the shared structure, if it is not given */
    if (!n) {
      for (i=0;i<nx;i++) {
        PetscCall(VecGetSize(x[i],&k));
        N+= k;
        PetscCall(VecGetLocalSize(x[i],&k));
        lN+= k;
      }
      PetscCall(PetscNew(&n));
      s->n = n;
      n->n = nx;
      n->N = N;
      n->lN = lN;
      n->friends = 1;
    } else { /* If not, check in the vector in the shared structure */
      s->n = n;
      s->n->friends++;
    }

    /* Set the virtual sizes as the real sizes of the vector */
    PetscCall(VecSetSizes(v,s->n->lN,s->n->N));
  }

  PetscCall(PetscObjectChangeTypeName((PetscObject)v,VECCOMP));
  PetscCall(PetscObjectComposeFunction((PetscObject)v,"VecCompSetSubVecs_C",VecCompSetSubVecs_Comp));
  PetscCall(PetscObjectComposeFunction((PetscObject)v,"VecCompGetSubVecs_C",VecCompGetSubVecs_Comp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode VecCreate_Comp(Vec V)
{
  PetscFunctionBegin;
  PetscCall(VecCreate_Comp_Private(V,NULL,0,PETSC_FALSE,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   VECCOMP - VECCOMP = "comp" - Vector type consisting of several subvectors,
   each stored separately.

   Level: developer

   Notes:
   This is similar to PETSc's `VECNEST` but customized for SLEPc's needs. In particular,
   the number of child vectors can be modified dynamically, with `VecCompSetSubVecs()`.

.seealso: `Vec`, `VecType`, `VecCreateComp()`, `VecCreateCompWithVecs()`
M*/

/*@
   VecCreateComp - Creates a new vector containing several subvectors,
   each stored separately.

   Collective

   Input Parameters:
+  comm - communicator for the new `Vec`
.  Nx   - array of (initial) global sizes of child vectors
.  n    - number of child vectors
.  t    - type of the child vectors
-  Vparent - (optional) template vector

   Output Parameter:
.  V - new vector

   Notes:
   This is similar to PETSc's `VECNEST` but customized for SLEPc's needs. In particular,
   the number of child vectors can be modified dynamically, with `VecCompSetSubVecs()`.

   Level: developer

.seealso: `VecCreateCompWithVecs()`, `VecCompSetSubVecs()`
@*/
PetscErrorCode VecCreateComp(MPI_Comm comm,PetscInt Nx[],PetscInt n,VecType t,Vec Vparent,Vec *V)
{
  Vec            *x;
  PetscInt       i;

  PetscFunctionBegin;
  PetscCall(VecCreate(comm,V));
  PetscCall(PetscMalloc1(n,&x));
  for (i=0;i<n;i++) {
    PetscCall(VecCreate(comm,&x[i]));
    PetscCall(VecSetSizes(x[i],PETSC_DECIDE,Nx[i]));
    PetscCall(VecSetType(x[i],t));
  }
  PetscCall(VecCreate_Comp_Private(*V,x,n,PETSC_TRUE,Vparent?((Vec_Comp*)Vparent->data)->n:NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecCreateCompWithVecs - Creates a new vector containing several subvectors,
   each stored separately, from an array of `Vec`s.

   Collective

   Input Parameters:
+  x - array of `Vec`s
.  n - number of child vectors
-  Vparent - (optional) template vector

   Output Parameter:
.  V - new vector

   Level: developer

.seealso: `VecCreateComp()`
@*/
PetscErrorCode VecCreateCompWithVecs(Vec x[],PetscInt n,Vec Vparent,Vec *V)
{
  PetscInt       i;

  PetscFunctionBegin;
  PetscAssertPointer(x,1);
  PetscValidHeaderSpecific(*x,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(*x,n,2);
  PetscCall(VecCreate(PetscObjectComm((PetscObject)x[0]),V));
  for (i=0;i<n;i++) PetscCall(PetscObjectReference((PetscObject)x[i]));
  PetscCall(VecCreate_Comp_Private(*V,x,n,PETSC_FALSE,Vparent?((Vec_Comp*)Vparent->data)->n:NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecDuplicate_Comp(Vec win,Vec *V)
{
  Vec            *x;
  PetscInt       i;
  Vec_Comp       *s = (Vec_Comp*)win->data;

  PetscFunctionBegin;
  SlepcValidVecComp(win,1);
  PetscCall(VecCreate(PetscObjectComm((PetscObject)win),V));
  PetscCall(PetscMalloc1(s->nx,&x));
  for (i=0;i<s->nx;i++) {
    if (s->x[i]) PetscCall(VecDuplicate(s->x[i],&x[i]));
    else x[i] = NULL;
  }
  PetscCall(VecCreate_Comp_Private(*V,x,s->nx,PETSC_TRUE,s->n));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode VecCompGetSubVecs_Comp(Vec win,PetscInt *n,const Vec **x)
{
  Vec_Comp *s = (Vec_Comp*)win->data;

  PetscFunctionBegin;
  if (x) *x = s->x;
  if (n) *n = s->n->n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
   VecCompGetSubVecs - Returns the entire array of vectors defining a
   compound vector.

   Collective

   Input Parameter:
.  win - compound vector

   Output Parameters:
+  n - number of child vectors
-  x - array of child vectors

   Level: developer

.seealso: `VecCreateComp()`
@*/
PetscErrorCode VecCompGetSubVecs(Vec win,PetscInt *n,const Vec *x[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(win,VEC_CLASSID,1);
  PetscUseMethod(win,"VecCompGetSubVecs_C",(Vec,PetscInt*,const Vec**),(win,n,x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode VecCompSetSubVecs_Comp(Vec win,PetscInt n,Vec *x)
{
  Vec_Comp       *s = (Vec_Comp*)win->data;
  PetscInt       i,N,nlocal;
  Vec_Comp_N     *nn;

  PetscFunctionBegin;
  PetscCheck(s,PetscObjectComm((PetscObject)win),PETSC_ERR_ORDER,"Must call VecSetSizes first");
  if (!s->nx) {
    /* vector has been created via VecCreate+VecSetType+VecSetSizes, so allocate data structures */
    PetscCall(PetscMalloc1(n,&s->x));
    PetscCall(VecGetSize(win,&N));
    PetscCheck(N%n==0,PetscObjectComm((PetscObject)win),PETSC_ERR_SUP,"Global dimension %" PetscInt_FMT " is not divisible by %" PetscInt_FMT,N,n);
    PetscCall(VecGetLocalSize(win,&nlocal));
    PetscCheck(nlocal%n==0,PetscObjectComm((PetscObject)win),PETSC_ERR_SUP,"Local dimension %" PetscInt_FMT " is not divisible by %" PetscInt_FMT,nlocal,n);
    s->nx = n;
    for (i=0;i<n;i++) {
      PetscCall(VecCreate(PetscObjectComm((PetscObject)win),&s->x[i]));
      PetscCall(VecSetSizes(s->x[i],nlocal/n,N/n));
      PetscCall(VecSetFromOptions(s->x[i]));
    }
    if (!s->n) {
      PetscCall(PetscNew(&nn));
      s->n = nn;
      nn->N = N;
      nn->lN = nlocal;
      nn->friends = 1;
    }
  } else PetscCheck(n<=s->nx,PetscObjectComm((PetscObject)win),PETSC_ERR_SUP,"Number of child vectors cannot be larger than %" PetscInt_FMT,s->nx);
  if (x) PetscCall(PetscArraycpy(s->x,x,n));
  s->n->n = n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecCompSetSubVecs - Resets the number of subvectors defining a compound vector,
   or replaces the subvectors.

   Collective

   Input Parameters:
+  win - compound vector
.  n - number of child vectors
-  x - array of child vectors

   Note:
   It is not possible to increase the number of subvectors with respect to the
   number set at its creation.

   Level: developer

.seealso: `VecCreateComp()`, `VecCompGetSubVecs()`
@*/
PetscErrorCode VecCompSetSubVecs(Vec win,PetscInt n,Vec x[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(win,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(win,n,2);
  PetscTryMethod(win,"VecCompSetSubVecs_C",(Vec,PetscInt,Vec*),(win,n,x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecAXPY_Comp(Vec v,PetscScalar alpha,Vec w)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  for (i=0;i<vs->n->n;i++) PetscCall(VecAXPY(vs->x[i],alpha,ws->x[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecAYPX_Comp(Vec v,PetscScalar alpha,Vec w)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  for (i=0;i<vs->n->n;i++) PetscCall(VecAYPX(vs->x[i],alpha,ws->x[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecAXPBY_Comp(Vec v,PetscScalar alpha,PetscScalar beta,Vec w)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,4);
  for (i=0;i<vs->n->n;i++) PetscCall(VecAXPBY(vs->x[i],alpha,beta,ws->x[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecMAXPY_Comp(Vec v,PetscInt n,const PetscScalar *alpha,Vec *w)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  Vec            *wx;
  PetscInt       i,j;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  for (i=0;i<n;i++) SlepcValidVecComp(w[i],4);

  PetscCall(PetscMalloc1(n,&wx));

  for (j=0;j<vs->n->n;j++) {
    for (i=0;i<n;i++) wx[i] = ((Vec_Comp*)w[i]->data)->x[j];
    PetscCall(VecMAXPY(vs->x[j],n,alpha,wx));
  }

  PetscCall(PetscFree(wx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecWAXPY_Comp(Vec v,PetscScalar alpha,Vec w,Vec z)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data,*zs = (Vec_Comp*)z->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  SlepcValidVecComp(z,4);
  for (i=0;i<vs->n->n;i++) PetscCall(VecWAXPY(vs->x[i],alpha,ws->x[i],zs->x[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecAXPBYPCZ_Comp(Vec v,PetscScalar alpha,PetscScalar beta,PetscScalar gamma,Vec w,Vec z)
{
  Vec_Comp        *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data,*zs = (Vec_Comp*)z->data;
  PetscInt        i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,5);
  SlepcValidVecComp(z,6);
  for (i=0;i<vs->n->n;i++) PetscCall(VecAXPBYPCZ(vs->x[i],alpha,beta,gamma,ws->x[i],zs->x[i]));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecGetSize_Comp(Vec v,PetscInt *size)
{
  Vec_Comp *vs = (Vec_Comp*)v->data;

  PetscFunctionBegin;
  PetscAssertPointer(size,2);
  if (vs->n) {
    SlepcValidVecComp(v,1);
    *size = vs->n->N;
  } else *size = v->map->N;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecGetLocalSize_Comp(Vec v,PetscInt *size)
{
  Vec_Comp *vs = (Vec_Comp*)v->data;

  PetscFunctionBegin;
  PetscAssertPointer(size,2);
  if (vs->n) {
    SlepcValidVecComp(v,1);
    *size = vs->n->lN;
  } else *size = v->map->n;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecMax_Comp(Vec v,PetscInt *idx,PetscReal *z)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       idxp,s=0,s0;
  PetscReal      zp,z0;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  if (!idx && !z) PetscFunctionReturn(PETSC_SUCCESS);

  if (vs->n->n > 0) PetscCall(VecMax(vs->x[0],idx?&idxp:NULL,&zp));
  else {
    zp = PETSC_MIN_REAL;
    if (idx) idxp = -1;
  }
  for (i=1;i<vs->n->n;i++) {
    PetscCall(VecGetSize(vs->x[i-1],&s0));
    s += s0;
    PetscCall(VecMax(vs->x[i],idx?&idxp:NULL,&z0));
    if (zp < z0) {
      if (idx) *idx = s+idxp;
      zp = z0;
    }
  }
  if (z) *z = zp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecMin_Comp(Vec v,PetscInt *idx,PetscReal *z)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       idxp,s=0,s0;
  PetscReal      zp,z0;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  if (!idx && !z) PetscFunctionReturn(PETSC_SUCCESS);

  if (vs->n->n > 0) PetscCall(VecMin(vs->x[0],idx?&idxp:NULL,&zp));
  else {
    zp = PETSC_MAX_REAL;
    if (idx) idxp = -1;
  }
  for (i=1;i<vs->n->n;i++) {
    PetscCall(VecGetSize(vs->x[i-1],&s0));
    s += s0;
    PetscCall(VecMin(vs->x[i],idx?&idxp:NULL,&z0));
    if (zp > z0) {
      if (idx) *idx = s+idxp;
      zp = z0;
    }
  }
  if (z) *z = zp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VecMaxPointwiseDivide_Comp(Vec v,Vec w,PetscReal *m)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscReal      work;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,2);
  if (!m || vs->n->n == 0) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecMaxPointwiseDivide(vs->x[0],ws->x[0],m));
  for (i=1;i<vs->n->n;i++) {
    PetscCall(VecMaxPointwiseDivide(vs->x[i],ws->x[i],&work));
    *m = PetscMax(*m,work);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#define __QUOTEME__(x) #x
#define __COMPOSE2__(A,B) A##B
#define __COMPOSE3__(A,B,C) A##B##C

#define __FUNC_TEMPLATE1__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v) \
{ \
  Vec_Comp        *vs = (Vec_Comp*)v->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  for (i=0;i<vs->n->n;i++) { \
    PetscCall(__COMPOSE2__(Vec,NAME)(vs->x[i])); \
  } \
  PetscFunctionReturn(PETSC_SUCCESS);\
}

__FUNC_TEMPLATE1__(Conjugate)
__FUNC_TEMPLATE1__(Reciprocal)
__FUNC_TEMPLATE1__(SqrtAbs)
__FUNC_TEMPLATE1__(Abs)
__FUNC_TEMPLATE1__(Exp)
__FUNC_TEMPLATE1__(Log)

#define __FUNC_TEMPLATE2__(NAME,T0) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,T0 __a) \
{ \
  Vec_Comp        *vs = (Vec_Comp*)v->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  for (i=0;i<vs->n->n;i++) { \
    PetscCall(__COMPOSE2__(Vec,NAME)(vs->x[i],__a)); \
  } \
  PetscFunctionReturn(PETSC_SUCCESS);\
}

__FUNC_TEMPLATE2__(Set,PetscScalar)
__FUNC_TEMPLATE2__(View,PetscViewer)
__FUNC_TEMPLATE2__(Scale,PetscScalar)
__FUNC_TEMPLATE2__(SetRandom,PetscRandom)
__FUNC_TEMPLATE2__(Shift,PetscScalar)

#define __FUNC_TEMPLATE3__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,Vec w) \
{ \
  Vec_Comp        *vs = (Vec_Comp*)v->data,\
                  *ws = (Vec_Comp*)w->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  SlepcValidVecComp(w,2); \
  for (i=0;i<vs->n->n;i++) { \
    PetscCall(__COMPOSE2__(Vec,NAME)(vs->x[i],ws->x[i])); \
  } \
  PetscFunctionReturn(PETSC_SUCCESS);\
}

__FUNC_TEMPLATE3__(Copy)
__FUNC_TEMPLATE3__(Swap)

#define __FUNC_TEMPLATE4__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,Vec w,Vec z) \
{ \
  Vec_Comp        *vs = (Vec_Comp*)v->data, \
                  *ws = (Vec_Comp*)w->data, \
                  *zs = (Vec_Comp*)z->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  SlepcValidVecComp(w,2); \
  SlepcValidVecComp(z,3); \
  for (i=0;i<vs->n->n;i++) { \
    PetscCall(__COMPOSE2__(Vec,NAME)(vs->x[i],ws->x[i],zs->x[i])); \
  } \
  PetscFunctionReturn(PETSC_SUCCESS);\
}

__FUNC_TEMPLATE4__(PointwiseMax)
__FUNC_TEMPLATE4__(PointwiseMaxAbs)
__FUNC_TEMPLATE4__(PointwiseMin)
__FUNC_TEMPLATE4__(PointwiseMult)
__FUNC_TEMPLATE4__(PointwiseDivide)
