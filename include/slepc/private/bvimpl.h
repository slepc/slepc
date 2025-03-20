/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcbv.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = BV */

SLEPC_EXTERN PetscBool BVRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode BVRegisterAll(void);

SLEPC_EXTERN PetscLogEvent BV_Create,BV_Copy,BV_Mult,BV_MultVec,BV_MultInPlace,BV_Dot,BV_DotVec,BV_Orthogonalize,BV_OrthogonalizeVec,BV_Scale,BV_Norm,BV_NormVec,BV_Normalize,BV_SetRandom,BV_MatMult,BV_MatMultVec,BV_MatProject,BV_SVDAndRank;

typedef struct _BVOps *BVOps;

struct _BVOps {
  PetscErrorCode (*mult)(BV,PetscScalar,PetscScalar,BV,Mat);
  PetscErrorCode (*multvec)(BV,PetscScalar,PetscScalar,Vec,PetscScalar*);
  PetscErrorCode (*multinplace)(BV,Mat,PetscInt,PetscInt);
  PetscErrorCode (*multinplacetrans)(BV,Mat,PetscInt,PetscInt);
  PetscErrorCode (*dot)(BV,BV,Mat);
  PetscErrorCode (*dotvec)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_local)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_begin)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_end)(BV,Vec,PetscScalar*);
  PetscErrorCode (*scale)(BV,PetscInt,PetscScalar);
  PetscErrorCode (*norm)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_local)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_begin)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_end)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*normalize)(BV,PetscScalar*);
  PetscErrorCode (*matmult)(BV,Mat,BV);
  PetscErrorCode (*copy)(BV,BV);
  PetscErrorCode (*copycolumn)(BV,PetscInt,PetscInt);
  PetscErrorCode (*resize)(BV,PetscInt,PetscBool);
  PetscErrorCode (*getcolumn)(BV,PetscInt,Vec*);
  PetscErrorCode (*restorecolumn)(BV,PetscInt,Vec*);
  PetscErrorCode (*getarray)(BV,PetscScalar**);
  PetscErrorCode (*restorearray)(BV,PetscScalar**);
  PetscErrorCode (*getarrayread)(BV,const PetscScalar**);
  PetscErrorCode (*restorearrayread)(BV,const PetscScalar**);
  PetscErrorCode (*restoresplit)(BV,BV*,BV*);
  PetscErrorCode (*restoresplitrows)(BV,IS,IS,BV*,BV*);
  PetscErrorCode (*gramschmidt)(BV,PetscInt,Vec,PetscBool*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*);
  PetscErrorCode (*getmat)(BV,Mat*);
  PetscErrorCode (*restoremat)(BV,Mat*);
  PetscErrorCode (*duplicate)(BV,BV);
  PetscErrorCode (*create)(BV);
  PetscErrorCode (*setfromoptions)(BV,PetscOptionItems);
  PetscErrorCode (*view)(BV,PetscViewer);
  PetscErrorCode (*destroy)(BV);
};

struct _p_BV {
  PETSCHEADER(struct _BVOps);
  /*------------------------- User parameters --------------------------*/
  PetscLayout        map;          /* layout of columns */
  VecType            vtype;        /* vector type */
  PetscInt           n,N;          /* dimensions of vectors (local, global) */
  PetscInt           m;            /* number of vectors */
  PetscInt           l;            /* number of leading columns */
  PetscInt           k;            /* number of active columns */
  PetscInt           nc;           /* number of constraints */
  PetscInt           ld;           /* leading dimension */
  BVOrthogType       orthog_type;  /* the method of vector orthogonalization */
  BVOrthogRefineType orthog_ref;   /* refinement method */
  PetscReal          orthog_eta;   /* refinement threshold */
  BVOrthogBlockType  orthog_block; /* the method of block orthogonalization */
  Mat                matrix;       /* inner product matrix */
  PetscBool          indef;        /* matrix is indefinite */
  BVMatMultType      vmm;          /* version of matmult operation */
  PetscBool          rrandom;      /* reproducible random vectors */
  PetscReal          deftol;       /* tolerance for BV_SafeSqrt */

  /*---------------------- Cached data and workspace -------------------*/
  Vec                buffer;       /* buffer vector used in orthogonalization */
  Mat                Abuffer;      /* auxiliary seqdense matrix that wraps the buffer */
  Vec                Bx;           /* result of matrix times a vector x */
  PetscObjectId      xid;          /* object id of vector x */
  PetscObjectState   xstate;       /* state of vector x */
  Vec                cv[2];        /* column vectors obtained with BVGetColumn() */
  PetscInt           ci[2];        /* column indices of obtained vectors */
  PetscObjectState   st[2];        /* state of obtained vectors */
  PetscObjectId      id[2];        /* object id of obtained vectors */
  PetscScalar        *h,*c;        /* orthogonalization coefficients */
  Vec                omega;        /* signature matrix values for indefinite case */
  PetscBool          defersfo;     /* deferred call to setfromoptions */
  BV                 cached;       /* cached BV to store result of matrix times BV */
  PetscObjectState   bvstate;      /* state of BV when BVApplyMatrixBV() was called */
  BV                 L,R;          /* BV objects obtained with BVGetSplit/Rows() */
  PetscObjectState   lstate,rstate;/* state of L and R when BVGetSplit/Rows() was called */
  PetscInt           lsplit;       /* value of l when BVGetSplit() was called (-1 if BVGetSplitRows()) */
  PetscInt           issplit;      /* !=0 if BV is from split (1=left, 2=right, -1=top, -2=bottom) */
  BV                 splitparent;  /* my parent if I am a split BV */
  PetscRandom        rand;         /* random number generator */
  Mat                Acreate;      /* matrix given at BVCreateFromMat() */
  Mat                Aget;         /* matrix returned for BVGetMat() */
  PetscBool          cuda;         /* true if NVIDIA GPU must be used */
  PetscBool          hip;          /* true if AMD GPU must be used */
  PetscBool          sfocalled;    /* setfromoptions has been called */
  PetscScalar        *work;
  PetscInt           lwork;
  void               *data;
};

/*
  BV_SafeSqrt - Computes the square root of a scalar value alpha, which is
  assumed to be z'*B*z. The result is
    if definite inner product:     res = sqrt(alpha)
    if indefinite inner product:   res = sgn(alpha)*sqrt(abs(alpha))
*/
static inline PetscErrorCode BV_SafeSqrt(BV bv,PetscScalar alpha,PetscReal *res)
{
  PetscReal      absal,realp;
  const char     *msg;

  PetscFunctionBegin;
  absal = PetscAbsScalar(alpha);
  realp = PetscRealPart(alpha);
  if (PetscUnlikely(absal<PETSC_MACHINE_EPSILON)) PetscCall(PetscInfo(bv,"Zero norm %g, either the vector is zero or a semi-inner product is being used\n",(double)absal));
#if defined(PETSC_USE_COMPLEX)
  PetscCheck(PetscAbsReal(PetscImaginaryPart(alpha))<bv->deftol || PetscAbsReal(PetscImaginaryPart(alpha))/absal<10*bv->deftol,PetscObjectComm((PetscObject)bv),PETSC_ERR_USER_INPUT,"The inner product is not well defined: nonzero imaginary part %g",(double)PetscImaginaryPart(alpha));
#endif
  if (PetscUnlikely(bv->indef)) {
    *res = (realp<0.0)? -PetscSqrtReal(-realp): PetscSqrtReal(realp);
  } else {
    msg = bv->matrix? "The inner product is not well defined: indefinite matrix %g": "Invalid inner product: %g";
    PetscCheck(realp>-bv->deftol,PetscObjectComm((PetscObject)bv),PETSC_ERR_USER_INPUT,msg,(double)realp);
    *res = (realp<0.0)? 0.0: PetscSqrtReal(realp);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_IPMatMult - Multiply a vector x by the inner-product matrix, cache the
  result in Bx.
*/
static inline PetscErrorCode BV_IPMatMult(BV bv,Vec x)
{
  PetscFunctionBegin;
  if (((PetscObject)x)->id != bv->xid || ((PetscObject)x)->state != bv->xstate) {
    if (PetscUnlikely(!bv->Bx)) PetscCall(MatCreateVecs(bv->matrix,&bv->Bx,NULL));
    PetscCall(MatMult(bv->matrix,x,bv->Bx));
    PetscCall(PetscObjectGetId((PetscObject)x,&bv->xid));
    PetscCall(VecGetState(x,&bv->xstate));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_IPMatMultBV - Multiply BV by the inner-product matrix, cache the
  result internally in bv->cached.
*/
static inline PetscErrorCode BV_IPMatMultBV(BV bv)
{
  PetscFunctionBegin;
  PetscCall(BVGetCachedBV(bv,&bv->cached));
  if (((PetscObject)bv)->state != bv->bvstate || bv->l != bv->cached->l || bv->k != bv->cached->k) {
    PetscCall(BVSetActiveColumns(bv->cached,bv->l,bv->k));
    if (bv->matrix) PetscCall(BVMatMult(bv,bv->matrix,bv->cached));
    else PetscCall(BVCopy(bv,bv->cached));
    bv->bvstate = ((PetscObject)bv)->state;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_AllocateCoeffs - Allocate orthogonalization coefficients if not done already.
*/
static inline PetscErrorCode BV_AllocateCoeffs(BV bv)
{
  PetscFunctionBegin;
  if (!bv->h) PetscCall(PetscMalloc2(bv->nc+bv->m,&bv->h,bv->nc+bv->m,&bv->c));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_AllocateSignature - Allocate signature coefficients if not done already.
*/
static inline PetscErrorCode BV_AllocateSignature(BV bv)
{
  PetscFunctionBegin;
  if (bv->indef && !bv->omega) {
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      PetscCall(VecCreateSeqCUDA(PETSC_COMM_SELF,bv->nc+bv->m,&bv->omega));
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_PLIB,"Something wrong happened");
#endif
    } else if (bv->hip) {
#if defined(PETSC_HAVE_HIP)
      PetscCall(VecCreateSeqHIP(PETSC_COMM_SELF,bv->nc+bv->m,&bv->omega));
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_PLIB,"Something wrong happened");
#endif
    } else PetscCall(VecCreateSeq(PETSC_COMM_SELF,bv->nc+bv->m,&bv->omega));
    PetscCall(VecSet(bv->omega,1.0));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_SetMatrixDiagonal - sets the inner product matrix for BV as a diagonal matrix
  with the diagonal specified by vector vomega, using the same matrix type as matrix M
*/
static inline PetscErrorCode BV_SetMatrixDiagonal(BV bv,Vec vomega,Mat M)
{
  Mat      Omega;
  MatType  Mtype;

  PetscFunctionBegin;
  PetscCall(MatGetType(M,&Mtype));
  PetscCall(MatCreate(PetscObjectComm((PetscObject)bv),&Omega));
  PetscCall(MatSetSizes(Omega,bv->n,bv->n,bv->N,bv->N));
  PetscCall(MatSetType(Omega,Mtype));
  PetscCall(MatDiagonalSet(Omega,vomega,INSERT_VALUES));
  PetscCall(BVSetMatrix(bv,Omega,PETSC_TRUE));
  PetscCall(MatDestroy(&Omega));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BVAvailableVec: First (0) or second (1) vector available for
  getcolumn operation (or -1 if both vectors already fetched).
*/
#define BVAvailableVec (((bv->ci[0]==-bv->nc-1)? 0: (bv->ci[1]==-bv->nc-1)? 1: -1))

/*
    Macros to test valid BV arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define BVCheckSizes(h,arg) do {(void)(h);} while (0)
#define BVCheckOp(h,arg,op) do {(void)(h);} while (0)

#else

#define BVCheckSizes(h,arg) \
  do { \
    PetscCheck((h)->m,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"BV sizes have not been defined: Parameter #%d",arg); \
  } while (0)

#define BVCheckOp(h,arg,op) \
  do { \
    PetscCheck((h)->ops->op,PetscObjectComm((PetscObject)(h)),PETSC_ERR_SUP,"Operation not implemented in this BV type: Parameter #%d",arg); \
  } while (0)

#endif

SLEPC_INTERN PetscErrorCode BVView_Vecs(BV,PetscViewer);

SLEPC_INTERN PetscErrorCode BVAllocateWork_Private(BV,PetscInt);

SLEPC_INTERN PetscErrorCode BVMult_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVMultVec_BLAS_Private(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVMultInPlace_Vecs_Private(BV,PetscInt,PetscInt,PetscInt,Vec*,const PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVAXPY_BLAS_Private(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVDot_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVDotVec_BLAS_Private(BV,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVScale_BLAS_Private(BV,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_LAPACK_Private(BV,PetscInt,PetscInt,const PetscScalar*,PetscInt,NormType,PetscReal*,PetscBool);
SLEPC_INTERN PetscErrorCode BVNormalize_LAPACK_Private(BV,PetscInt,PetscInt,const PetscScalar*,PetscInt,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVGetMat_Default(BV,Mat*);
SLEPC_INTERN PetscErrorCode BVRestoreMat_Default(BV,Mat*);
SLEPC_INTERN PetscErrorCode BVMatCholInv_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVMatTriInv_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVMatSVQB_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVOrthogonalize_LAPACK_TSQR(BV,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVOrthogonalize_LAPACK_TSQR_OnlyR(BV,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);

/* reduction operations used in BVOrthogonalize and BVNormalize */
SLEPC_EXTERN MPI_Op MPIU_TSQR, MPIU_LAPY2;
SLEPC_EXTERN void MPIAPI SlepcGivensPacked(void*,void*,PetscMPIInt*,MPI_Datatype*);
SLEPC_EXTERN void MPIAPI SlepcPythag(void*,void*,PetscMPIInt*,MPI_Datatype*);

/*
   BV_CleanCoefficients_Default - Sets to zero all entries of column j of the bv buffer
*/
static inline PetscErrorCode BV_CleanCoefficients_Default(BV bv,PetscInt j,PetscScalar *h)
{
  PetscScalar    *hh=h,*a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecGetArray(bv->buffer,&a));
    hh = a + j*(bv->nc+bv->m);
  }
  for (i=0;i<bv->nc+j;i++) hh[i] = 0.0;
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_AddCoefficients_Default - Add the contents of the scratch (0-th column) of the bv buffer
   into column j of the bv buffer
*/
static inline PetscErrorCode BV_AddCoefficients_Default(BV bv,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscScalar    *hh=h,*cc=c;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecGetArray(bv->buffer,&cc));
    hh = cc + j*(bv->nc+bv->m);
  }
  for (i=0;i<bv->nc+j;i++) hh[i] += cc[i];
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&cc));
  PetscCall(PetscLogFlops(1.0*(bv->nc+j)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SetValue_Default - Sets value in row j (counted after the constraints) of column k
   of the coefficients array
*/
static inline PetscErrorCode BV_SetValue_Default(BV bv,PetscInt j,PetscInt k,PetscScalar *h,PetscScalar value)
{
  PetscScalar    *hh=h,*a;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecGetArray(bv->buffer,&a));
    hh = a + k*(bv->nc+bv->m);
  }
  hh[bv->nc+j] = value;
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SquareSum_Default - Returns the value h'*h, where h represents the contents of the
   coefficients array (up to position j)
*/
static inline PetscErrorCode BV_SquareSum_Default(BV bv,PetscInt j,PetscScalar *h,PetscReal *sum)
{
  PetscScalar    *hh=h;
  PetscInt       i;

  PetscFunctionBegin;
  *sum = 0.0;
  if (!h) PetscCall(VecGetArray(bv->buffer,&hh));
  for (i=0;i<bv->nc+j;i++) *sum += PetscRealPart(hh[i]*PetscConj(hh[i]));
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&hh));
  PetscCall(PetscLogFlops(2.0*(bv->nc+j)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_ApplySignature_Default - Computes the pointwise product h*omega, where h represents
   the contents of the coefficients array (up to position j) and omega is the signature;
   if inverse=TRUE then the operation is h/omega
*/
static inline PetscErrorCode BV_ApplySignature_Default(BV bv,PetscInt j,PetscScalar *h,PetscBool inverse)
{
  PetscScalar       *hh=h;
  PetscInt          i;
  const PetscScalar *omega;

  PetscFunctionBegin;
  if (PetscUnlikely(!(bv->nc+j))) PetscFunctionReturn(PETSC_SUCCESS);
  if (!h) PetscCall(VecGetArray(bv->buffer,&hh));
  PetscCall(VecGetArrayRead(bv->omega,&omega));
  if (inverse) for (i=0;i<bv->nc+j;i++) hh[i] /= PetscRealPart(omega[i]);
  else for (i=0;i<bv->nc+j;i++) hh[i] *= PetscRealPart(omega[i]);
  PetscCall(VecRestoreArrayRead(bv->omega,&omega));
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&hh));
  PetscCall(PetscLogFlops(1.0*(bv->nc+j)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SquareRoot_Default - Returns the square root of position j (counted after the constraints)
   of the coefficients array
*/
static inline PetscErrorCode BV_SquareRoot_Default(BV bv,PetscInt j,PetscScalar *h,PetscReal *beta)
{
  PetscScalar    *hh=h;

  PetscFunctionBegin;
  if (!h) PetscCall(VecGetArray(bv->buffer,&hh));
  PetscCall(BV_SafeSqrt(bv,hh[bv->nc+j],beta));
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&hh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_StoreCoefficients_Default - Copy the contents of the coefficients array to an array dest
   provided by the caller (only values from l to j are copied)
*/
static inline PetscErrorCode BV_StoreCoefficients_Default(BV bv,PetscInt j,PetscScalar *h,PetscScalar *dest)
{
  PetscScalar    *hh=h,*a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    PetscCall(VecGetArray(bv->buffer,&a));
    hh = a + j*(bv->nc+bv->m);
  }
  for (i=bv->l;i<j;i++) dest[i-bv->l] = hh[bv->nc+i];
  if (!h) PetscCall(VecRestoreArray(bv->buffer,&a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_GetEigenvector - retrieves k-th eigenvector from basis vectors V.
  The argument eigi is the imaginary part of the corresponding eigenvalue.
*/
static inline PetscErrorCode BV_GetEigenvector(BV V,PetscInt k,PetscScalar eigi,Vec Vr,Vec Vi)
{
  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  (void)eigi;
  if (Vr) PetscCall(BVCopyVec(V,k,Vr));
  if (Vi) PetscCall(VecSet(Vi,0.0));
#else
  if (eigi > 0.0) { /* first value of conjugate pair */
    if (Vr) PetscCall(BVCopyVec(V,k,Vr));
    if (Vi) PetscCall(BVCopyVec(V,k+1,Vi));
  } else if (eigi < 0.0) { /* second value of conjugate pair */
    if (Vr) PetscCall(BVCopyVec(V,k-1,Vr));
    if (Vi) {
      PetscCall(BVCopyVec(V,k,Vi));
      PetscCall(VecScale(Vi,-1.0));
    }
  } else { /* real eigenvalue */
    if (Vr) PetscCall(BVCopyVec(V,k,Vr));
    if (Vi) PetscCall(VecSet(Vi,0.0));
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_OrthogonalizeColumn_Safe - this is intended for cases where we know that
   the resulting vector is going to be numerically zero, so normalization or
   iterative refinement may cause problems in parallel (collective operation
   not being called by all processes)
*/
static inline PetscErrorCode BV_OrthogonalizeColumn_Safe(BV bv,PetscInt j,PetscScalar *H,PetscReal *norm,PetscBool *lindep)
{
  BVOrthogRefineType orthog_ref;

  PetscFunctionBegin;
  PetscCall(PetscInfo(bv,"Orthogonalizing column %" PetscInt_FMT " without refinement\n",j));
  orthog_ref     = bv->orthog_ref;
  bv->orthog_ref = BV_ORTHOG_REFINE_NEVER;  /* avoid refinement */
  PetscCall(BVOrthogonalizeColumn(bv,j,H,NULL,NULL));
  bv->orthog_ref = orthog_ref;  /* restore refinement setting */
  if (norm)   *norm  = 0.0;
  if (lindep) *lindep = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_SetDefaultLD - set the default value of the leading dimension, based on
   the local size.
*/
static inline PetscErrorCode BV_SetDefaultLD(BV bv,PetscInt nloc)
{
  size_t bytes,align;

  PetscFunctionBegin;
  if (bv->ld) {   /* set by user */
    PetscCheck(bv->ld>=nloc,PetscObjectComm((PetscObject)bv),PETSC_ERR_USER_INPUT,"The leading dimension %" PetscInt_FMT " should be larger or equal to the local number of rows %" PetscInt_FMT,bv->ld,nloc);
  } else {
    align = PetscMax(PETSC_MEMALIGN,16);   /* assume that CUDA requires 16-byte alignment */
    bytes = (nloc*sizeof(PetscScalar) + align - 1) & ~(align - 1);
    PetscCall(PetscIntCast(bytes/sizeof(PetscScalar),&bv->ld));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#if defined(PETSC_HAVE_CUDA)
/*
   BV_MatDenseCUDAGetArrayRead - if Q is MATSEQDENSE it will allocate memory on the
   GPU and copy the contents; otherwise, calls MatDenseCUDAGetArrayRead()
*/
static inline PetscErrorCode BV_MatDenseCUDAGetArrayRead(BV bv,Mat Q,const PetscScalar **d_q)
{
  const PetscScalar *q;
  PetscInt          ldq,mq;
  PetscCuBLASInt    ldq_=0;
  PetscBool         matiscuda;

  PetscFunctionBegin;
  (void)bv; // avoid unused parameter warning
  PetscCall(MatGetSize(Q,NULL,&mq));
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(PetscCuBLASIntCast(ldq,&ldq_));
  PetscCall(PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSECUDA,&matiscuda));
  if (matiscuda) PetscCall(MatDenseCUDAGetArrayRead(Q,d_q));
  else {
    PetscCall(MatDenseGetArrayRead(Q,&q));
    PetscCallCUDA(cudaMalloc((void**)d_q,ldq*mq*sizeof(PetscScalar)));
    PetscCallCUDA(cudaMemcpy((void*)*d_q,q,ldq*mq*sizeof(PetscScalar),cudaMemcpyHostToDevice));
    PetscCall(PetscLogCpuToGpu(ldq*mq*sizeof(PetscScalar)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_MatDenseCUDARestoreArrayRead - restores the pointer obtained with BV_MatDenseCUDAGetArrayRead(),
   freeing the GPU memory in case of MATSEQDENSE
*/
static inline PetscErrorCode BV_MatDenseCUDARestoreArrayRead(BV bv,Mat Q,const PetscScalar **d_q)
{
  PetscBool matiscuda;

  PetscFunctionBegin;
  (void)bv; // avoid unused parameter warning
  PetscCall(PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSECUDA,&matiscuda));
  if (matiscuda) PetscCall(MatDenseCUDARestoreArrayRead(Q,d_q));
  else {
    PetscCall(MatDenseRestoreArrayRead(Q,NULL));
    PetscCallCUDA(cudaFree((void*)*d_q));
    *d_q = NULL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode BVMult_BLAS_CUDA(BV,PetscInt,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVMultVec_BLAS_CUDA(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_BLAS_CUDA(BV,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVAXPY_BLAS_CUDA(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVDot_BLAS_CUDA(BV,PetscInt,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVDotVec_BLAS_CUDA(BV,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVScale_BLAS_CUDA(BV,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_BLAS_CUDA(BV,PetscInt,const PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNormalize_BLAS_CUDA(BV,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*);

SLEPC_INTERN PetscErrorCode BV_CleanCoefficients_CUDA(BV,PetscInt,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_AddCoefficients_CUDA(BV,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_SetValue_CUDA(BV,PetscInt,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BV_SquareSum_CUDA(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_ApplySignature_CUDA(BV,PetscInt,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BV_SquareRoot_CUDA(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_StoreCoefficients_CUDA(BV,PetscInt,PetscScalar*,PetscScalar*);
#define BV_CleanCoefficients(a,b,c)   ((a)->cuda?BV_CleanCoefficients_CUDA:BV_CleanCoefficients_Default)((a),(b),(c))
#define BV_AddCoefficients(a,b,c,d)   ((a)->cuda?BV_AddCoefficients_CUDA:BV_AddCoefficients_Default)((a),(b),(c),(d))
#define BV_SetValue(a,b,c,d,e)        ((a)->cuda?BV_SetValue_CUDA:BV_SetValue_Default)((a),(b),(c),(d),(e))
#define BV_SquareSum(a,b,c,d)         ((a)->cuda?BV_SquareSum_CUDA:BV_SquareSum_Default)((a),(b),(c),(d))
#define BV_ApplySignature(a,b,c,d)    ((a)->cuda?BV_ApplySignature_CUDA:BV_ApplySignature_Default)((a),(b),(c),(d))
#define BV_SquareRoot(a,b,c,d)        ((a)->cuda?BV_SquareRoot_CUDA:BV_SquareRoot_Default)((a),(b),(c),(d))
#define BV_StoreCoefficients(a,b,c,d) ((a)->cuda?BV_StoreCoefficients_CUDA:BV_StoreCoefficients_Default)((a),(b),(c),(d))

#elif defined(PETSC_HAVE_HIP)
#include <petscdevice_cupm.h>
/*
   BV_MatDenseHIPGetArrayRead - if Q is MATSEQDENSE it will allocate memory on the
   GPU and copy the contents; otherwise, calls MatDenseHIPGetArrayRead()
*/
static inline PetscErrorCode BV_MatDenseHIPGetArrayRead(BV bv,Mat Q,const PetscScalar **d_q)
{
  const PetscScalar *q;
  PetscInt          ldq,mq;
  PetscCuBLASInt    ldq_=0;
  PetscBool         matiship;

  PetscFunctionBegin;
  (void)bv; // avoid unused parameter warning
  PetscCall(MatGetSize(Q,NULL,&mq));
  PetscCall(MatDenseGetLDA(Q,&ldq));
  PetscCall(PetscHipBLASIntCast(ldq,&ldq_));
  PetscCall(PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSEHIP,&matiship));
  if (matiship) PetscCall(MatDenseHIPGetArrayRead(Q,d_q));
  else {
    PetscCall(MatDenseGetArrayRead(Q,&q));
    PetscCallHIP(hipMalloc((void**)d_q,ldq*mq*sizeof(PetscScalar)));
    PetscCallHIP(hipMemcpy((void*)*d_q,q,ldq*mq*sizeof(PetscScalar),hipMemcpyHostToDevice));
    PetscCall(PetscLogCpuToGpu(ldq*mq*sizeof(PetscScalar)));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_MatDenseHIPRestoreArrayRead - restores the pointer obtained with BV_MatDenseHIPGetArrayRead(),
   freeing the GPU memory in case of MATSEQDENSE
*/
static inline PetscErrorCode BV_MatDenseHIPRestoreArrayRead(BV bv,Mat Q,const PetscScalar **d_q)
{
  PetscBool matiship;

  PetscFunctionBegin;
  (void)bv; // avoid unused parameter warning
  PetscCall(PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSEHIP,&matiship));
  if (matiship) PetscCall(MatDenseHIPRestoreArrayRead(Q,d_q));
  else {
    PetscCall(MatDenseRestoreArrayRead(Q,NULL));
    PetscCallHIP(hipFree((void*)*d_q));
    *d_q = NULL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode BVMult_BLAS_HIP(BV,PetscInt,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVMultVec_BLAS_HIP(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_BLAS_HIP(BV,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVAXPY_BLAS_HIP(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscInt,PetscScalar,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVDot_BLAS_HIP(BV,PetscInt,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVDotVec_BLAS_HIP(BV,PetscInt,PetscInt,const PetscScalar*,PetscInt,const PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVScale_BLAS_HIP(BV,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_BLAS_HIP(BV,PetscInt,const PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNormalize_BLAS_HIP(BV,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*);

SLEPC_INTERN PetscErrorCode BV_CleanCoefficients_HIP(BV,PetscInt,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_AddCoefficients_HIP(BV,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_SetValue_HIP(BV,PetscInt,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BV_SquareSum_HIP(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_ApplySignature_HIP(BV,PetscInt,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BV_SquareRoot_HIP(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_StoreCoefficients_HIP(BV,PetscInt,PetscScalar*,PetscScalar*);
#define BV_CleanCoefficients(a,b,c)   ((a)->hip?BV_CleanCoefficients_HIP:BV_CleanCoefficients_Default)((a),(b),(c))
#define BV_AddCoefficients(a,b,c,d)   ((a)->hip?BV_AddCoefficients_HIP:BV_AddCoefficients_Default)((a),(b),(c),(d))
#define BV_SetValue(a,b,c,d,e)        ((a)->hip?BV_SetValue_HIP:BV_SetValue_Default)((a),(b),(c),(d),(e))
#define BV_SquareSum(a,b,c,d)         ((a)->hip?BV_SquareSum_HIP:BV_SquareSum_Default)((a),(b),(c),(d))
#define BV_ApplySignature(a,b,c,d)    ((a)->hip?BV_ApplySignature_HIP:BV_ApplySignature_Default)((a),(b),(c),(d))
#define BV_SquareRoot(a,b,c,d)        ((a)->hip?BV_SquareRoot_HIP:BV_SquareRoot_Default)((a),(b),(c),(d))
#define BV_StoreCoefficients(a,b,c,d) ((a)->hip?BV_StoreCoefficients_HIP:BV_StoreCoefficients_Default)((a),(b),(c),(d))

#else /* CPU */
#define BV_CleanCoefficients(a,b,c)   BV_CleanCoefficients_Default((a),(b),(c))
#define BV_AddCoefficients(a,b,c,d)   BV_AddCoefficients_Default((a),(b),(c),(d))
#define BV_SetValue(a,b,c,d,e)        BV_SetValue_Default((a),(b),(c),(d),(e))
#define BV_SquareSum(a,b,c,d)         BV_SquareSum_Default((a),(b),(c),(d))
#define BV_ApplySignature(a,b,c,d)    BV_ApplySignature_Default((a),(b),(c),(d))
#define BV_SquareRoot(a,b,c,d)        BV_SquareRoot_Default((a),(b),(c),(d))
#define BV_StoreCoefficients(a,b,c,d) BV_StoreCoefficients_Default((a),(b),(c),(d))
#endif /* PETSC_HAVE_CUDA */
