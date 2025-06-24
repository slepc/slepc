/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcvec.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = Sys */

#if !defined(PETSC_USE_DEBUG)

#define SlepcValidVecComp(y,arg) do {(void)(y);} while (0)

#else

#define SlepcValidVecComp(y,arg) \
  do { \
    PetscCheck(((Vec_Comp*)(y)->data)->nx>=((Vec_Comp*)(y)->data)->n->n,PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Invalid number of subvectors required: Parameter #%d",arg); \
  } while (0)

#endif

/* Contexts for VecComp */
typedef struct {
  PetscInt      n;        /* number of active subvectors */
  PetscInt      N;        /* virtual global size */
  PetscInt      lN;       /* virtual local size */
  PetscInt      friends;  /* number of vectors sharing this structure */
} Vec_Comp_N;

typedef struct {
  Vec           *x;       /* the vectors */
  PetscInt      nx;       /* number of available subvectors */
  Vec_Comp_N    *n;       /* structure shared by friend vectors */
} Vec_Comp;

/* Operations implemented in VecComp */
SLEPC_INTERN PetscErrorCode VecDuplicateVecs_Comp(Vec,PetscInt,Vec*[]);
SLEPC_INTERN PetscErrorCode VecDestroyVecs_Comp(PetscInt,Vec[]);
SLEPC_INTERN PetscErrorCode VecDuplicate_Comp(Vec,Vec*);
SLEPC_INTERN PetscErrorCode VecDestroy_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecSet_Comp(Vec,PetscScalar);
SLEPC_INTERN PetscErrorCode VecView_Comp(Vec,PetscViewer);
SLEPC_INTERN PetscErrorCode VecScale_Comp(Vec,PetscScalar);
SLEPC_INTERN PetscErrorCode VecCopy_Comp(Vec,Vec);
SLEPC_INTERN PetscErrorCode VecSwap_Comp(Vec,Vec);
SLEPC_INTERN PetscErrorCode VecAXPY_Comp(Vec,PetscScalar,Vec);
SLEPC_INTERN PetscErrorCode VecAYPX_Comp(Vec,PetscScalar,Vec);
SLEPC_INTERN PetscErrorCode VecAXPBY_Comp(Vec,PetscScalar,PetscScalar,Vec);
SLEPC_INTERN PetscErrorCode VecMAXPY_Comp(Vec,PetscInt,const PetscScalar*,Vec*);
SLEPC_INTERN PetscErrorCode VecWAXPY_Comp(Vec,PetscScalar,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecAXPBYPCZ_Comp(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecPointwiseMult_Comp(Vec,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecPointwiseDivide_Comp(Vec,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecGetSize_Comp(Vec,PetscInt*);
SLEPC_INTERN PetscErrorCode VecGetLocalSize_Comp(Vec,PetscInt*);
SLEPC_INTERN PetscErrorCode VecMax_Comp(Vec,PetscInt*,PetscReal*);
SLEPC_INTERN PetscErrorCode VecMin_Comp(Vec,PetscInt*,PetscReal*);
SLEPC_INTERN PetscErrorCode VecSetRandom_Comp(Vec,PetscRandom);
SLEPC_INTERN PetscErrorCode VecConjugate_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecReciprocal_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecMaxPointwiseDivide_Comp(Vec,Vec,PetscReal*);
SLEPC_INTERN PetscErrorCode VecPointwiseMax_Comp(Vec,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecPointwiseMaxAbs_Comp(Vec,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecPointwiseMin_Comp(Vec,Vec,Vec);
SLEPC_INTERN PetscErrorCode VecDotNorm2_Comp_Seq(Vec,Vec,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode VecDotNorm2_Comp_MPI(Vec,Vec,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode VecSqrtAbs_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecAbs_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecExp_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecLog_Comp(Vec);
SLEPC_INTERN PetscErrorCode VecShift_Comp(Vec,PetscScalar);
SLEPC_EXTERN PetscErrorCode VecCreate_Comp(Vec);

/* VecPool */
typedef struct VecPool_ {
  Vec      v;              /* template vector */
  Vec      *vecs;          /* pool of vectors */
  PetscInt n;              /* size of vecs */
  PetscInt used;           /* number of already used vectors */
  PetscInt guess;          /* expected maximum number of vectors */
  struct VecPool_ *next;   /* list of pool of vectors */
} VecPool_;
typedef VecPool_* VecPool;

SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcVecPoolCreate(Vec,PetscInt,VecPool*);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcVecPoolDestroy(VecPool*);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcVecPoolGetVecs(VecPool,PetscInt,Vec**);
SLEPC_SINGLE_LIBRARY_INTERN PetscErrorCode SlepcVecPoolRestoreVecs(VecPool,PetscInt,Vec**);
