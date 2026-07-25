/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

typedef struct {
  Vec       v;
  PetscBool mpi;    /* true if either VECMPI, VECMPICUDA, or VECMPIHIP */
} BV_SVEC;

#if PetscDefined(HAVE_CUDA)
SLEPC_INTERN PetscErrorCode BVMult_Svec_CUDA(BV,PetscScalar,PetscScalar,BV,Mat);
SLEPC_INTERN PetscErrorCode BVMultVec_Svec_CUDA(BV,PetscScalar,PetscScalar,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_Svec_CUDA(BV,Mat,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVMultInPlaceHermitianTranspose_Svec_CUDA(BV,Mat,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVDot_Svec_CUDA(BV,BV,Mat);
SLEPC_INTERN PetscErrorCode BVDotVec_Svec_CUDA(BV,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVDotVec_Local_Svec_CUDA(BV,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVScale_Svec_CUDA(BV,PetscInt,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_Svec_CUDA(BV,PetscInt,NormType,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNorm_Local_Svec_CUDA(BV,PetscInt,NormType,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNormalize_Svec_CUDA(BV,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMatMult_Svec_CUDA(BV,Mat,BV);
SLEPC_INTERN PetscErrorCode BVCopy_Svec_CUDA(BV,BV);
SLEPC_INTERN PetscErrorCode BVCopyColumn_Svec_CUDA(BV,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVResize_Svec_CUDA(BV,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVGetColumn_Svec_CUDA(BV,PetscInt,Vec*);
SLEPC_INTERN PetscErrorCode BVRestoreColumn_Svec_CUDA(BV,PetscInt,Vec*);
SLEPC_INTERN PetscErrorCode BVRestoreSplit_Svec_CUDA(BV,BV*,BV*);
SLEPC_INTERN PetscErrorCode BVRestoreSplitRows_Svec_CUDA(BV,IS,IS,BV*,BV*);
SLEPC_INTERN PetscErrorCode BVGetMat_Svec_CUDA(BV,Mat*);
SLEPC_INTERN PetscErrorCode BVRestoreMat_Svec_CUDA(BV,Mat*);
#endif

#if PetscDefined(HAVE_HIP)
SLEPC_INTERN PetscErrorCode BVMult_Svec_HIP(BV,PetscScalar,PetscScalar,BV,Mat);
SLEPC_INTERN PetscErrorCode BVMultVec_Svec_HIP(BV,PetscScalar,PetscScalar,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_Svec_HIP(BV,Mat,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVMultInPlaceHermitianTranspose_Svec_HIP(BV,Mat,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVDot_Svec_HIP(BV,BV,Mat);
SLEPC_INTERN PetscErrorCode BVDotVec_Svec_HIP(BV,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVDotVec_Local_Svec_HIP(BV,Vec,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVScale_Svec_HIP(BV,PetscInt,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_Svec_HIP(BV,PetscInt,NormType,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNorm_Local_Svec_HIP(BV,PetscInt,NormType,PetscReal*);
SLEPC_INTERN PetscErrorCode BVNormalize_Svec_HIP(BV,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMatMult_Svec_HIP(BV,Mat,BV);
SLEPC_INTERN PetscErrorCode BVCopy_Svec_HIP(BV,BV);
SLEPC_INTERN PetscErrorCode BVCopyColumn_Svec_HIP(BV,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode BVResize_Svec_HIP(BV,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode BVGetColumn_Svec_HIP(BV,PetscInt,Vec*);
SLEPC_INTERN PetscErrorCode BVRestoreColumn_Svec_HIP(BV,PetscInt,Vec*);
SLEPC_INTERN PetscErrorCode BVRestoreSplit_Svec_HIP(BV,BV*,BV*);
SLEPC_INTERN PetscErrorCode BVRestoreSplitRows_Svec_HIP(BV,IS,IS,BV*,BV*);
SLEPC_INTERN PetscErrorCode BVGetMat_Svec_HIP(BV,Mat*);
SLEPC_INTERN PetscErrorCode BVRestoreMat_Svec_HIP(BV,Mat*);
#endif
