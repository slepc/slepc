/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Sorting criterion for various solvers
*/

#pragma once

#include <slepcsctypes.h>

/* SUBMANSEC = Sys */

SLEPC_EXTERN PetscErrorCode SlepcSCCompare(SlepcSC,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*);
SLEPC_EXTERN PetscErrorCode SlepcSortEigenvalues(SlepcSC,PetscInt,PetscScalar[],PetscScalar[],PetscInt[]);
SLEPC_EXTERN PetscErrorCode SlepcSortEigenvaluesSpecial(SlepcSC,PetscInt,PetscScalar[],PetscScalar[],PetscInt[]);

SLEPC_EXTERN PetscErrorCode SlepcMap_ST(PetscObject,PetscInt,PetscScalar*,PetscScalar*);

SLEPC_EXTERN PetscErrorCode SlepcCompareLargestMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareLargestReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareLargestImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestPosReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscCtx);
