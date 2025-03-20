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

#include <slepcsys.h>
#include <slepcrgtypes.h>

/* SUBMANSEC = sys */

/*S
  SlepcEigenvalueComparisonFn - A prototype of an eigenvalue comparison function that would be passed to EPSSetEigenvalueComparison() and analogue functions in other solver types

  Calling Sequence:
+   ar     - real part of the 1st eigenvalue
.   ai     - imaginary part of the 1st eigenvalue
.   br     - real part of the 2nd eigenvalue
.   bi     - imaginary part of the 2nd eigenvalue
.   res    - [output] result of comparison
-   ctx    - [optional] user-defined context for private data for the
             eigenvalue comparison routine (may be NULL)

  Level: advanced

.seealso: EPSSetEigenvalueComparison()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(SlepcEigenvalueComparisonFn)(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx);

/*S
  SlepcArbitrarySelectionFn - A prototype of an arbitrary selection function that would be passed to EPSSetArbitrarySelection() and analogue functions in other solver types

  Calling Sequence:
+   er     - real part of the current eigenvalue approximation
.   ei     - imaginary part of the current eigenvalue approximation
.   xr     - real part of the current eigenvector approximation
.   xi     - imaginary part of the current eigenvector approximation
.   rr     - result of evaluation (real part)
.   ri     - result of evaluation (imaginary part)
-   ctx    - [optional] user-defined context for private data for the
             arbitrary selection routine (may be NULL)

  Level: advanced

.seealso: EPSSetArbitrarySelection()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(SlepcArbitrarySelectionFn)(PetscScalar er,PetscScalar ei,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx);

/*S
    SlepcSC - Data structure (C struct) for storing information about
        the sorting criterion used by different eigensolver objects.

   Notes:
   The SlepcSC structure contains a mapping function and a comparison
   function (with associated contexts).
   The mapping function usually calls ST's backtransform.
   An optional region can also be used to give higher priority to values inside it.

   The comparison function must have the following calling sequence

$  comparison(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)

+  ar  - real part of the 1st eigenvalue
.  ai  - imaginary part of the 1st eigenvalue
.  br  - real part of the 2nd eigenvalue
.  bi  - imaginary part of the 2nd eigenvalue
.  res - result of comparison
-  ctx - optional context, stored in comparisonctx

   The returning parameter 'res' can be
+  negative - if the 1st value is preferred to the 2st one
.  zero     - if both values are equally preferred
-  positive - if the 2st value is preferred to the 1st one

   Fortran usage is not supported.

   Level: developer

.seealso: SlepcSCCompare()
S*/
struct _n_SlepcSC {
  /* map values before sorting, typically a call to STBackTransform (mapctx=ST) */
  PetscErrorCode (*map)(PetscObject,PetscInt,PetscScalar*,PetscScalar*);
  PetscObject    mapobj;
  /* comparison function such as SlepcCompareLargestMagnitude */
  SlepcEigenvalueComparisonFn *comparison;
  void           *comparisonctx;
  /* optional region for filtering */
  RG             rg;
};
typedef struct _n_SlepcSC* SlepcSC;

SLEPC_EXTERN PetscErrorCode SlepcSCCompare(SlepcSC,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*);
SLEPC_EXTERN PetscErrorCode SlepcSortEigenvalues(SlepcSC,PetscInt n,PetscScalar *eigr,PetscScalar *eigi,PetscInt *perm);

SLEPC_EXTERN PetscErrorCode SlepcMap_ST(PetscObject,PetscInt,PetscScalar*,PetscScalar*);

SLEPC_EXTERN PetscErrorCode SlepcCompareLargestMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareLargestReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareLargestImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetMagnitude(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareTargetImaginary(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
SLEPC_EXTERN PetscErrorCode SlepcCompareSmallestPosReal(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*);
