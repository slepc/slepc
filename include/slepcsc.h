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

/* SUBMANSEC = Sys */

/*S
   SlepcEigenvalueComparisonFn - A prototype of an eigenvalue comparison function that would
   be passed to `EPSSetEigenvalueComparison()` and analogue functions in other solver types.

   Calling Sequence:
+  ar     - real part of the 1st eigenvalue
.  ai     - imaginary part of the 1st eigenvalue
.  br     - real part of the 2nd eigenvalue
.  bi     - imaginary part of the 2nd eigenvalue
.  res    - [output] result of comparison
-  ctx    - [optional] user-defined context for private data for the
            eigenvalue comparison routine (may be `NULL`)

   Note:
   The return value `res` can be
+  negative - if the 1st value is preferred to the 2st one
.  zero     - if both values are equally preferred
-  positive - if the 2st value is preferred to the 1st one

   Level: advanced

.seealso: `SlepcSC`, `EPSSetEigenvalueComparison()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SlepcEigenvalueComparisonFn(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx);

/*S
   SlepcArbitrarySelectionFn - A prototype of an arbitrary selection function that would be
   passed to `EPSSetArbitrarySelection()` and analogue functions in other solver types.

   Calling Sequence:
+  er     - real part of the current eigenvalue approximation
.  ei     - imaginary part of the current eigenvalue approximation
.  xr     - real part of the current eigenvector approximation
.  xi     - imaginary part of the current eigenvector approximation
.  rr     - result of evaluation (real part)
.  ri     - result of evaluation (imaginary part)
-  ctx    - [optional] user-defined context for private data for the
            arbitrary selection routine (may be `NULL`)

   Notes:
   Given a scalar and a vector (usually an approximate eigenvalue and the corresponding
   eigenvector) this function evaluates the result `rr`, `ri` which is the value that
   will be used in a subsequent sorting.

   This evaluation function is collective, that is, all processes call it and
   it can use collective operations; furthermore, the computed result must
   be the same in all processes.

   The result is expressed as a complex number so that it is possible to
   use the standard eigenvalue sorting functions, but normally only `rr` is used.
   Set `ri` to zero unless it is meaningful in your application.

   Level: advanced

.seealso: `SlepcSC`, `SlepcEigenvalueComparisonFn`, `EPSSetArbitrarySelection()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode SlepcArbitrarySelectionFn(PetscScalar er,PetscScalar ei,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx);

/*S
   SlepcSC - Data structure (a C `struct`) for storing information about the sorting
   criterion used by different eigensolver objects.

   Notes:
   The `SlepcSC` structure contains a mapping function and a comparison
   function (with associated contexts).
   The mapping function usually calls the backtransform operation of an `ST` object.
   An optional region can also be used to give higher priority to values inside it.

   The member `comparison` typically points to a predefined comparison function
   such as `SlepcCompareLargestMagnitude`, though it can also be user-defined.

   Fortran Note:
   Fortran usage is not supported.

   Developer Notes:
   This is a low-level data structure common to all solver classes for the task of
   sorting a list of scalars, typically computed eigenvalues.

   The main operation associated to this data structure is `SlepcSCCompare()` to
   compare two scalar values, which in turn is called from `SlepcSortEigenvalues()`.

   Level: developer

.seealso: `SlepcEigenvalueComparisonFn`, `SlepcSortEigenvalues()`, `SlepcSCCompare()`
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
SLEPC_EXTERN PetscErrorCode SlepcSortEigenvalues(SlepcSC,PetscInt,PetscScalar[],PetscScalar[],PetscInt[]);
SLEPC_EXTERN PetscErrorCode SlepcSortEigenvaluesSpecial(SlepcSC,PetscInt,PetscScalar[],PetscScalar[],PetscInt[]);

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
