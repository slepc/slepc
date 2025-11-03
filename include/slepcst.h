/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Spectral transformation module for eigenvalue problems
*/

#pragma once

#include <slepcsys.h>
#include <slepcbv.h>
#include <petscksp.h>

/* SUBMANSEC = ST */

SLEPC_EXTERN PetscErrorCode STInitializePackage(void);
SLEPC_EXTERN PetscErrorCode STFinalizePackage(void);

/*S
   ST - Spectral transformation, encapsulates the functionality required
   for acceleration techniques based on the transformation of the spectrum,
   e.g., shift-and-invert.

   Level: beginner

.seealso: [](ch:st), `STCreate()`
S*/
typedef struct _p_ST* ST;

/*J
   STType - String with the name of the spectral transformation type.

   Level: beginner

.seealso: [](ch:st), `ST`, `STSetType()`
J*/
typedef const char *STType;
#define STSHIFT     "shift"
#define STSINVERT   "sinvert"
#define STCAYLEY    "cayley"
#define STPRECOND   "precond"
#define STFILTER    "filter"
#define STSHELL     "shell"

/* Logging support */
SLEPC_EXTERN PetscClassId ST_CLASSID;

SLEPC_EXTERN PetscErrorCode STCreate(MPI_Comm,ST*);
SLEPC_EXTERN PetscErrorCode STDestroy(ST*);
SLEPC_EXTERN PetscErrorCode STReset(ST);
SLEPC_EXTERN PetscErrorCode STSetType(ST,STType);
SLEPC_EXTERN PetscErrorCode STGetType(ST,STType*);
SLEPC_EXTERN PetscErrorCode STSetMatrices(ST,PetscInt,Mat[]);
SLEPC_EXTERN PetscErrorCode STGetMatrix(ST,PetscInt,Mat*);
SLEPC_EXTERN PetscErrorCode STGetMatrixTransformed(ST,PetscInt,Mat*);
SLEPC_EXTERN PetscErrorCode STGetNumMatrices(ST,PetscInt*);
SLEPC_EXTERN PetscErrorCode STGetOperator(ST,Mat*);
SLEPC_EXTERN PetscErrorCode STRestoreOperator(ST,Mat*);
SLEPC_EXTERN PetscErrorCode STSetUp(ST);
SLEPC_EXTERN PetscErrorCode STSetFromOptions(ST);
SLEPC_EXTERN PetscErrorCode STView(ST,PetscViewer);
SLEPC_EXTERN PetscErrorCode STViewFromOptions(ST,PetscObject,const char[]);

PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STSetMatrices()", ) static inline PetscErrorCode STSetOperators(ST st,PetscInt n,Mat *A) {return STSetMatrices(st,n,A);}
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STGetMatrix()", ) static inline PetscErrorCode STGetOperators(ST st,PetscInt k,Mat *A) {return STGetMatrix(st,k,A);}
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STGetMatrixTransformed()", ) static inline PetscErrorCode STGetTOperators(ST st,PetscInt k,Mat *A) {return STGetMatrixTransformed(st,k,A);}
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STGetOperator() followed by MatComputeOperator()", ) static inline PetscErrorCode STComputeExplicitOperator(ST st,Mat *A)
{
  Mat Op;

  PetscFunctionBegin;
  PetscCall(STGetOperator(st,&Op));
  PetscCall(MatComputeOperator(Op,MATAIJ,A));
  PetscCall(STRestoreOperator(st,&Op));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode STApply(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STApplyMat(ST,Mat,Mat);
SLEPC_EXTERN PetscErrorCode STApplyTranspose(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STApplyHermitianTranspose(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatMult(ST,PetscInt,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatMultTranspose(ST,PetscInt,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatMultHermitianTranspose(ST,PetscInt,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatSolve(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatSolveTranspose(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatSolveHermitianTranspose(ST,Vec,Vec);
SLEPC_EXTERN PetscErrorCode STMatMatSolve(ST,Mat,Mat);
SLEPC_EXTERN PetscErrorCode STGetBilinearForm(ST,Mat*);
SLEPC_EXTERN PetscErrorCode STMatSetUp(ST,PetscScalar,PetscScalar[]);
SLEPC_EXTERN PetscErrorCode STPostSolve(ST);
SLEPC_EXTERN PetscErrorCode STResetMatrixState(ST);
SLEPC_EXTERN PetscErrorCode STSetWorkVecs(ST,PetscInt);

SLEPC_EXTERN PetscErrorCode STSetKSP(ST,KSP);
SLEPC_EXTERN PetscErrorCode STGetKSP(ST,KSP*);
SLEPC_EXTERN PetscErrorCode STSetShift(ST,PetscScalar);
SLEPC_EXTERN PetscErrorCode STGetShift(ST,PetscScalar*);
SLEPC_EXTERN PetscErrorCode STSetDefaultShift(ST,PetscScalar);
SLEPC_EXTERN PetscErrorCode STScaleShift(ST,PetscScalar);
SLEPC_EXTERN PetscErrorCode STSetBalanceMatrix(ST,Vec);
SLEPC_EXTERN PetscErrorCode STGetBalanceMatrix(ST,Vec*);
SLEPC_EXTERN PetscErrorCode STSetTransform(ST,PetscBool);
SLEPC_EXTERN PetscErrorCode STGetTransform(ST,PetscBool*);
SLEPC_EXTERN PetscErrorCode STSetStructured(ST,PetscBool);
SLEPC_EXTERN PetscErrorCode STGetStructured(ST,PetscBool*);

SLEPC_EXTERN PetscErrorCode STSetOptionsPrefix(ST,const char[]);
SLEPC_EXTERN PetscErrorCode STAppendOptionsPrefix(ST,const char[]);
SLEPC_EXTERN PetscErrorCode STGetOptionsPrefix(ST,const char*[]);

SLEPC_EXTERN PetscErrorCode STBackTransform(ST,PetscInt,PetscScalar[],PetscScalar[]);
SLEPC_EXTERN PetscErrorCode STIsInjective(ST,PetscBool*);

SLEPC_EXTERN PetscErrorCode STCheckNullSpace(ST,BV);

SLEPC_EXTERN PetscErrorCode STSetPreconditionerMat(ST,Mat);
SLEPC_EXTERN PetscErrorCode STGetPreconditionerMat(ST,Mat*);
SLEPC_EXTERN PetscErrorCode STSetSplitPreconditioner(ST,PetscInt,Mat[],MatStructure);
SLEPC_EXTERN PetscErrorCode STGetSplitPreconditionerTerm(ST,PetscInt,Mat*);
SLEPC_EXTERN PetscErrorCode STGetSplitPreconditionerInfo(ST,PetscInt*,MatStructure*);

SLEPC_EXTERN PetscErrorCode STMatCreateVecs(ST,Vec*,Vec*);
SLEPC_EXTERN PetscErrorCode STMatCreateVecsEmpty(ST,Vec*,Vec*);
SLEPC_EXTERN PetscErrorCode STMatGetSize(ST,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode STMatGetLocalSize(ST,PetscInt*,PetscInt*);

/*E
   STMatMode - Determines how to handle the coefficient matrix of the linear
   system associated with the spectral transformation.

   Values:
+  `ST_MATMODE_COPY`    - the coefficient matrix is built explicitly on a copy of $A$
.  `ST_MATMODE_INPLACE` - the coefficient matrix is built explicitly overwriting $A$
-  `ST_MATMODE_SHELL`   - the coefficient matrix is handled implicitly

   Level: intermediate

.seealso: [](ch:st), `STSetMatMode()`, `STGetMatMode()`
E*/
typedef enum { ST_MATMODE_COPY,
               ST_MATMODE_INPLACE,
               ST_MATMODE_SHELL } STMatMode;
SLEPC_EXTERN const char *STMatModes[];

/*MC
   ST_MATMODE_COPY - The coefficient matrix of the linear system, $A-\sigma B$, is
   built explicitly on a copy of $A$.

   Note:
   If memory is an issue, one may prefer one of the other two modes.

   Level: intermediate

.seealso: [](ch:st), `STMatMode`, `STSetMatMode()`, `ST_MATMODE_INPLACE`, `ST_MATMODE_SHELL`
M*/

/*MC
   ST_MATMODE_INPLACE - The coefficient matrix of the linear system, $A-\sigma B$, is
   built explicitly overwritting $A$.

   Note:
   This mode uses less memory than `ST_MATMODE_COPY`, but it modifies $A$. This
   alteration of $A$ is reversed after the eigensolution process has finished,
   but due to roundoff, the result might not be exactly equal to the original $A$.

   Level: intermediate

.seealso: [](ch:st), `STMatMode`, `STSetMatMode()`, `ST_MATMODE_COPY`, `ST_MATMODE_SHELL`
M*/

/*MC
   ST_MATMODE_SHELL - The coefficient matrix of the linear system, $A-\sigma B$, is
   not built explicitly, and instead it is handled implicitly via a `MATSHELL`.

   Note:
   This mode severely restricts the number of possibilities available for solving
   the linear system via `KSP`.

   Level: intermediate

.seealso: [](ch:st), `STMatMode`, `STSetMatMode()`, `ST_MATMODE_COPY`, `ST_MATMODE_INPLACE`
M*/

SLEPC_EXTERN PetscErrorCode STSetMatMode(ST,STMatMode);
SLEPC_EXTERN PetscErrorCode STGetMatMode(ST,STMatMode*);
SLEPC_EXTERN PetscErrorCode STSetMatStructure(ST,MatStructure);
SLEPC_EXTERN PetscErrorCode STGetMatStructure(ST,MatStructure*);

SLEPC_EXTERN PetscFunctionList STList;
SLEPC_EXTERN PetscErrorCode STRegister(const char[],PetscErrorCode(*)(ST));

/* --------- options specific to particular spectral transformations-------- */

/*S
   STShellApplyFn - A prototype of a function for the user-defined `STApply()`
   operation in an `STSHELL`.

   Calling Sequence:
+  st   - the spectral transformation context
.  xin  - input vector
-  xout - output vector

   Level: advanced

.seealso: [](ch:st), `STShellSetApply()`, `STApply()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode STShellApplyFn(ST st,Vec xin,Vec xout);

/*S
   STShellApplyTransposeFn - A prototype of a function for the user-defined `STApplyTranspose()`
   operation in an `STSHELL`.

   Calling Sequence:
+  st   - the spectral transformation context
.  xin  - input vector
-  xout - output vector

   Level: advanced

.seealso: [](ch:st), `STShellSetApplyTranspose()`, `STApplyTranspose()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode STShellApplyTransposeFn(ST st,Vec xin,Vec xout);

/*S
   STShellApplyHermitianTransposeFn - A prototype of a function for the user-defined
   `STApplyHermitianTranspose()` operation in an `STSHELL`.

   Calling Sequence:
+  st   - the spectral transformation context
.  xin  - input vector
-  xout - output vector

   Level: advanced

.seealso: [](ch:st), `STShellSetApplyHermitianTranspose()`, `STApplyHermitianTranspose()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode STShellApplyHermitianTransposeFn(ST st,Vec xin,Vec xout);

/*S
   STShellBackTransformFn - A prototype of a function for the user-defined `STBackTransform()`
   operation in an `STSHELL`.

   Calling Sequence:
+  st   - the spectral transformation context
.  n    - number of eigenvalues to be backtransformed
.  eigr - pointer to the real parts of the eigenvalues to transform back
-  eigi - pointer to the imaginary parts

   Level: advanced

.seealso: [](ch:st), `STShellSetBackTransform()`, `STBackTransform()`
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode STShellBackTransformFn(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi);

SLEPC_EXTERN PetscErrorCode STShellGetContext(ST,void*);
SLEPC_EXTERN PetscErrorCode STShellSetContext(ST,void*);
SLEPC_EXTERN PetscErrorCode STShellSetApply(ST,STShellApplyFn*);
SLEPC_EXTERN PetscErrorCode STShellSetApplyTranspose(ST,STShellApplyTransposeFn*);
SLEPC_EXTERN PetscErrorCode STShellSetApplyHermitianTranspose(ST,STShellApplyHermitianTransposeFn*);
SLEPC_EXTERN PetscErrorCode STShellSetBackTransform(ST,STShellBackTransformFn*);

SLEPC_EXTERN PetscErrorCode STCayleyGetAntishift(ST,PetscScalar*);
SLEPC_EXTERN PetscErrorCode STCayleySetAntishift(ST,PetscScalar);

PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STGetPreconditionerMat()", ) static inline PetscErrorCode STPrecondGetMatForPC(ST st,Mat *A) {return STGetPreconditionerMat(st,A);}
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STSetPreconditionerMat()", ) static inline PetscErrorCode STPrecondSetMatForPC(ST st,Mat A) {return STSetPreconditionerMat(st,A);}
SLEPC_EXTERN PetscErrorCode STPrecondGetKSPHasMat(ST,PetscBool*);
SLEPC_EXTERN PetscErrorCode STPrecondSetKSPHasMat(ST,PetscBool);

/*E
   STFilterType - Selects the method used to build the filter.

   Values:
+  `ST_FILTER_FILTLAN`   - the filter is built with FILTLAN (Filtered Lanczos)
-  `ST_FILTER_CHEBYSHEV` - the filter is built with a Chebyshev series

   Level: intermediate

.seealso: [](ch:st), `STFilterSetType()`, `STGetFilterType()`
E*/
typedef enum { ST_FILTER_FILTLAN   = 1,
               ST_FILTER_CHEBYSHEV = 2 } STFilterType;
SLEPC_EXTERN const char *STFilterTypes[];

/*MC
   ST_FILTER_FILTLAN - The polynomial filter is built with FILTLAN (Filtered Lanczos).

   Note:
   This filter implements the Filtered Lanczos method {cite:p}`Fan12`. In fact,
   the implementation adapts files from the FILTLAN package by the same authors,
   which have been converted for a native integration in SLEPc.

   Level: intermediate

.seealso: [](ch:st), `STFilterType`, `STFilterSetType()`, `ST_FILTER_CHEBYSHEV`
M*/

/*MC
   ST_FILTER_CHEBYSHEV - The filter is based on a truncated Chebyshev series.

   Note:
   The application of the filter is implemented by means of a Clenshaw
   algorithm. This filter also supports using damping, see `STFilterSetDamping()`.

   Level: intermediate

.seealso: [](ch:st), `STFilterType`, `STSetFilterType()`, `STFilterSetDamping()`, `ST_FILTER_FILTLAN`
M*/

/*E
   STFilterDamping - The damping type used to build the filter.

   Values:
+  `ST_FILTER_DAMPING_NONE`    - no damping
.  `ST_FILTER_DAMPING_JACKSON` - Jackson damping
.  `ST_FILTER_DAMPING_LANCZOS` - Lanczos damping
-  `ST_FILTER_DAMPING_FEJER`   - Fejer damping

   Notes:
   The default is no damping.

   For the definition of the damping coefficients, see for instance {cite:p}`Pie16`

   Level: advanced

.seealso: [](ch:st), `STSetFilterDamping()`, `STGetFilterDamping()`
E*/
typedef enum { ST_FILTER_DAMPING_NONE,
               ST_FILTER_DAMPING_JACKSON,
               ST_FILTER_DAMPING_LANCZOS,
               ST_FILTER_DAMPING_FEJER } STFilterDamping;
SLEPC_EXTERN const char *STFilterDampings[];

/*MC
   ST_MATMODE_COPY - The coefficient matrix of the linear system, $A-\sigma B$, is
   built explicitly on a copy of $A$.

   Note:
   If memory is an issue, one may prefer one of the other two modes.

   Level: intermediate

.seealso: [](ch:st), `STMatMode`, `STSetMatMode()`, `ST_MATMODE_INPLACE`, `ST_MATMODE_SHELL`
M*/

SLEPC_EXTERN PetscErrorCode STFilterSetType(ST,STFilterType);
SLEPC_EXTERN PetscErrorCode STFilterGetType(ST,STFilterType*);
SLEPC_EXTERN PetscErrorCode STFilterSetInterval(ST,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode STFilterGetInterval(ST,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode STFilterSetRange(ST,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode STFilterGetRange(ST,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode STFilterSetDegree(ST,PetscInt);
SLEPC_EXTERN PetscErrorCode STFilterGetDegree(ST,PetscInt*);
SLEPC_EXTERN PetscErrorCode STFilterGetThreshold(ST,PetscReal*);
SLEPC_EXTERN PetscErrorCode STFilterSetDamping(ST,STFilterDamping);
SLEPC_EXTERN PetscErrorCode STFilterGetDamping(ST,STFilterDamping*);
