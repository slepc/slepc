/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the mathematical function object in SLEPc
*/

#pragma once

#include <slepcsys.h>

/* SUBMANSEC = FN */

#define FN_MAX_SOLVE 16

SLEPC_EXTERN PetscErrorCode FNInitializePackage(void);
SLEPC_EXTERN PetscErrorCode FNFinalizePackage(void);

/*S
   FN - Abstraction of a mathematical function.

   Level: beginner

.seealso: [](sec:fn), `FNCreate()`
S*/
typedef struct _p_FN* FN;

/*J
   FNType - String with the name of the mathematical function.

   Level: beginner

.seealso: [](sec:fn), `FNSetType()`, `FN`
J*/
typedef const char *FNType;
#define FNRATIONAL "rational"
#define FNEXP      "exp"
#define FNLOG      "log"
#define FNPHI      "phi"
#define FNSQRT     "sqrt"
#define FNINVSQRT  "invsqrt"
#define FNCOMBINE  "combine"

/* Logging support */
SLEPC_EXTERN PetscClassId FN_CLASSID;

/*E
   FNCombineType - Determines how two functions are combined.

   Values:
+  `FN_COMBINE_ADD`      - add the functions
.  `FN_COMBINE_MULTIPLY` - multiply the functions
.  `FN_COMBINE_DIVIDE`   - divide the functions
-  `FN_COMBINE_COMPOSE`  - use function composition

   Level: intermediate

.seealso: [](sec:fn), `FNCombineSetChildren()`
E*/
typedef enum { FN_COMBINE_ADD,
               FN_COMBINE_MULTIPLY,
               FN_COMBINE_DIVIDE,
               FN_COMBINE_COMPOSE } FNCombineType;

/*MC
   FN_COMBINE_ADD - In functions of type `FNCOMBINE`, add the two child functions
   together.

   Level: intermediate

.seealso: [](sec:fn), `FNCombineType`, `FNCombineSetChildren()`, `FN_COMBINE_MULTIPLY`, `FN_COMBINE_DIVIDE`, `FN_COMBINE_COMPOSE`
M*/

/*MC
   FN_COMBINE_MULTIPLY - In functions of type `FNCOMBINE`, multiply the two child functions
   together.

   Level: intermediate

.seealso: [](sec:fn), `FNCombineType`, `FNCombineSetChildren()`, `FN_COMBINE_ADD`, `FN_COMBINE_DIVIDE`, `FN_COMBINE_COMPOSE`
M*/

/*MC
   FN_COMBINE_DIVIDE - In functions of type `FNCOMBINE`, compute the ratio of the
   two child functions (the first one in the numerator).

   Level: intermediate

.seealso: [](sec:fn), `FNCombineType`, `FNCombineSetChildren()`, `FN_COMBINE_ADD`, `FN_COMBINE_MULTIPLY`, `FN_COMBINE_COMPOSE`
M*/

/*MC
   FN_COMBINE_COMPOSE - In functions of type `FNCOMBINE`, compose the two child functions,
   i.e., evaluate the second function on the result of evaluating the first function.

   Level: intermediate

.seealso: [](sec:fn), `FNCombineType`, `FNCombineSetChildren()`, `FN_COMBINE_ADD`, `FN_COMBINE_MULTIPLY`, `FN_COMBINE_DIVIDE`
M*/

/*E
   FNParallelType - Indicates the parallel mode that will be used for matrix
   function evaluation.

   Values:
+  `FN_PARALLEL_REDUNDANT`    - all processes compute redundantly
-  `FN_PARALLEL_SYNCHRONIZED` - only the first MPI process performs the computation

   Level: advanced

.seealso: [](sec:fn), `FNSetParallel()`
E*/
typedef enum { FN_PARALLEL_REDUNDANT,
               FN_PARALLEL_SYNCHRONIZED } FNParallelType;
SLEPC_EXTERN const char *FNParallelTypes[];

/*MC
   FN_PARALLEL_REDUNDANT - In matrix function evaluation, all processes compute
   redundantly.

   Note:
   When this parallel mode is selected, all processes will make the computation
   redundantly, starting from the same data, and producing the same result.
   This result may be slightly different in the different processes if using a
   multithreaded BLAS library, which may cause issues in ill-conditioned problems.

   Level: advanced

.seealso: [](sec:fn), `FNParallelType`, `FNSetParallel()`, `FN_PARALLEL_SYNCHRONIZED`
M*/

/*MC
   FN_PARALLEL_SYNCHRONIZED - In matrix function evaluation, only the first MPI
   process performs the computation.

   Note:
   When this parallel mode is selected, only the first MPI process performs the
   computation and then the computed matrix is broadcast to the other
   processes in the communicator. This communication is done automatically at
   the end of `FNEvaluateFunctionMat()` or `FNEvaluateFunctionMatVec()`.

   Level: advanced

.seealso: [](sec:fn), `FNParallelType`, `FNSetParallel()`, `FN_PARALLEL_REDUNDANT`
M*/

SLEPC_EXTERN PetscErrorCode FNCreate(MPI_Comm,FN*);
SLEPC_EXTERN PetscErrorCode FNSetType(FN,FNType);
SLEPC_EXTERN PetscErrorCode FNGetType(FN,FNType*);
SLEPC_EXTERN PetscErrorCode FNSetOptionsPrefix(FN,const char[]);
SLEPC_EXTERN PetscErrorCode FNAppendOptionsPrefix(FN,const char[]);
SLEPC_EXTERN PetscErrorCode FNGetOptionsPrefix(FN,const char*[]);
SLEPC_EXTERN PetscErrorCode FNSetFromOptions(FN);
SLEPC_EXTERN PetscErrorCode FNView(FN,PetscViewer);
SLEPC_EXTERN PetscErrorCode FNViewFromOptions(FN,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode FNDestroy(FN*);
SLEPC_EXTERN PetscErrorCode FNDuplicate(FN,MPI_Comm,FN*);

SLEPC_EXTERN PetscErrorCode FNSetScale(FN,PetscScalar,PetscScalar);
SLEPC_EXTERN PetscErrorCode FNGetScale(FN,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode FNSetMethod(FN,PetscInt);
SLEPC_EXTERN PetscErrorCode FNGetMethod(FN,PetscInt*);
SLEPC_EXTERN PetscErrorCode FNSetParallel(FN,FNParallelType);
SLEPC_EXTERN PetscErrorCode FNGetParallel(FN,FNParallelType*);

SLEPC_EXTERN PetscErrorCode FNEvaluateFunction(FN,PetscScalar,PetscScalar*);
SLEPC_EXTERN PetscErrorCode FNEvaluateDerivative(FN,PetscScalar,PetscScalar*);
SLEPC_EXTERN PetscErrorCode FNEvaluateFunctionMat(FN,Mat,Mat);
SLEPC_EXTERN PetscErrorCode FNEvaluateFunctionMatVec(FN,Mat,Vec);

SLEPC_EXTERN PetscFunctionList FNList;
SLEPC_EXTERN PetscErrorCode FNRegister(const char[],PetscErrorCode(*)(FN));

/* --------- options specific to particular functions -------- */

SLEPC_EXTERN PetscErrorCode FNRationalSetNumerator(FN,PetscInt,PetscScalar[]);
SLEPC_EXTERN PetscErrorCode FNRationalGetNumerator(FN,PetscInt*,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode FNRationalSetDenominator(FN,PetscInt,PetscScalar[]);
SLEPC_EXTERN PetscErrorCode FNRationalGetDenominator(FN,PetscInt*,PetscScalar*[]);

SLEPC_EXTERN PetscErrorCode FNCombineSetChildren(FN,FNCombineType,FN,FN);
SLEPC_EXTERN PetscErrorCode FNCombineGetChildren(FN,FNCombineType*,FN*,FN*);

SLEPC_EXTERN PetscErrorCode FNPhiSetIndex(FN,PetscInt);
SLEPC_EXTERN PetscErrorCode FNPhiGetIndex(FN,PetscInt*);
