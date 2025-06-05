/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc object for solving linear matrix equations
*/

#pragma once

#include <slepcbv.h>

/* SUBMANSEC = LME */

SLEPC_EXTERN PetscErrorCode LMEInitializePackage(void);
SLEPC_EXTERN PetscErrorCode LMEFinalizePackage(void);

/*S
    LME - SLEPc object that encapsulates functionality for linear matrix equations

    Level: beginner

.seealso:  LMECreate()
S*/
typedef struct _p_LME* LME;

/*J
    LMEType - String with the name of a method for solving linear matrix equations

    Level: beginner

.seealso: LMESetType(), LME
J*/
typedef const char *LMEType;
#define LMEKRYLOV   "krylov"

/* Logging support */
SLEPC_EXTERN PetscClassId LME_CLASSID;

/*E
    LMEProblemType - Determines the type of linear matrix equation

    Level: beginner

.seealso: LMESetProblemType(), LMEGetProblemType()
E*/
typedef enum { LME_LYAPUNOV,
               LME_SYLVESTER,
               LME_GEN_LYAPUNOV,
               LME_GEN_SYLVESTER,
               LME_DT_LYAPUNOV ,
               LME_STEIN} LMEProblemType;
SLEPC_EXTERN const char *LMEProblemTypes[];

SLEPC_EXTERN PetscErrorCode LMECreate(MPI_Comm,LME *);
SLEPC_EXTERN PetscErrorCode LMEDestroy(LME*);
SLEPC_EXTERN PetscErrorCode LMEReset(LME);
SLEPC_EXTERN PetscErrorCode LMESetType(LME,LMEType);
SLEPC_EXTERN PetscErrorCode LMEGetType(LME,LMEType*);
SLEPC_EXTERN PetscErrorCode LMESetProblemType(LME,LMEProblemType);
SLEPC_EXTERN PetscErrorCode LMEGetProblemType(LME,LMEProblemType*);
SLEPC_EXTERN PetscErrorCode LMESetCoefficients(LME,Mat,Mat,Mat,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetCoefficients(LME,Mat*,Mat*,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetRHS(LME,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetRHS(LME,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetSolution(LME,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetSolution(LME,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetFromOptions(LME);
SLEPC_EXTERN PetscErrorCode LMESetUp(LME);
SLEPC_EXTERN PetscErrorCode LMESolve(LME);
SLEPC_EXTERN PetscErrorCode LMEView(LME,PetscViewer);
SLEPC_EXTERN PetscErrorCode LMEViewFromOptions(LME,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode LMEConvergedReasonView(LME,PetscViewer);
SLEPC_EXTERN PetscErrorCode LMEConvergedReasonViewFromOptions(LME);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "LMEConvergedReasonView()", ) static inline PetscErrorCode LMEReasonView(LME lme,PetscViewer v) {return LMEConvergedReasonView(lme,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "LMEConvergedReasonViewFromOptions()", ) static inline PetscErrorCode LMEReasonViewFromOptions(LME lme) {return LMEConvergedReasonViewFromOptions(lme);}

SLEPC_EXTERN PetscErrorCode LMESetBV(LME,BV);
SLEPC_EXTERN PetscErrorCode LMEGetBV(LME,BV*);
SLEPC_EXTERN PetscErrorCode LMESetTolerances(LME,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEGetTolerances(LME,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode LMESetDimensions(LME,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEGetDimensions(LME,PetscInt*);
SLEPC_EXTERN PetscErrorCode LMEGetIterationNumber(LME,PetscInt*);

SLEPC_EXTERN PetscErrorCode LMEGetErrorEstimate(LME,PetscReal*);
SLEPC_EXTERN PetscErrorCode LMEComputeError(LME,PetscReal*);
SLEPC_EXTERN PetscErrorCode LMESetErrorIfNotConverged(LME,PetscBool);
SLEPC_EXTERN PetscErrorCode LMEGetErrorIfNotConverged(LME,PetscBool*);

SLEPC_EXTERN PetscErrorCode LMEDenseLyapunov(LME,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEDenseHessLyapunovChol(LME,PetscInt,PetscScalar*,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscReal*);
PETSC_DEPRECATED_FUNCTION(3, 8, 0, "LMEDenseHessLyapunovChol()", ) static inline PetscErrorCode LMEDenseLyapunovChol(LME lme,PetscScalar *H,PetscInt m,PetscInt ldh,PetscScalar *r,PetscScalar *L,PetscInt ldl,PetscReal *res) {return LMEDenseHessLyapunovChol(lme,m,H,ldh,1,r,m,L,ldl,res);}

/*S
  LMEMonitorFn - A function prototype for functions provided to LMEMonitorSet()

  Calling Sequence:
+   lme    - linear matrix equation solver context obtained from LMECreate()
.   its    - iteration number
.   errest - error estimate
-   ctx    - optional monitoring context, as provided with LMEMonitorSet()

  Level: beginner

.seealso: LMEMonitorSet()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode LMEMonitorFn(LME lme,PetscInt its,PetscReal errest,void *ctx);

/*S
  LMEMonitorRegisterFn - A function prototype for functions provided to LMEMonitorRegister()

  Calling Sequence:
+   lme    - linear matrix equation solver context obtained from LMECreate()
.   its    - iteration number
.   errest - error estimate
-   ctx    - PetscViewerAndFormat object

  Level: beginner

  Note:
  This is an LMEMonitorFn specialized for a context of PetscViewerAndFormat.

.seealso: LMEMonitorSet(), LMEMonitorRegister(), LMEMonitorFn, LMEMonitorRegisterCreateFn, LMEMonitorRegisterDestroyFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode LMEMonitorRegisterFn(LME lme,PetscInt its,PetscReal errest,PetscViewerAndFormat *ctx);

/*S
  LMEMonitorRegisterCreateFn - A function prototype for functions that do the creation when provided to LMEMonitorRegister()

  Calling Sequence:
+   viewer - the viewer to be used with the LMEMonitorRegisterFn
.   format - the format of the viewer
.   ctx    - a context for the monitor
-   result - a PetscViewerAndFormat object

  Level: beginner

.seealso: LMEMonitorRegisterFn, LMEMonitorSet(), LMEMonitorRegister(), LMEMonitorFn, LMEMonitorRegisterDestroyFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode LMEMonitorRegisterCreateFn(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **result);

/*S
  LMEMonitorRegisterDestroyFn - A function prototype for functions that do the after use destruction when provided to LMEMonitorRegister()

  Calling Sequence:
.   vf - a PetscViewerAndFormat object to be destroyed, including any context

  Level: beginner

.seealso: LMEMonitorRegisterFn, LMEMonitorSet(), LMEMonitorRegister(), LMEMonitorFn, LMEMonitorRegisterCreateFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode LMEMonitorRegisterDestroyFn(PetscViewerAndFormat **result);

SLEPC_EXTERN PetscErrorCode LMEMonitor(LME,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode LMEMonitorSet(LME,LMEMonitorFn,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode LMEMonitorCancel(LME);
SLEPC_EXTERN PetscErrorCode LMEGetMonitorContext(LME,void*);

SLEPC_EXTERN PetscErrorCode LMEMonitorSetFromOptions(LME,const char[],const char[],void*);
SLEPC_EXTERN LMEMonitorRegisterFn       LMEMonitorDefault;
SLEPC_EXTERN LMEMonitorRegisterFn       LMEMonitorDefaultDrawLG;
SLEPC_EXTERN LMEMonitorRegisterCreateFn LMEMonitorDefaultDrawLGCreate;

SLEPC_EXTERN PetscErrorCode LMESetOptionsPrefix(LME,const char*);
SLEPC_EXTERN PetscErrorCode LMEAppendOptionsPrefix(LME,const char*);
SLEPC_EXTERN PetscErrorCode LMEGetOptionsPrefix(LME,const char*[]);

/*E
    LMEConvergedReason - reason a matrix function iteration was said to
         have converged or diverged

    Level: intermediate

.seealso: LMESolve(), LMEGetConvergedReason(), LMESetTolerances()
E*/
typedef enum {/* converged */
              LME_CONVERGED_TOL                =  1,
              /* diverged */
              LME_DIVERGED_ITS                 = -1,
              LME_DIVERGED_BREAKDOWN           = -2,
              LME_CONVERGED_ITERATING          =  0} LMEConvergedReason;
SLEPC_EXTERN const char *const*LMEConvergedReasons;

SLEPC_EXTERN PetscErrorCode LMEGetConvergedReason(LME,LMEConvergedReason *);

SLEPC_EXTERN PetscFunctionList LMEList;
SLEPC_EXTERN PetscFunctionList LMEMonitorList;
SLEPC_EXTERN PetscFunctionList LMEMonitorCreateList;
SLEPC_EXTERN PetscFunctionList LMEMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode LMERegister(const char[],PetscErrorCode(*)(LME));
SLEPC_EXTERN PetscErrorCode LMEMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,LMEMonitorRegisterFn*,LMEMonitorRegisterCreateFn*,LMEMonitorRegisterDestroyFn*);

SLEPC_EXTERN PetscErrorCode LMEAllocateSolution(LME,PetscInt);
