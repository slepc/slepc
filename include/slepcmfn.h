/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc matrix function solver object
*/

#pragma once

#include <slepcbv.h>
#include <slepcfn.h>

/* SUBMANSEC = MFN */

SLEPC_EXTERN PetscErrorCode MFNInitializePackage(void);
SLEPC_EXTERN PetscErrorCode MFNFinalizePackage(void);

/*S
    MFN - SLEPc object that encapsulates functionality for matrix functions.

    Level: beginner

.seealso:  MFNCreate()
S*/
typedef struct _p_MFN* MFN;

/*J
    MFNType - String with the name of a method for computing matrix functions.

    Level: beginner

.seealso: MFNSetType(), MFN
J*/
typedef const char *MFNType;
#define MFNKRYLOV   "krylov"
#define MFNEXPOKIT  "expokit"

/* Logging support */
SLEPC_EXTERN PetscClassId MFN_CLASSID;

SLEPC_EXTERN PetscErrorCode MFNCreate(MPI_Comm,MFN *);
SLEPC_EXTERN PetscErrorCode MFNDestroy(MFN*);
SLEPC_EXTERN PetscErrorCode MFNReset(MFN);
SLEPC_EXTERN PetscErrorCode MFNSetType(MFN,MFNType);
SLEPC_EXTERN PetscErrorCode MFNGetType(MFN,MFNType*);
SLEPC_EXTERN PetscErrorCode MFNSetOperator(MFN,Mat);
SLEPC_EXTERN PetscErrorCode MFNGetOperator(MFN,Mat*);
SLEPC_EXTERN PetscErrorCode MFNSetFromOptions(MFN);
SLEPC_EXTERN PetscErrorCode MFNSetUp(MFN);
SLEPC_EXTERN PetscErrorCode MFNSolve(MFN,Vec,Vec);
SLEPC_EXTERN PetscErrorCode MFNSolveTranspose(MFN,Vec,Vec);
SLEPC_EXTERN PetscErrorCode MFNView(MFN,PetscViewer);
SLEPC_EXTERN PetscErrorCode MFNViewFromOptions(MFN,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode MFNConvergedReasonView(MFN,PetscViewer);
SLEPC_EXTERN PetscErrorCode MFNConvergedReasonViewFromOptions(MFN);
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "MFNConvergedReasonView()", ) static inline PetscErrorCode MFNReasonView(MFN mfn,PetscViewer v) {return MFNConvergedReasonView(mfn,v);}
PETSC_DEPRECATED_FUNCTION(3, 14, 0, "MFNConvergedReasonViewFromOptions()", ) static inline PetscErrorCode MFNReasonViewFromOptions(MFN mfn) {return MFNConvergedReasonViewFromOptions(mfn);}

SLEPC_EXTERN PetscErrorCode MFNSetBV(MFN,BV);
SLEPC_EXTERN PetscErrorCode MFNGetBV(MFN,BV*);
SLEPC_EXTERN PetscErrorCode MFNSetFN(MFN,FN);
SLEPC_EXTERN PetscErrorCode MFNGetFN(MFN,FN*);
SLEPC_EXTERN PetscErrorCode MFNSetTolerances(MFN,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode MFNGetTolerances(MFN,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode MFNSetDimensions(MFN,PetscInt);
SLEPC_EXTERN PetscErrorCode MFNGetDimensions(MFN,PetscInt*);
SLEPC_EXTERN PetscErrorCode MFNGetIterationNumber(MFN,PetscInt*);

SLEPC_EXTERN PetscErrorCode MFNSetErrorIfNotConverged(MFN,PetscBool);
SLEPC_EXTERN PetscErrorCode MFNGetErrorIfNotConverged(MFN,PetscBool*);

/*S
  MFNMonitorFn - A function prototype for functions provided to MFNMonitorSet()

  Calling Sequence:
+   mfn    - matrix function context obtained from MFNCreate()
.   its    - iteration number
.   errest - error estimate
-   ctx    - optional monitoring context, as provided with MFNMonitorSet()

  Level: beginner

.seealso: MFNMonitorSet()
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode MFNMonitorFn(MFN mfn,PetscInt its,PetscReal errest,void *ctx);

/*S
  MFNMonitorRegisterFn - A function prototype for functions provided to MFNMonitorRegister()

  Calling Sequence:
+   mfn    - matrix function context obtained from MFNCreate()
.   its    - iteration number
.   errest - error estimate
-   ctx    - PetscViewerAndFormat object

  Level: beginner

  Note:
  This is an MFNMonitorFn specialized for a context of PetscViewerAndFormat.

.seealso: MFNMonitorSet(), MFNMonitorRegister(), MFNMonitorFn, MFNMonitorRegisterCreateFn, MFNMonitorRegisterDestroyFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode MFNMonitorRegisterFn(MFN mfn,PetscInt its,PetscReal errest,PetscViewerAndFormat *ctx);

/*S
  MFNMonitorRegisterCreateFn - A function prototype for functions that do the creation when provided to MFNMonitorRegister()

  Calling Sequence:
+   viewer - the viewer to be used with the MFNMonitorRegisterFn
.   format - the format of the viewer
.   ctx    - a context for the monitor
-   result - a PetscViewerAndFormat object

  Level: beginner

.seealso: MFNMonitorRegisterFn, MFNMonitorSet(), MFNMonitorRegister(), MFNMonitorFn, MFNMonitorRegisterDestroyFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode MFNMonitorRegisterCreateFn(PetscViewer viewer,PetscViewerFormat format,void *ctx,PetscViewerAndFormat **result);

/*S
  MFNMonitorRegisterDestroyFn - A function prototype for functions that do the after use destruction when provided to MFNMonitorRegister()

  Calling Sequence:
.   vf - a PetscViewerAndFormat object to be destroyed, including any context

  Level: beginner

.seealso: MFNMonitorRegisterFn, MFNMonitorSet(), MFNMonitorRegister(), MFNMonitorFn, MFNMonitorRegisterCreateFn
S*/
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode MFNMonitorRegisterDestroyFn(PetscViewerAndFormat **result);

SLEPC_EXTERN PetscErrorCode MFNMonitor(MFN,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode MFNMonitorSet(MFN,MFNMonitorFn,void*,PetscCtxDestroyFn*);
SLEPC_EXTERN PetscErrorCode MFNMonitorCancel(MFN);
SLEPC_EXTERN PetscErrorCode MFNGetMonitorContext(MFN,void*);

SLEPC_EXTERN PetscErrorCode MFNMonitorSetFromOptions(MFN,const char[],const char[],void*);
SLEPC_EXTERN MFNMonitorRegisterFn       MFNMonitorDefault;
SLEPC_EXTERN MFNMonitorRegisterFn       MFNMonitorDefaultDrawLG;
SLEPC_EXTERN MFNMonitorRegisterCreateFn MFNMonitorDefaultDrawLGCreate;

SLEPC_EXTERN PetscErrorCode MFNSetOptionsPrefix(MFN,const char*);
SLEPC_EXTERN PetscErrorCode MFNAppendOptionsPrefix(MFN,const char*);
SLEPC_EXTERN PetscErrorCode MFNGetOptionsPrefix(MFN,const char*[]);

/*E
    MFNConvergedReason - reason a matrix function iteration was said to
         have converged or diverged

    Level: intermediate

.seealso: MFNSolve(), MFNGetConvergedReason(), MFNSetTolerances()
E*/
typedef enum {/* converged */
              MFN_CONVERGED_TOL                =  1,
              MFN_CONVERGED_ITS                =  2,
              /* diverged */
              MFN_DIVERGED_ITS                 = -1,
              MFN_DIVERGED_BREAKDOWN           = -2,
              MFN_CONVERGED_ITERATING          =  0} MFNConvergedReason;
SLEPC_EXTERN const char *const*MFNConvergedReasons;

SLEPC_EXTERN PetscErrorCode MFNGetConvergedReason(MFN,MFNConvergedReason *);

SLEPC_EXTERN PetscFunctionList MFNList;
SLEPC_EXTERN PetscFunctionList MFNMonitorList;
SLEPC_EXTERN PetscFunctionList MFNMonitorCreateList;
SLEPC_EXTERN PetscFunctionList MFNMonitorDestroyList;
SLEPC_EXTERN PetscErrorCode MFNRegister(const char[],PetscErrorCode(*)(MFN));
SLEPC_EXTERN PetscErrorCode MFNMonitorRegister(const char[],PetscViewerType,PetscViewerFormat,MFNMonitorRegisterFn*,MFNMonitorRegisterCreateFn*,MFNMonitorRegisterDestroyFn*);

SLEPC_EXTERN PetscErrorCode MFNAllocateSolution(MFN,PetscInt);
