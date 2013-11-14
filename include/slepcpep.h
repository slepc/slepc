/*
   User interface for SLEPc's polynomial eigenvalue solvers.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2013, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.

   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(__SLEPCPEP_H)
#define __SLEPCPEP_H
#include <slepceps.h>

PETSC_EXTERN PetscErrorCode PEPInitializePackage(void);

/*S
     PEP - Abstract SLEPc object that manages all the polynomial eigenvalue
     problem solvers.

   Level: beginner

.seealso:  PEPCreate()
S*/
typedef struct _p_PEP* PEP;

/*J
    PEPType - String with the name of a polynomial eigensolver

   Level: beginner

.seealso: PEPSetType(), PEP
J*/
typedef const char* PEPType;
#define PEPLINEAR    "linear"
#define PEPPARNOLDI  "parnoldi"
#define PEPPLANCZOS  "planczos"
#define PEPTOAR      "toar"
#define PEPSTOAR     "stoar"

/* Logging support */
PETSC_EXTERN PetscClassId PEP_CLASSID;

/*E
    PEPProblemType - Determines the type of the polynomial eigenproblem

    Level: intermediate

.seealso: PEPSetProblemType(), PEPGetProblemType()
E*/
typedef enum { PEP_GENERAL=1,
               PEP_HERMITIAN,   /* All A_i  Hermitian */
               PEP_GYROSCOPIC   /* QEP with M, K  Hermitian, M>0, C skew-Hermitian */
             } PEPProblemType;

/*E
    PEPWhich - Determines which part of the spectrum is requested

    Level: intermediate

.seealso: PEPSetWhichEigenpairs(), PEPGetWhichEigenpairs()
E*/
typedef enum { PEP_LARGEST_MAGNITUDE=1,
               PEP_SMALLEST_MAGNITUDE,
               PEP_LARGEST_REAL,
               PEP_SMALLEST_REAL,
               PEP_LARGEST_IMAGINARY,
               PEP_SMALLEST_IMAGINARY,
               PEP_TARGET_MAGNITUDE,
               PEP_TARGET_REAL,
               PEP_TARGET_IMAGINARY} PEPWhich;

PETSC_EXTERN PetscErrorCode PEPCreate(MPI_Comm,PEP*);
PETSC_EXTERN PetscErrorCode PEPDestroy(PEP*);
PETSC_EXTERN PetscErrorCode PEPReset(PEP);
PETSC_EXTERN PetscErrorCode PEPSetType(PEP,PEPType);
PETSC_EXTERN PetscErrorCode PEPGetType(PEP,PEPType*);
PETSC_EXTERN PetscErrorCode PEPSetProblemType(PEP,PEPProblemType);
PETSC_EXTERN PetscErrorCode PEPGetProblemType(PEP,PEPProblemType*);
PETSC_EXTERN PetscErrorCode PEPSetOperators(PEP,PetscInt,Mat[]);
PETSC_EXTERN PetscErrorCode PEPGetOperators(PEP,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode PEPGetNumMatrices(PEP,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPSetTarget(PEP,PetscScalar);
PETSC_EXTERN PetscErrorCode PEPGetTarget(PEP,PetscScalar*);
PETSC_EXTERN PetscErrorCode PEPSetFromOptions(PEP);
PETSC_EXTERN PetscErrorCode PEPSetUp(PEP);
PETSC_EXTERN PetscErrorCode PEPSolve(PEP);
PETSC_EXTERN PetscErrorCode PEPView(PEP,PetscViewer);
PETSC_EXTERN PetscErrorCode PEPPrintSolution(PEP,PetscViewer);

PETSC_EXTERN PetscErrorCode PEPSetIP(PEP,IP);
PETSC_EXTERN PetscErrorCode PEPGetIP(PEP,IP*);
PETSC_EXTERN PetscErrorCode PEPSetDS(PEP,DS);
PETSC_EXTERN PetscErrorCode PEPGetDS(PEP,DS*);
PETSC_EXTERN PetscErrorCode PEPSetST(PEP,ST);
PETSC_EXTERN PetscErrorCode PEPGetST(PEP,ST*);

PETSC_EXTERN PetscErrorCode PEPSetTolerances(PEP,PetscReal,PetscInt);
PETSC_EXTERN PetscErrorCode PEPGetTolerances(PEP,PetscReal*,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPSetConvergenceTest(PEP,PetscErrorCode (*)(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*),void*);
PETSC_EXTERN PetscErrorCode PEPConvergedDefault(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
PETSC_EXTERN PetscErrorCode PEPConvergedAbsolute(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
PETSC_EXTERN PetscErrorCode PEPSetDimensions(PEP,PetscInt,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode PEPGetDimensions(PEP,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPSetScaleFactor(PEP,PetscReal);
PETSC_EXTERN PetscErrorCode PEPGetScaleFactor(PEP,PetscReal*);

PETSC_EXTERN PetscErrorCode PEPGetConverged(PEP,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPGetEigenpair(PEP,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
PETSC_EXTERN PetscErrorCode PEPComputeRelativeError(PEP,PetscInt,PetscReal*);
PETSC_EXTERN PetscErrorCode PEPComputeResidualNorm(PEP,PetscInt,PetscReal*);
PETSC_EXTERN PetscErrorCode PEPGetErrorEstimate(PEP,PetscInt,PetscReal*);

PETSC_EXTERN PetscErrorCode PEPMonitor(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt);
PETSC_EXTERN PetscErrorCode PEPMonitorSet(PEP,PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void*,PetscErrorCode (*)(void**));
PETSC_EXTERN PetscErrorCode PEPMonitorCancel(PEP);
PETSC_EXTERN PetscErrorCode PEPGetMonitorContext(PEP,void **);
PETSC_EXTERN PetscErrorCode PEPGetIterationNumber(PEP,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPGetOperationCounters(PEP,PetscInt*,PetscInt*,PetscInt*);

PETSC_EXTERN PetscErrorCode PEPSetInitialSpace(PEP,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode PEPSetInitialSpaceLeft(PEP,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode PEPSetWhichEigenpairs(PEP,PEPWhich);
PETSC_EXTERN PetscErrorCode PEPGetWhichEigenpairs(PEP,PEPWhich*);
PETSC_EXTERN PetscErrorCode PEPSetLeftVectorsWanted(PEP,PetscBool);
PETSC_EXTERN PetscErrorCode PEPGetLeftVectorsWanted(PEP,PetscBool*);
PETSC_EXTERN PetscErrorCode PEPSetEigenvalueComparison(PEP,PetscErrorCode (*func)(PEP,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*),void*);

PETSC_EXTERN PetscErrorCode PEPMonitorAll(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
PETSC_EXTERN PetscErrorCode PEPMonitorFirst(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
PETSC_EXTERN PetscErrorCode PEPMonitorConverged(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
PETSC_EXTERN PetscErrorCode PEPMonitorLG(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
PETSC_EXTERN PetscErrorCode PEPMonitorLGAll(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);

PETSC_EXTERN PetscErrorCode PEPSetTrackAll(PEP,PetscBool);
PETSC_EXTERN PetscErrorCode PEPGetTrackAll(PEP,PetscBool*);

PETSC_EXTERN PetscErrorCode PEPSetOptionsPrefix(PEP,const char*);
PETSC_EXTERN PetscErrorCode PEPAppendOptionsPrefix(PEP,const char*);
PETSC_EXTERN PetscErrorCode PEPGetOptionsPrefix(PEP,const char*[]);

/*E
    PEPConvergedReason - Reason an eigensolver was said to
         have converged or diverged

    Level: beginner

.seealso: PEPSolve(), PEPGetConvergedReason(), PEPSetTolerances()
E*/
typedef enum {/* converged */
              PEP_CONVERGED_TOL                =  2,
              /* diverged */
              PEP_DIVERGED_ITS                 = -3,
              PEP_DIVERGED_BREAKDOWN           = -4,
              PEP_CONVERGED_ITERATING          =  0} PEPConvergedReason;

PETSC_EXTERN PetscErrorCode PEPGetConvergedReason(PEP,PEPConvergedReason *);

PETSC_EXTERN PetscErrorCode PEPSortEigenvalues(PEP,PetscInt,PetscScalar*,PetscScalar*,PetscInt*);
PETSC_EXTERN PetscErrorCode PEPCompareEigenvalues(PEP,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*);

PETSC_EXTERN PetscFunctionList PEPList;
PETSC_EXTERN PetscBool         PEPRegisterAllCalled;
PETSC_EXTERN PetscErrorCode PEPRegisterAll(void);
PETSC_EXTERN PetscErrorCode PEPRegister(const char[],PetscErrorCode(*)(PEP));

PETSC_EXTERN PetscErrorCode PEPSetWorkVecs(PEP,PetscInt);

/* --------- options specific to particular eigensolvers -------- */


#endif

