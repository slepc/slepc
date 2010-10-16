/*
   User interface for the SLEPC singular value solvers. 

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2010, Universidad Politecnica de Valencia, Spain

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

#if !defined(__SLEPCSVD_H)
#define __SLEPCSVD_H
#include "slepcsys.h"
#include "slepceps.h"
PETSC_EXTERN_CXX_BEGIN

extern PetscClassId SVD_CLASSID;

/*S
     SVD - Abstract SLEPc object that manages all the singular value 
     problem solvers.

   Level: beginner

.seealso:  SVDCreate()
S*/
typedef struct _p_SVD* SVD;

/*E
    SVDType - String with the name of a SLEPc singular value solver

   Level: beginner

.seealso: SVDSetType(), SVD
E*/
#define SVDType        char*
#define SVDCROSS       "cross"
#define SVDCYCLIC      "cyclic"
#define SVDLAPACK      "lapack"
#define SVDLANCZOS     "lanczos"
#define SVDTRLANCZOS   "trlanczos"

/*E
    SVDTransposeMode - determines how to handle the transpose of the matrix

    Level: advanced

.seealso: SVDSetTransposeMode(), SVDGetTransposeMode()
E*/
typedef enum { SVD_TRANSPOSE_EXPLICIT,
               SVD_TRANSPOSE_IMPLICIT } SVDTransposeMode;

/*E
    SVDWhich - determines whether largest or smallest singular triplets
    are to be computed

    Level: intermediate

.seealso: SVDSetWhichSingularTriplets(), SVDGetWhichSingularTriplets()
E*/
typedef enum { SVD_LARGEST,
               SVD_SMALLEST } SVDWhich;

/*E
    SVDConvergedReason - reason a singular value solver was said to 
         have converged or diverged

   Level: beginner

.seealso: SVDSolve(), SVDGetConvergedReason(), SVDSetTolerances()
E*/
typedef enum {/* converged */
              SVD_CONVERGED_TOL                =  2,
              /* diverged */
              SVD_DIVERGED_ITS                 = -3,
              SVD_DIVERGED_BREAKDOWN           = -4,
              SVD_CONVERGED_ITERATING          =  0 } SVDConvergedReason;

EXTERN PetscErrorCode SVDCreate(MPI_Comm,SVD*);
EXTERN PetscErrorCode SVDSetIP(SVD,IP);
EXTERN PetscErrorCode SVDGetIP(SVD,IP*);
EXTERN PetscErrorCode SVDSetType(SVD,const SVDType);
EXTERN PetscErrorCode SVDGetType(SVD,const SVDType*);
EXTERN PetscErrorCode SVDSetOperator(SVD,Mat);
EXTERN PetscErrorCode SVDGetOperator(SVD,Mat*);
EXTERN PetscErrorCode SVDSetInitialSpace(SVD,PetscInt,Vec*);
EXTERN PetscErrorCode SVDSetTransposeMode(SVD,SVDTransposeMode);
EXTERN PetscErrorCode SVDGetTransposeMode(SVD,SVDTransposeMode*);
EXTERN PetscErrorCode SVDSetDimensions(SVD,PetscInt,PetscInt,PetscInt);
EXTERN PetscErrorCode SVDGetDimensions(SVD,PetscInt*,PetscInt*,PetscInt*);
EXTERN PetscErrorCode SVDSetTolerances(SVD,PetscReal,PetscInt);
EXTERN PetscErrorCode SVDGetTolerances(SVD,PetscReal*,PetscInt*);
EXTERN PetscErrorCode SVDSetWhichSingularTriplets(SVD,SVDWhich);
EXTERN PetscErrorCode SVDGetWhichSingularTriplets(SVD,SVDWhich*);
EXTERN PetscErrorCode SVDSetFromOptions(SVD);
EXTERN PetscErrorCode SVDSetOptionsPrefix(SVD,const char*);
EXTERN PetscErrorCode SVDAppendOptionsPrefix(SVD,const char*);
EXTERN PetscErrorCode SVDGetOptionsPrefix(SVD,const char*[]);
EXTERN PetscErrorCode SVDSetUp(SVD);
EXTERN PetscErrorCode SVDSolve(SVD);
EXTERN PetscErrorCode SVDGetIterationNumber(SVD,PetscInt*);
EXTERN PetscErrorCode SVDGetConvergedReason(SVD,SVDConvergedReason*);
EXTERN PetscErrorCode SVDGetConverged(SVD,PetscInt*);
EXTERN PetscErrorCode SVDGetSingularTriplet(SVD,PetscInt,PetscReal*,Vec,Vec);
EXTERN PetscErrorCode SVDComputeResidualNorms(SVD,PetscInt,PetscReal*,PetscReal*);
EXTERN PetscErrorCode SVDComputeRelativeError(SVD,PetscInt,PetscReal*);
EXTERN PetscErrorCode SVDGetOperationCounters(SVD,PetscInt*,PetscInt*);
EXTERN PetscErrorCode SVDView(SVD,PetscViewer);
EXTERN PetscErrorCode SVDDestroy(SVD);

EXTERN PetscErrorCode SVDMonitorSet(SVD,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*),
                                    void*,PetscErrorCode (*monitordestroy)(void*));
EXTERN PetscErrorCode SVDMonitorCancel(SVD);
EXTERN PetscErrorCode SVDGetMonitorContext(SVD,void **);
EXTERN PetscErrorCode SVDMonitorAll(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
EXTERN PetscErrorCode SVDMonitorFirst(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
EXTERN PetscErrorCode SVDMonitorConverged(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
EXTERN PetscErrorCode SVDMonitorLG(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
EXTERN PetscErrorCode SVDMonitorLGAll(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);

EXTERN PetscErrorCode SVDSetTrackAll(SVD,PetscBool);
EXTERN PetscErrorCode SVDGetTrackAll(SVD,PetscBool*);

EXTERN PetscErrorCode SVDDense(PetscInt,PetscInt,PetscScalar*,PetscReal*,PetscScalar*,PetscScalar*);

EXTERN PetscErrorCode SVDCrossSetEPS(SVD,EPS);
EXTERN PetscErrorCode SVDCrossGetEPS(SVD,EPS*);

EXTERN PetscErrorCode SVDCyclicSetExplicitMatrix(SVD,PetscBool);
EXTERN PetscErrorCode SVDCyclicGetExplicitMatrix(SVD,PetscBool*);
EXTERN PetscErrorCode SVDCyclicSetEPS(SVD,EPS);
EXTERN PetscErrorCode SVDCyclicGetEPS(SVD,EPS*);

EXTERN PetscErrorCode SVDLanczosSetOneSide(SVD,PetscBool);
EXTERN PetscErrorCode SVDLanczosGetOneSide(SVD,PetscBool*);

EXTERN PetscErrorCode SVDTRLanczosSetOneSide(SVD,PetscBool);
EXTERN PetscErrorCode SVDTRLanczosGetOneSide(SVD,PetscBool*);

EXTERN PetscErrorCode SVDRegister(const char*,const char*,const char*,PetscErrorCode(*)(SVD));
#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
#define SVDRegisterDynamic(a,b,c,d) SVDRegister(a,b,c,0)
#else
#define SVDRegisterDynamic(a,b,c,d) SVDRegister(a,b,c,d)
#endif
EXTERN PetscErrorCode SVDRegisterDestroy(void);

PETSC_EXTERN_CXX_END
#endif
