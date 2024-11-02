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
    ST - Abstract SLEPc object that manages spectral transformations.
    This object is accessed only in advanced applications.

    Level: beginner

.seealso:  STCreate(), EPS
S*/
typedef struct _p_ST* ST;

/*J
    STType - String with the name of a SLEPc spectral transformation

    Level: beginner

.seealso: STSetType(), ST
J*/
typedef const char *STType;
#define STSHELL     "shell"
#define STSHIFT     "shift"
#define STSINVERT   "sinvert"
#define STCAYLEY    "cayley"
#define STPRECOND   "precond"
#define STFILTER    "filter"

/* Logging support */
SLEPC_EXTERN PetscClassId ST_CLASSID;

SLEPC_EXTERN PetscErrorCode STCreate(MPI_Comm,ST*);
SLEPC_EXTERN PetscErrorCode STDestroy(ST*);
SLEPC_EXTERN PetscErrorCode STReset(ST);
SLEPC_EXTERN PetscErrorCode STSetType(ST,STType);
SLEPC_EXTERN PetscErrorCode STGetType(ST,STType*);
SLEPC_EXTERN PetscErrorCode STSetMatrices(ST,PetscInt,Mat*);
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
SLEPC_EXTERN PetscErrorCode STMatSetUp(ST,PetscScalar,PetscScalar*);
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

SLEPC_EXTERN PetscErrorCode STSetOptionsPrefix(ST,const char*);
SLEPC_EXTERN PetscErrorCode STAppendOptionsPrefix(ST,const char*);
SLEPC_EXTERN PetscErrorCode STGetOptionsPrefix(ST,const char*[]);

SLEPC_EXTERN PetscErrorCode STBackTransform(ST,PetscInt,PetscScalar*,PetscScalar*);
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
    STMatMode - Determines how to handle the coefficient matrix associated
    to the spectral transformation

    Level: intermediate

.seealso: STSetMatMode(), STGetMatMode()
E*/
typedef enum { ST_MATMODE_COPY,
               ST_MATMODE_INPLACE,
               ST_MATMODE_SHELL } STMatMode;
SLEPC_EXTERN const char *STMatModes[];

SLEPC_EXTERN PetscErrorCode STSetMatMode(ST,STMatMode);
SLEPC_EXTERN PetscErrorCode STGetMatMode(ST,STMatMode*);
SLEPC_EXTERN PetscErrorCode STSetMatStructure(ST,MatStructure);
SLEPC_EXTERN PetscErrorCode STGetMatStructure(ST,MatStructure*);

SLEPC_EXTERN PetscFunctionList STList;
SLEPC_EXTERN PetscErrorCode STRegister(const char[],PetscErrorCode(*)(ST));

/* --------- options specific to particular spectral transformations-------- */

SLEPC_EXTERN PetscErrorCode STShellGetContext(ST,void*);
SLEPC_EXTERN PetscErrorCode STShellSetContext(ST,void*);
SLEPC_EXTERN PetscErrorCode STShellSetApply(ST,PetscErrorCode (*)(ST,Vec,Vec));
SLEPC_EXTERN PetscErrorCode STShellSetApplyTranspose(ST,PetscErrorCode (*)(ST,Vec,Vec));
SLEPC_EXTERN PetscErrorCode STShellSetApplyHermitianTranspose(ST,PetscErrorCode (*)(ST,Vec,Vec));
SLEPC_EXTERN PetscErrorCode STShellSetBackTransform(ST,PetscErrorCode (*)(ST,PetscInt,PetscScalar*,PetscScalar*));

SLEPC_EXTERN PetscErrorCode STCayleyGetAntishift(ST,PetscScalar*);
SLEPC_EXTERN PetscErrorCode STCayleySetAntishift(ST,PetscScalar);

PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STGetPreconditionerMat()", ) static inline PetscErrorCode STPrecondGetMatForPC(ST st,Mat *A) {return STGetPreconditionerMat(st,A);}
PETSC_DEPRECATED_FUNCTION(3, 15, 0, "STSetPreconditionerMat()", ) static inline PetscErrorCode STPrecondSetMatForPC(ST st,Mat A) {return STSetPreconditionerMat(st,A);}
SLEPC_EXTERN PetscErrorCode STPrecondGetKSPHasMat(ST,PetscBool*);
SLEPC_EXTERN PetscErrorCode STPrecondSetKSPHasMat(ST,PetscBool);

SLEPC_EXTERN PetscErrorCode STFilterSetInterval(ST,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode STFilterGetInterval(ST,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode STFilterSetRange(ST,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode STFilterGetRange(ST,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode STFilterSetDegree(ST,PetscInt);
SLEPC_EXTERN PetscErrorCode STFilterGetDegree(ST,PetscInt*);
SLEPC_EXTERN PetscErrorCode STFilterGetThreshold(ST,PetscReal*);
