/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcnep.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = NEP */

SLEPC_EXTERN PetscBool NEPRegisterAllCalled;
SLEPC_EXTERN PetscBool NEPMonitorRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode NEPRegisterAll(void);
SLEPC_EXTERN PetscErrorCode NEPMonitorRegisterAll(void);
SLEPC_EXTERN PetscLogEvent NEP_SetUp,NEP_Solve,NEP_Refine,NEP_FunctionEval,NEP_JacobianEval,NEP_Resolvent,NEP_CISS_SVD;

typedef struct _NEPOps *NEPOps;

struct _NEPOps {
  PetscErrorCode (*solve)(NEP);
  PetscErrorCode (*setup)(NEP);
  PetscErrorCode (*setfromoptions)(NEP,PetscOptionItems);
  PetscErrorCode (*publishoptions)(NEP);
  PetscErrorCode (*destroy)(NEP);
  PetscErrorCode (*reset)(NEP);
  PetscErrorCode (*view)(NEP,PetscViewer);
  PetscErrorCode (*computevectors)(NEP);
  PetscErrorCode (*setdstype)(NEP);
};

/*
     Maximum number of monitors you can run with a single NEP
*/
#define MAXNEPMONITORS 5

typedef enum { NEP_STATE_INITIAL,
               NEP_STATE_SETUP,
               NEP_STATE_SOLVED,
               NEP_STATE_EIGENVECTORS } NEPStateType;

/*
     How the problem function T(lambda) has been defined by the user
     - Callback: one callback to build the function matrix, another one for the Jacobian
     - Split: in split form sum_j(A_j*f_j(lambda))
*/
typedef enum { NEP_USER_INTERFACE_CALLBACK=1,
               NEP_USER_INTERFACE_SPLIT } NEPUserInterface;

/*
   To check for unsupported features at NEPSetUp_XXX()
*/
typedef enum { NEP_FEATURE_CALLBACK=1,      /* callback user interface */
               NEP_FEATURE_REGION=4,        /* nontrivial region for filtering */
               NEP_FEATURE_CONVERGENCE=16,  /* convergence test selected by user */
               NEP_FEATURE_STOPPING=32,     /* stopping test */
               NEP_FEATURE_TWOSIDED=64      /* two-sided variant */
             } NEPFeatureType;

/*
   Defines the NEP data structure.
*/
struct _p_NEP {
  PETSCHEADER(struct _NEPOps);
  /*------------------------- User parameters ---------------------------*/
  PetscInt       max_it;           /* maximum number of iterations */
  PetscInt       nev;              /* number of eigenvalues to compute */
  PetscInt       ncv;              /* number of basis vectors */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini;             /* number of initial vectors (negative means not copied yet) */
  PetscScalar    target;           /* target value */
  PetscReal      tol;              /* tolerance */
  NEPConv        conv;             /* convergence test */
  NEPStop        stop;             /* stopping test */
  NEPWhich       which;            /* which part of the spectrum to be sought */
  NEPProblemType problem_type;     /* which kind of problem to be solved */
  NEPRefine      refine;           /* type of refinement to be applied after solve */
  PetscInt       npart;            /* number of partitions of the communicator */
  PetscReal      rtol;             /* tolerance for refinement */
  PetscInt       rits;             /* number of iterations of the refinement method */
  NEPRefineScheme scheme;          /* scheme for solving linear systems within refinement */
  PetscBool      trackall;         /* whether all the residuals must be computed */
  PetscBool      twosided;         /* whether to compute left eigenvectors (two-sided solver) */

  /*-------------- User-provided functions and contexts -----------------*/
  NEPFunctionFn        *computefunction;
  NEPJacobianFn        *computejacobian;
  void                 *functionctx;
  void                 *jacobianctx;
  NEPConvergenceTestFn *converged;
  NEPConvergenceTestFn *convergeduser;
  PetscCtxDestroyFn    *convergeddestroy;
  NEPStoppingTestFn    *stopping;
  NEPStoppingTestFn    *stoppinguser;
  PetscCtxDestroyFn    *stoppingdestroy;
  void                 *convergedctx;
  void                 *stoppingctx;
  PetscErrorCode       (*monitor[MAXNEPMONITORS])(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscCtxDestroyFn    *monitordestroy[MAXNEPMONITORS];
  void                 *monitorcontext[MAXNEPMONITORS];
  PetscInt             numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  DS             ds;               /* direct solver object */
  BV             V;                /* set of basis vectors and computed eigenvectors */
  BV             W;                /* left basis vectors (if left eigenvectors requested) */
  RG             rg;               /* optional region for filtering */
  SlepcSC        sc;               /* sorting criterion data */
  Mat            function;         /* function matrix */
  Mat            function_pre;     /* function matrix (preconditioner) */
  Mat            jacobian;         /* Jacobian matrix */
  Mat            *A;               /* matrix coefficients of split form */
  FN             *f;               /* matrix functions of split form */
  PetscInt       nt;               /* number of terms in split form */
  MatStructure   mstr;             /* pattern of split matrices */
  Mat            *P;               /* matrix coefficients of split form (preconditioner) */
  MatStructure   mstrp;            /* pattern of split matrices (preconditioner) */
  Vec            *IS;              /* references to user-provided initial space */
  PetscScalar    *eigr,*eigi;      /* real and imaginary parts of eigenvalues */
  PetscReal      *errest;          /* error estimates */
  PetscInt       *perm;            /* permutation for eigenvalue ordering */
  PetscInt       nwork;            /* number of work vectors */
  Vec            *work;            /* work vectors */
  KSP            refineksp;        /* ksp used in refinement */
  PetscSubcomm   refinesubc;       /* context for sub-communicators */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables --------------------------*/
  NEPStateType   state;            /* initial -> setup -> solved -> eigenvectors */
  PetscInt       nconv;            /* number of converged eigenvalues */
  PetscInt       its;              /* number of iterations so far computed */
  PetscInt       n,nloc;           /* problem dimensions (global, local) */
  PetscReal      *nrma;            /* computed matrix norms */
  NEPUserInterface fui;            /* how the user has defined the nonlinear operator */
  PetscBool      useds;            /* whether the solver uses the DS object or not */
  Mat            resolvent;        /* shell matrix to be used in NEPApplyResolvent */
  NEPConvergedReason reason;
};

/*
    Macros to test valid NEP arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define NEPCheckProblem(h,arg) do {(void)(h);} while (0)
#define NEPCheckCallback(h,arg) do {(void)(h);} while (0)
#define NEPCheckSplit(h,arg) do {(void)(h);} while (0)
#define NEPCheckDerivatives(h,arg) do {(void)(h);} while (0)
#define NEPCheckSolved(h,arg) do {(void)(h);} while (0)

#else

#define NEPCheckProblem(h,arg) \
  do { \
    PetscCheck(((h)->fui),PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"The nonlinear eigenproblem has not been specified yet. Parameter #%d",arg); \
  } while (0)

#define NEPCheckCallback(h,arg) \
  do { \
    PetscCheck((h)->fui==NEP_USER_INTERFACE_CALLBACK,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"This operation requires the nonlinear eigenproblem specified with callbacks. Parameter #%d",arg); \
  } while (0)

#define NEPCheckSplit(h,arg) \
  do { \
    PetscCheck((h)->fui==NEP_USER_INTERFACE_SPLIT,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"This operation requires the nonlinear eigenproblem in split form. Parameter #%d",arg); \
  } while (0)

#define NEPCheckSolved(h,arg) \
  do { \
    PetscCheck((h)->state>=NEP_STATE_SOLVED,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call NEPSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

/* Check for unsupported features */
#define NEPCheckUnsupportedCondition(nep,mask,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!((mask) & NEP_FEATURE_CALLBACK) || (nep)->fui!=NEP_USER_INTERFACE_CALLBACK,PetscObjectComm((PetscObject)(nep)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used with callback functions (use the split operator)",((PetscObject)(nep))->type_name,(msg)); \
      if ((mask) & NEP_FEATURE_REGION) { \
        PetscBool      __istrivial; \
        PetscCall(RGIsTrivial((nep)->rg,&__istrivial)); \
        PetscCheck(__istrivial,PetscObjectComm((PetscObject)(nep)),PETSC_ERR_SUP,"The solver '%s'%s does not support region filtering",((PetscObject)(nep))->type_name,(msg)); \
      } \
      PetscCheck(!((mask) & NEP_FEATURE_CONVERGENCE) || (nep)->converged==NEPConvergedRelative,PetscObjectComm((PetscObject)(nep)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default convergence test",((PetscObject)(nep))->type_name,(msg)); \
      PetscCheck(!((mask) & NEP_FEATURE_STOPPING) || (nep)->stopping==NEPStoppingBasic,PetscObjectComm((PetscObject)(nep)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default stopping test",((PetscObject)(nep))->type_name,(msg)); \
      PetscCheck(!((mask) & NEP_FEATURE_TWOSIDED) || !(nep)->twosided,PetscObjectComm((PetscObject)(nep)),PETSC_ERR_SUP,"The solver '%s'%s cannot compute left eigenvectors (no two-sided variant)",((PetscObject)(nep))->type_name,(msg)); \
    } \
  } while (0)
#define NEPCheckUnsupported(nep,mask) NEPCheckUnsupportedCondition(nep,mask,PETSC_TRUE,"")

/* Check for ignored features */
#define NEPCheckIgnoredCondition(nep,mask,condition,msg) \
  do { \
    if (condition) { \
      if (((mask) & NEP_FEATURE_CALLBACK) && (nep)->fui==NEP_USER_INTERFACE_CALLBACK) PetscCall(PetscInfo((nep),"The solver '%s'%s ignores the user interface settings\n",((PetscObject)(nep))->type_name,(msg))); \
      if ((mask) & NEP_FEATURE_REGION) { \
        PetscBool __istrivial; \
        PetscCall(RGIsTrivial((nep)->rg,&__istrivial)); \
        if (!__istrivial) PetscCall(PetscInfo((nep),"The solver '%s'%s ignores the specified region\n",((PetscObject)(nep))->type_name,(msg))); \
      } \
      if (((mask) & NEP_FEATURE_CONVERGENCE) && (nep)->converged!=NEPConvergedRelative) PetscCall(PetscInfo((nep),"The solver '%s'%s ignores the convergence test settings\n",((PetscObject)(nep))->type_name,(msg))); \
      if (((mask) & NEP_FEATURE_STOPPING) && (nep)->stopping!=NEPStoppingBasic) PetscCall(PetscInfo((nep),"The solver '%s'%s ignores the stopping test settings\n",((PetscObject)(nep))->type_name,(msg))); \
      if (((mask) & NEP_FEATURE_TWOSIDED) && (nep)->twosided) PetscCall(PetscInfo((nep),"The solver '%s'%s ignores the two-sided flag\n",((PetscObject)(nep))->type_name,(msg))); \
    } \
  } while (0)
#define NEPCheckIgnored(nep,mask) NEPCheckIgnoredCondition(nep,mask,PETSC_TRUE,"")

/*
  NEP_KSPSetOperators - Sets the KSP matrices
*/
static inline PetscErrorCode NEP_KSPSetOperators(KSP ksp,Mat A,Mat B)
{
  const char     *prefix;

  PetscFunctionBegin;
  PetscCall(KSPSetOperators(ksp,A,B));
  PetscCall(MatGetOptionsPrefix(B,&prefix));
  if (!prefix) {
    /* set Mat prefix to be the same as KSP to enable setting command-line options (e.g. MUMPS)
       only applies if the Mat has no user-defined prefix */
    PetscCall(KSPGetOptionsPrefix(ksp,&prefix));
    PetscCall(MatSetOptionsPrefix(B,prefix));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode NEPSetDimensions_Default(NEP,PetscInt,PetscInt*,PetscInt*);
SLEPC_INTERN PetscErrorCode NEPComputeVectors(NEP);
SLEPC_INTERN PetscErrorCode NEPReset_Problem(NEP);
SLEPC_INTERN PetscErrorCode NEPGetDefaultShift(NEP,PetscScalar*);
SLEPC_INTERN PetscErrorCode NEPComputeVectors_Schur(NEP);
SLEPC_INTERN PetscErrorCode NEPComputeResidualNorm_Private(NEP,PetscBool,PetscScalar,Vec,Vec*,PetscReal*);
SLEPC_INTERN PetscErrorCode NEPNewtonRefinementSimple(NEP,PetscInt*,PetscReal,PetscInt);
