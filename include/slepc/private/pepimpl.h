/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcpep.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = PEP */

SLEPC_EXTERN PetscBool PEPRegisterAllCalled;
SLEPC_EXTERN PetscBool PEPMonitorRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode PEPRegisterAll(void);
SLEPC_EXTERN PetscErrorCode PEPMonitorRegisterAll(void);
SLEPC_EXTERN PetscLogEvent PEP_SetUp,PEP_Solve,PEP_Refine,PEP_CISS_SVD;

typedef struct _PEPOps *PEPOps;

struct _PEPOps {
  PetscErrorCode (*solve)(PEP);
  PetscErrorCode (*setup)(PEP);
  PetscErrorCode (*setfromoptions)(PEP,PetscOptionItems);
  PetscErrorCode (*publishoptions)(PEP);
  PetscErrorCode (*destroy)(PEP);
  PetscErrorCode (*reset)(PEP);
  PetscErrorCode (*view)(PEP,PetscViewer);
  PetscErrorCode (*backtransform)(PEP);
  PetscErrorCode (*computevectors)(PEP);
  PetscErrorCode (*extractvectors)(PEP);
  PetscErrorCode (*setdefaultst)(PEP);
  PetscErrorCode (*setdstype)(PEP);
};

/*
     Maximum number of monitors you can run with a single PEP
*/
#define MAXPEPMONITORS 5

typedef enum { PEP_STATE_INITIAL,
               PEP_STATE_SETUP,
               PEP_STATE_SOLVED,
               PEP_STATE_EIGENVECTORS } PEPStateType;

/*
   To check for unsupported features at PEPSetUp_XXX()
*/
typedef enum { PEP_FEATURE_NONMONOMIAL=1,   /* non-monomial bases */
               PEP_FEATURE_REGION=4,        /* nontrivial region for filtering */
               PEP_FEATURE_EXTRACT=8,       /* eigenvector extraction */
               PEP_FEATURE_CONVERGENCE=16,  /* convergence test selected by user */
               PEP_FEATURE_STOPPING=32,     /* stopping test */
               PEP_FEATURE_SCALE=64         /* scaling */
             } PEPFeatureType;

/*
   Defines the PEP data structure.
*/
struct _p_PEP {
  PETSCHEADER(struct _PEPOps);
  /*------------------------- User parameters ---------------------------*/
  PetscInt       max_it;           /* maximum number of iterations */
  PetscInt       nev;              /* number of eigenvalues to compute */
  PetscInt       ncv;              /* number of basis vectors */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini;             /* number of initial vectors (negative means not copied yet) */
  PetscScalar    target;           /* target value */
  PetscReal      tol;              /* tolerance */
  PEPConv        conv;             /* convergence test */
  PEPStop        stop;             /* stopping test */
  PEPWhich       which;            /* which part of the spectrum to be sought */
  PetscReal      inta,intb;        /* interval [a,b] for spectrum slicing */
  PEPBasis       basis;            /* polynomial basis used to represent the problem */
  PEPProblemType problem_type;     /* which kind of problem to be solved */
  PEPScale       scale;            /* scaling strategy to be used */
  PetscReal      sfactor,dsfactor; /* scaling factors */
  PetscInt       sits;             /* number of iterations of the scaling method */
  PetscReal      slambda;          /* norm eigenvalue approximation for scaling */
  PEPRefine      refine;           /* type of refinement to be applied after solve */
  PetscInt       npart;            /* number of partitions of the communicator */
  PetscReal      rtol;             /* tolerance for refinement */
  PetscInt       rits;             /* number of iterations of the refinement method */
  PEPRefineScheme scheme;          /* scheme for solving linear systems within refinement */
  PEPExtract     extract;          /* type of extraction used */
  PetscBool      trackall;         /* whether all the residuals must be computed */

  /*-------------- User-provided functions and contexts -----------------*/
  PEPConvergenceTestFn *converged;
  PEPConvergenceTestFn *convergeduser;
  PetscCtxDestroyFn    *convergeddestroy;
  PEPStoppingTestFn    *stopping;
  PEPStoppingTestFn    *stoppinguser;
  PetscCtxDestroyFn    *stoppingdestroy;
  void                 *convergedctx;
  void                 *stoppingctx;
  PetscErrorCode       (*monitor[MAXPEPMONITORS])(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscCtxDestroyFn    *monitordestroy[MAXPEPMONITORS];
  void                 *monitorcontext[MAXPEPMONITORS];
  PetscInt             numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  ST             st;               /* spectral transformation object */
  DS             ds;               /* direct solver object */
  BV             V;                /* set of basis vectors and computed eigenvectors */
  RG             rg;               /* optional region for filtering */
  SlepcSC        sc;               /* sorting criterion data */
  Mat            *A;               /* coefficient matrices of the polynomial */
  PetscInt       nmat;             /* number of matrices */
  Vec            Dl,Dr;            /* diagonal matrices for balancing */
  Vec            *IS;              /* references to user-provided initial space */
  PetscScalar    *eigr,*eigi;      /* real and imaginary parts of eigenvalues */
  PetscReal      *errest;          /* error estimates */
  PetscInt       *perm;            /* permutation for eigenvalue ordering */
  PetscReal      *pbc;             /* coefficients defining the polynomial basis */
  PetscScalar    *solvematcoeffs;  /* coefficients to compute the matrix to be inverted */
  PetscInt       nwork;            /* number of work vectors */
  Vec            *work;            /* work vectors */
  KSP            refineksp;        /* ksp used in refinement */
  PetscSubcomm   refinesubc;       /* context for sub-communicators */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables --------------------------*/
  PEPStateType   state;            /* initial -> setup -> solved -> eigenvectors */
  PetscInt       nconv;            /* number of converged eigenvalues */
  PetscInt       its;              /* number of iterations so far computed */
  PetscInt       n,nloc;           /* problem dimensions (global, local) */
  PetscReal      *nrma;            /* computed matrix norms */
  PetscReal      nrml[2];          /* computed matrix norms for the linearization */
  PetscBool      sfactor_set;      /* flag to indicate the user gave sfactor */
  PetscBool      lineariz;         /* current solver is based on linearization */
  PEPConvergedReason reason;
};

/*
    Macros to test valid PEP arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define PEPCheckSolved(h,arg) do {(void)(h);} while (0)

#else

#define PEPCheckSolved(h,arg) \
  do { \
    PetscCheck((h)->state>=PEP_STATE_SOLVED,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

/*
    Macros to check settings at PEPSetUp()
*/

/* PEPCheckHermitian: the problem is Hermitian or Hyperbolic */
#define PEPCheckHermitianCondition(pep,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck((pep)->problem_type==PEP_HERMITIAN || (pep)->problem_type==PEP_HYPERBOLIC,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s can only be used for Hermitian (or hyperbolic) problems",((PetscObject)(pep))->type_name,(msg)); \
    } \
  } while (0)
#define PEPCheckHermitian(pep) PEPCheckHermitianCondition(pep,PETSC_TRUE,"")

/* PEPCheckQuadratic: the polynomial has degree 2 */
#define PEPCheckQuadraticCondition(pep,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck((pep)->nmat==3,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s is only available for quadratic problems",((PetscObject)(pep))->type_name,(msg)); \
    } \
  } while (0)
#define PEPCheckQuadratic(pep) PEPCheckQuadraticCondition(pep,PETSC_TRUE,"")

/* PEPCheckShiftSinvert: shift or shift-and-invert ST */
#define PEPCheckShiftSinvertCondition(pep,condition,msg) \
  do { \
    if (condition) { \
      PetscBool __flg; \
      PetscCall(PetscObjectTypeCompareAny((PetscObject)(pep)->st,&__flg,STSINVERT,STSHIFT,"")); \
      PetscCheck(__flg,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s requires shift or shift-and-invert spectral transform",((PetscObject)(pep))->type_name,(msg)); \
    } \
  } while (0)
#define PEPCheckShiftSinvert(pep) PEPCheckShiftSinvertCondition(pep,PETSC_TRUE,"")

/* PEPCheckSinvertCayley: shift-and-invert or Cayley ST */
#define PEPCheckSinvertCayleyCondition(pep,condition,msg) \
  do { \
    if (condition) { \
      PetscBool __flg; \
      PetscCall(PetscObjectTypeCompareAny((PetscObject)(pep)->st,&__flg,STSINVERT,STCAYLEY,"")); \
      PetscCheck(__flg,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s requires shift-and-invert or Cayley transform",((PetscObject)(pep))->type_name,(msg)); \
    } \
  } while (0)
#define PEPCheckSinvertCayley(pep) PEPCheckSinvertCayleyCondition(pep,PETSC_TRUE,"")

/* Check for unsupported features */
#define PEPCheckUnsupportedCondition(pep,mask,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!((mask) & PEP_FEATURE_NONMONOMIAL) || (pep)->basis==PEP_BASIS_MONOMIAL,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s is not implemented for non-monomial bases",((PetscObject)(pep))->type_name,(msg)); \
      if ((mask) & PEP_FEATURE_REGION) { \
        PetscBool      __istrivial; \
        PetscCall(RGIsTrivial((pep)->rg,&__istrivial)); \
        PetscCheck(__istrivial,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s does not support region filtering",((PetscObject)(pep))->type_name,(msg)); \
      } \
      PetscCheck(!((mask) & PEP_FEATURE_EXTRACT) || !(pep)->extract || (pep)->extract==PEP_EXTRACT_NONE,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s does not support extraction variants",((PetscObject)(pep))->type_name,(msg)); \
      PetscCheck(!((mask) & PEP_FEATURE_CONVERGENCE) || (pep)->converged==PEPConvergedRelative,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default convergence test",((PetscObject)(pep))->type_name,(msg)); \
      PetscCheck(!((mask) & PEP_FEATURE_STOPPING) || (pep)->stopping==PEPStoppingBasic,PetscObjectComm((PetscObject)(pep)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default stopping test",((PetscObject)(pep))->type_name,(msg)); \
    } \
  } while (0)
#define PEPCheckUnsupported(pep,mask) PEPCheckUnsupportedCondition(pep,mask,PETSC_TRUE,"")

/* Check for ignored features */
#define PEPCheckIgnoredCondition(pep,mask,condition,msg) \
  do { \
    if (condition) { \
      if (((mask) & PEP_FEATURE_NONMONOMIAL) && (pep)->basis!=PEP_BASIS_MONOMIAL) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the basis settings\n",((PetscObject)(pep))->type_name,(msg))); \
      if ((mask) & PEP_FEATURE_REGION) { \
        PetscBool __istrivial; \
        PetscCall(RGIsTrivial((pep)->rg,&__istrivial)); \
        if (!__istrivial) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the specified region\n",((PetscObject)(pep))->type_name,(msg))); \
      } \
      if (((mask) & PEP_FEATURE_EXTRACT) && (pep)->extract && (pep)->extract!=PEP_EXTRACT_NONE) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the extract settings\n",((PetscObject)(pep))->type_name,(msg))); \
      if (((mask) & PEP_FEATURE_CONVERGENCE) && (pep)->converged!=PEPConvergedRelative) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the convergence test settings\n",((PetscObject)(pep))->type_name,(msg))); \
      if (((mask) & PEP_FEATURE_STOPPING) && (pep)->stopping!=PEPStoppingBasic) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the stopping test settings\n",((PetscObject)(pep))->type_name,(msg))); \
      if (((mask) & PEP_FEATURE_SCALE) && (pep)->scale!=PEP_SCALE_NONE) PetscCall(PetscInfo((pep),"The solver '%s'%s ignores the scaling settings\n",((PetscObject)(pep))->type_name,(msg))); \
    } \
  } while (0)
#define PEPCheckIgnored(pep,mask) PEPCheckIgnoredCondition(pep,mask,PETSC_TRUE,"")

/*
  PEP_KSPSetOperators - Sets the KSP matrices
*/
static inline PetscErrorCode PEP_KSPSetOperators(KSP ksp,Mat A,Mat B)
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

SLEPC_INTERN PetscErrorCode PEPSetWhichEigenpairs_Default(PEP);
SLEPC_INTERN PetscErrorCode PEPSetDimensions_Default(PEP,PetscInt,PetscInt*,PetscInt*);
SLEPC_INTERN PetscErrorCode PEPExtractVectors(PEP);
SLEPC_INTERN PetscErrorCode PEPBackTransform_Default(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors_Default(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors_Indefinite(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeResidualNorm_Private(PEP,PetscScalar,PetscScalar,Vec,Vec,Vec*,PetscReal*);
SLEPC_INTERN PetscErrorCode PEPKrylovConvergence(PEP,PetscBool,PetscInt,PetscInt,PetscReal,PetscInt*);
SLEPC_INTERN PetscErrorCode PEPComputeScaleFactor(PEP);
SLEPC_INTERN PetscErrorCode PEPBuildDiagonalScaling(PEP);
SLEPC_INTERN PetscErrorCode PEPBasisCoefficients(PEP,PetscReal*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasis(PEP,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasisDerivative(PEP,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasisMat(PEP,PetscInt,PetscScalar*,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode PEPNewtonRefinement_TOAR(PEP,PetscScalar,PetscInt*,PetscReal*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode PEPNewtonRefinementSimple(PEP,PetscInt*,PetscReal,PetscInt);
SLEPC_INTERN PetscErrorCode PEPSetDefaultST(PEP);
SLEPC_INTERN PetscErrorCode PEPSetDefaultST_Transform(PEP);
