/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepcsvd.h>
#include <slepc/private/slepcimpl.h>

/* SUBMANSEC = SVD */

SLEPC_EXTERN PetscBool SVDRegisterAllCalled;
SLEPC_EXTERN PetscBool SVDMonitorRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode SVDRegisterAll(void);
SLEPC_EXTERN PetscErrorCode SVDMonitorRegisterAll(void);
SLEPC_EXTERN PetscLogEvent SVD_SetUp,SVD_Solve;

typedef struct _SVDOps *SVDOps;

struct _SVDOps {
  PetscErrorCode (*solve)(SVD);
  PetscErrorCode (*solveg)(SVD);
  PetscErrorCode (*solveh)(SVD);
  PetscErrorCode (*setup)(SVD);
  PetscErrorCode (*setfromoptions)(SVD,PetscOptionItems*);
  PetscErrorCode (*publishoptions)(SVD);
  PetscErrorCode (*destroy)(SVD);
  PetscErrorCode (*reset)(SVD);
  PetscErrorCode (*view)(SVD,PetscViewer);
  PetscErrorCode (*computevectors)(SVD);
  PetscErrorCode (*setdstype)(SVD);
};

/*
     Maximum number of monitors you can run with a single SVD
*/
#define MAXSVDMONITORS 5

typedef enum { SVD_STATE_INITIAL,
               SVD_STATE_SETUP,
               SVD_STATE_SOLVED,
               SVD_STATE_VECTORS } SVDStateType;

/*
   To check for unsupported features at SVDSetUp_XXX()
*/
typedef enum { SVD_FEATURE_CONVERGENCE=16,  /* convergence test selected by user */
               SVD_FEATURE_STOPPING=32,     /* stopping test */
               SVD_FEATURE_THRESHOLD=64     /* threshold stopping test */
             } SVDFeatureType;

/*
   Defines the SVD data structure.
*/
struct _p_SVD {
  PETSCHEADER(struct _SVDOps);
  /*------------------------- User parameters ---------------------------*/
  Mat            OP,OPb;           /* problem matrices */
  Vec            omega;            /* signature for hyperbolic problems */
  PetscInt       max_it;           /* max iterations */
  PetscInt       nsv;              /* number of requested values */
  PetscInt       ncv;              /* basis size */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini,ninil;       /* number of initial vecs (negative means not copied yet) */
  PetscReal      tol;              /* tolerance */
  PetscReal      thres;            /* threshold value */
  PetscBool      threlative;       /* threshold is relative */
  SVDConv        conv;             /* convergence test */
  SVDStop        stop;             /* stopping test */
  SVDWhich       which;            /* which singular values are computed */
  SVDProblemType problem_type;     /* which kind of problem to be solved */
  PetscBool      impltrans;        /* implicit transpose mode */
  PetscBool      trackall;         /* whether all the residuals must be computed */

  /*-------------- User-provided functions and contexts -----------------*/
  SVDConvergenceTestFn *converged;
  SVDConvergenceTestFn *convergeduser;
  PetscCtxDestroyFn    *convergeddestroy;
  SVDStoppingTestFn    *stopping;
  SVDStoppingTestFn    *stoppinguser;
  PetscCtxDestroyFn    *stoppingdestroy;
  void                 *convergedctx;
  void                 *stoppingctx;
  PetscErrorCode       (*monitor[MAXSVDMONITORS])(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
  PetscCtxDestroyFn    *monitordestroy[MAXSVDMONITORS];
  void                 *monitorcontext[MAXSVDMONITORS];
  PetscInt             numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  DS             ds;               /* direct solver object */
  BV             U,V;              /* left and right singular vectors */
  SlepcSC        sc;               /* sorting criterion data */
  Mat            A,B;              /* problem matrices */
  Mat            AT,BT;            /* transposed matrices */
  Vec            *IS,*ISL;         /* placeholder for references to user initial space */
  PetscReal      *sigma;           /* singular values */
  PetscReal      *errest;          /* error estimates */
  PetscReal      *sign;            /* +-1 for each singular value in hyperbolic problems=U'*Omega*U */
  PetscInt       *perm;            /* permutation for singular value ordering */
  PetscInt       nworkl,nworkr;    /* number of work vectors */
  Vec            *workl,*workr;    /* work vectors */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables -------------------------- */
  SVDStateType   state;            /* initial -> setup -> solved -> vectors */
  PetscInt       nconv;            /* number of converged values */
  PetscInt       its;              /* iteration counter */
  PetscBool      leftbasis;        /* if U is filled by the solver */
  PetscBool      swapped;          /* the U and V bases have been swapped (M<N) */
  PetscBool      expltrans;        /* explicit transpose created */
  PetscReal      nrma,nrmb;        /* computed matrix norms */
  PetscBool      isgeneralized;
  PetscBool      ishyperbolic;
  SVDConvergedReason reason;
};

/*
    Macros to test valid SVD arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define SVDCheckSolved(h,arg) do {(void)(h);} while (0)

#else

#define SVDCheckSolved(h,arg) \
  do { \
    PetscCheck((h)->state>=SVD_STATE_SOLVED,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call SVDSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

/*
    Macros to check settings at SVDSetUp()
*/

/* SVDCheckStandard: the problem is not GSVD */
#define SVDCheckStandardCondition(svd,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!(svd)->isgeneralized,PetscObjectComm((PetscObject)(svd)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for generalized problems",((PetscObject)(svd))->type_name,(msg)); \
    } \
  } while (0)
#define SVDCheckStandard(svd) SVDCheckStandardCondition(svd,PETSC_TRUE,"")

/* SVDCheckDefinite: the problem is not hyperbolic */
#define SVDCheckDefiniteCondition(svd,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!(svd)->ishyperbolic,PetscObjectComm((PetscObject)(svd)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for hyperbolic problems",((PetscObject)(svd))->type_name,(msg)); \
    } \
  } while (0)
#define SVDCheckDefinite(svd) SVDCheckDefiniteCondition(svd,PETSC_TRUE,"")

/* Check for unsupported features */
#define SVDCheckUnsupportedCondition(svd,mask,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!((mask) & SVD_FEATURE_CONVERGENCE) || (svd)->converged==SVDConvergedRelative,PetscObjectComm((PetscObject)(svd)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default convergence test",((PetscObject)(svd))->type_name,(msg)); \
      PetscCheck(!((mask) & SVD_FEATURE_STOPPING) || (svd)->stopping==SVDStoppingBasic,PetscObjectComm((PetscObject)(svd)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default stopping test",((PetscObject)(svd))->type_name,(msg)); \
      PetscCheck(!((mask) & SVD_FEATURE_THRESHOLD) || (svd)->stopping!=SVDStoppingThreshold,PetscObjectComm((PetscObject)(svd)),PETSC_ERR_SUP,"The solver '%s'%s does not support the threshold stopping test",((PetscObject)(svd))->type_name,(msg)); \
    } \
  } while (0)
#define SVDCheckUnsupported(svd,mask) SVDCheckUnsupportedCondition(svd,mask,PETSC_TRUE,"")

/* Check for ignored features */
#define SVDCheckIgnoredCondition(svd,mask,condition,msg) \
  do { \
    if (condition) { \
      if (((mask) & SVD_FEATURE_CONVERGENCE) && (svd)->converged!=SVDConvergedRelative) PetscCall(PetscInfo((svd),"The solver '%s'%s ignores the convergence test settings\n",((PetscObject)(svd))->type_name,(msg))); \
      if (((mask) & SVD_FEATURE_STOPPING) && (svd)->stopping!=SVDStoppingBasic) PetscCall(PetscInfo((svd),"The solver '%s'%s ignores the stopping test settings\n",((PetscObject)(svd))->type_name,(msg))); \
    } \
  } while (0)
#define SVDCheckIgnored(svd,mask) SVDCheckIgnoredCondition(svd,mask,PETSC_TRUE,"")

/*
    SVDSetCtxThreshold - Fills SVDStoppingCtx with data needed for the threshold stopping test
*/
#define SVDSetCtxThreshold(svd,sigma,k) \
  do { \
    if (svd->stop==SVD_STOP_THRESHOLD && k) { \
      ((SVDStoppingCtx)svd->stoppingctx)->firstsv = sigma[0]; \
      ((SVDStoppingCtx)svd->stoppingctx)->lastsv  = sigma[k-1]; \
    } \
  } while (0)

/*
  SVD_KSPSetOperators - Sets the KSP matrices
*/
static inline PetscErrorCode SVD_KSPSetOperators(KSP ksp,Mat A,Mat B)
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

/*
   Create the template vector for the left basis in GSVD, as in
   MatCreateVecsEmpty(Z,NULL,&t) for Z=[A;B] without forming Z.
 */
static inline PetscErrorCode SVDCreateLeftTemplate(SVD svd,Vec *t)
{
  PetscInt       M,P,m,p;
  Vec            v1,v2;
  VecType        vec_type;

  PetscFunctionBegin;
  PetscCall(MatCreateVecsEmpty(svd->OP,NULL,&v1));
  PetscCall(VecGetSize(v1,&M));
  PetscCall(VecGetLocalSize(v1,&m));
  PetscCall(VecGetType(v1,&vec_type));
  PetscCall(MatCreateVecsEmpty(svd->OPb,NULL,&v2));
  PetscCall(VecGetSize(v2,&P));
  PetscCall(VecGetLocalSize(v2,&p));
  PetscCall(VecCreate(PetscObjectComm((PetscObject)(v1)),t));
  PetscCall(VecSetType(*t,vec_type));
  PetscCall(VecSetSizes(*t,m+p,M+P));
  PetscCall(VecSetUp(*t));
  PetscCall(VecDestroy(&v1));
  PetscCall(VecDestroy(&v2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode SVDKrylovConvergence(SVD,PetscBool,PetscInt,PetscInt,PetscReal,PetscInt*);
SLEPC_INTERN PetscErrorCode SVDTwoSideLanczos(SVD,PetscReal*,PetscReal*,BV,BV,PetscInt,PetscInt*,PetscBool*);
SLEPC_INTERN PetscErrorCode SVDSetDimensions_Default(SVD);
SLEPC_INTERN PetscErrorCode SVDComputeVectors(SVD);
SLEPC_INTERN PetscErrorCode SVDComputeVectors_Left(SVD);
