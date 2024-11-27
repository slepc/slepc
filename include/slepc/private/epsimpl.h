/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#pragma once

#include <slepceps.h>
#include <slepc/private/bvimpl.h>

/* SUBMANSEC = EPS */

SLEPC_EXTERN PetscBool EPSRegisterAllCalled;
SLEPC_EXTERN PetscBool EPSMonitorRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode EPSRegisterAll(void);
SLEPC_EXTERN PetscErrorCode EPSMonitorRegisterAll(void);
SLEPC_EXTERN PetscLogEvent EPS_SetUp,EPS_Solve,EPS_CISS_SVD;

typedef struct _EPSOps *EPSOps;

struct _EPSOps {
  PetscErrorCode (*solve)(EPS);
  PetscErrorCode (*setup)(EPS);
  PetscErrorCode (*setupsort)(EPS);
  PetscErrorCode (*setfromoptions)(EPS,PetscOptionItems*);
  PetscErrorCode (*publishoptions)(EPS);
  PetscErrorCode (*destroy)(EPS);
  PetscErrorCode (*reset)(EPS);
  PetscErrorCode (*view)(EPS,PetscViewer);
  PetscErrorCode (*backtransform)(EPS);
  PetscErrorCode (*computevectors)(EPS);
  PetscErrorCode (*setdefaultst)(EPS);
  PetscErrorCode (*setdstype)(EPS);
};

/*
   Maximum number of monitors you can run with a single EPS
*/
#define MAXEPSMONITORS 5

/*
   The solution process goes through several states
*/
typedef enum { EPS_STATE_INITIAL,
               EPS_STATE_SETUP,
               EPS_STATE_SOLVED,
               EPS_STATE_EIGENVECTORS } EPSStateType;

/*
   To classify the different solvers into categories
*/
typedef enum { EPS_CATEGORY_KRYLOV,      /* Krylov solver: relies on STApply and STBackTransform (same as OTHER) */
               EPS_CATEGORY_PRECOND,     /* Preconditioned solver: uses ST only to manage preconditioner */
               EPS_CATEGORY_CONTOUR,     /* Contour integral: ST used to solve linear systems at integration points */
               EPS_CATEGORY_OTHER } EPSSolverType;

/*
   To check for unsupported features at EPSSetUp_XXX()
*/
typedef enum { EPS_FEATURE_BALANCE=1,       /* balancing */
               EPS_FEATURE_ARBITRARY=2,     /* arbitrary selection of eigepairs */
               EPS_FEATURE_REGION=4,        /* nontrivial region for filtering */
               EPS_FEATURE_EXTRACTION=8,    /* extraction technique different from Ritz */
               EPS_FEATURE_CONVERGENCE=16,  /* convergence test selected by user */
               EPS_FEATURE_STOPPING=32,     /* stopping test */
               EPS_FEATURE_THRESHOLD=64,    /* threshold stopping test */
               EPS_FEATURE_TWOSIDED=128     /* two-sided variant */
             } EPSFeatureType;

/*
   Defines the EPS data structure
*/
struct _p_EPS {
  PETSCHEADER(struct _EPSOps);
  /*------------------------- User parameters ---------------------------*/
  PetscInt       max_it;           /* maximum number of iterations */
  PetscInt       nev;              /* number of eigenvalues to compute */
  PetscInt       ncv;              /* number of basis vectors */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini,ninil;       /* number of initial vectors (negative means not copied yet) */
  PetscInt       nds;              /* number of basis vectors of deflation space */
  PetscScalar    target;           /* target value */
  PetscReal      tol;              /* tolerance */
  PetscReal      thres;            /* threshold */
  PetscBool      threlative;       /* threshold is relative */
  EPSConv        conv;             /* convergence test */
  EPSStop        stop;             /* stopping test */
  EPSWhich       which;            /* which part of the spectrum to be sought */
  PetscReal      inta,intb;        /* interval [a,b] for spectrum slicing */
  EPSProblemType problem_type;     /* which kind of problem to be solved */
  EPSExtraction  extraction;       /* which kind of extraction to be applied */
  EPSBalance     balance;          /* the balancing method */
  PetscInt       balance_its;      /* number of iterations of the balancing method */
  PetscReal      balance_cutoff;   /* cutoff value for balancing */
  PetscBool      trueres;          /* whether the true residual norm must be computed */
  PetscBool      trackall;         /* whether all the residuals must be computed */
  PetscBool      purify;           /* whether eigenvectors need to be purified */
  PetscBool      twosided;         /* whether to compute left eigenvectors (two-sided solver) */

  /*-------------- User-provided functions and contexts -----------------*/
  EPSConvergenceTestFn      *converged;
  EPSConvergenceTestFn      *convergeduser;
  PetscCtxDestroyFn         *convergeddestroy;
  EPSStoppingTestFn         *stopping;
  EPSStoppingTestFn         *stoppinguser;
  PetscCtxDestroyFn         *stoppingdestroy;
  SlepcArbitrarySelectionFn *arbitrary;
  void                      *convergedctx;
  void                      *stoppingctx;
  void                      *arbitraryctx;
  PetscErrorCode            (*monitor[MAXEPSMONITORS])(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscCtxDestroyFn         *monitordestroy[MAXEPSMONITORS];
  void                      *monitorcontext[MAXEPSMONITORS];
  PetscInt                  numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  ST             st;               /* spectral transformation object */
  DS             ds;               /* direct solver object */
  BV             V;                /* set of basis vectors and computed eigenvectors */
  BV             W;                /* left basis vectors (if left eigenvectors requested) */
  RG             rg;               /* optional region for filtering */
  SlepcSC        sc;               /* sorting criterion data */
  Vec            D;                /* diagonal matrix for balancing */
  Vec            *IS,*ISL;         /* references to user-provided initial spaces */
  Vec            *defl;            /* references to user-provided deflation space */
  PetscScalar    *eigr,*eigi;      /* real and imaginary parts of eigenvalues */
  PetscReal      *errest;          /* error estimates */
  PetscScalar    *rr,*ri;          /* values computed by user's arbitrary selection function */
  PetscInt       *perm;            /* permutation for eigenvalue ordering */
  PetscInt       nwork;            /* number of work vectors */
  Vec            *work;            /* work vectors */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables --------------------------*/
  EPSStateType   state;            /* initial -> setup -> solved -> eigenvectors */
  EPSSolverType  categ;            /* solver category */
  PetscInt       nconv;            /* number of converged eigenvalues */
  PetscInt       its;              /* number of iterations so far computed */
  PetscInt       n,nloc;           /* problem dimensions (global, local) */
  PetscReal      nrma,nrmb;        /* computed matrix norms */
  PetscBool      useds;            /* whether the solver uses the DS object or not */
  PetscBool      isgeneralized;
  PetscBool      ispositive;
  PetscBool      ishermitian;
  PetscBool      isstructured;
  EPSConvergedReason reason;
};

/*
    Macros to test valid EPS arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define EPSCheckSolved(h,arg) do {(void)(h);} while (0)

#else

#define EPSCheckSolved(h,arg) \
  do { \
    PetscCheck((h)->state>=EPS_STATE_SOLVED,PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

/*
    Macros to check settings at EPSSetUp()
*/

/* EPSCheckHermitianDefinite: the problem is HEP or GHEP */
#define EPSCheckHermitianDefiniteCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck((eps)->ishermitian,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for non-%s problems",((PetscObject)(eps))->type_name,(msg),SLEPC_STRING_HERMITIAN); \
      PetscCheck(!(eps)->isgeneralized || (eps)->ispositive,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s requires that the problem is %s-definite",((PetscObject)(eps))->type_name,(msg),SLEPC_STRING_HERMITIAN); \
    } \
  } while (0)
#define EPSCheckHermitianDefinite(eps) EPSCheckHermitianDefiniteCondition(eps,PETSC_TRUE,"")

/* EPSCheckHermitian: the problem is HEP, GHEP, or GHIEP */
#define EPSCheckHermitianCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck((eps)->ishermitian,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for non-%s problems",((PetscObject)(eps))->type_name,(msg),SLEPC_STRING_HERMITIAN); \
    } \
  } while (0)
#define EPSCheckHermitian(eps) EPSCheckHermitianCondition(eps,PETSC_TRUE,"")

/* EPSCheckDefinite: the problem is not GHIEP */
#define EPSCheckDefiniteCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!(eps)->isgeneralized || !(eps)->ishermitian || (eps)->ispositive,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for %s-indefinite problems",((PetscObject)(eps))->type_name,(msg),SLEPC_STRING_HERMITIAN); \
    } \
  } while (0)
#define EPSCheckDefinite(eps) EPSCheckDefiniteCondition(eps,PETSC_TRUE,"")

/* EPSCheckStandard: the problem is HEP or NHEP */
#define EPSCheckStandardCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!(eps)->isgeneralized,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s cannot be used for generalized problems",((PetscObject)(eps))->type_name,(msg)); \
    } \
  } while (0)
#define EPSCheckStandard(eps) EPSCheckStandardCondition(eps,PETSC_TRUE,"")

/* EPSCheckNotStructured: the problem is not structured */
#define EPSCheckNotStructuredCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!(eps)->isstructured,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s does not provide support for structured eigenproblems",((PetscObject)(eps))->type_name,(msg)); \
    } \
  } while (0)
#define EPSCheckNotStructured(eps) EPSCheckNotStructuredCondition(eps,PETSC_TRUE,"")

/* EPSCheckSinvert: shift-and-invert ST */
#define EPSCheckSinvertCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscBool __flg; \
      PetscCall(PetscObjectTypeCompare((PetscObject)(eps)->st,STSINVERT,&__flg)); \
      PetscCheck(__flg,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s requires a shift-and-invert spectral transform",((PetscObject)(eps))->type_name,(msg)); \
    } \
  } while (0)
#define EPSCheckSinvert(eps) EPSCheckSinvertCondition(eps,PETSC_TRUE,"")

/* EPSCheckSinvertCayley: shift-and-invert or Cayley ST */
#define EPSCheckSinvertCayleyCondition(eps,condition,msg) \
  do { \
    if (condition) { \
      PetscBool __flg; \
      PetscCall(PetscObjectTypeCompareAny((PetscObject)(eps)->st,&__flg,STSINVERT,STCAYLEY,"")); \
      PetscCheck(__flg,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s requires shift-and-invert or Cayley transform",((PetscObject)(eps))->type_name,(msg)); \
    } \
  } while (0)
#define EPSCheckSinvertCayley(eps) EPSCheckSinvertCayleyCondition(eps,PETSC_TRUE,"")

/* Check for unsupported features */
#define EPSCheckUnsupportedCondition(eps,mask,condition,msg) \
  do { \
    if (condition) { \
      PetscCheck(!((mask) & EPS_FEATURE_BALANCE) || (eps)->balance==EPS_BALANCE_NONE,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s does not support balancing",((PetscObject)(eps))->type_name,(msg)); \
      PetscCheck(!((mask) & EPS_FEATURE_ARBITRARY) || !(eps)->arbitrary,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s does not support arbitrary selection of eigenpairs",((PetscObject)(eps))->type_name,(msg)); \
      if ((mask) & EPS_FEATURE_REGION) { \
        PetscBool      __istrivial; \
        PetscCall(RGIsTrivial((eps)->rg,&__istrivial)); \
        PetscCheck(__istrivial,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s does not support region filtering",((PetscObject)(eps))->type_name,(msg)); \
      } \
      PetscCheck(!((mask) & EPS_FEATURE_EXTRACTION) || (eps)->extraction==EPS_RITZ,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s only supports Ritz extraction",((PetscObject)(eps))->type_name,(msg)); \
      PetscCheck(!((mask) & EPS_FEATURE_CONVERGENCE) || (eps)->converged==EPSConvergedRelative,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default convergence test",((PetscObject)(eps))->type_name,(msg)); \
      PetscCheck(!((mask) & EPS_FEATURE_STOPPING) || (eps)->stopping==EPSStoppingBasic,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s only supports the default stopping test",((PetscObject)(eps))->type_name,(msg)); \
      PetscCheck(!((mask) & EPS_FEATURE_THRESHOLD) || (eps)->stopping!=EPSStoppingThreshold,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s does not support the threshold stopping test",((PetscObject)(eps))->type_name,(msg)); \
      PetscCheck(!((mask) & EPS_FEATURE_TWOSIDED) || !(eps)->twosided,PetscObjectComm((PetscObject)(eps)),PETSC_ERR_SUP,"The solver '%s'%s cannot compute left eigenvectors (no two-sided variant)",((PetscObject)(eps))->type_name,(msg)); \
    } \
  } while (0)
#define EPSCheckUnsupported(eps,mask) EPSCheckUnsupportedCondition(eps,mask,PETSC_TRUE,"")

/* Check for ignored features */
#define EPSCheckIgnoredCondition(eps,mask,condition,msg) \
  do { \
    if (condition) { \
      if (((mask) & EPS_FEATURE_BALANCE) && (eps)->balance!=EPS_BALANCE_NONE) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the balancing settings\n",((PetscObject)(eps))->type_name,(msg))); \
      if (((mask) & EPS_FEATURE_ARBITRARY) && (eps)->arbitrary) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the settings for arbitrary selection of eigenpairs\n",((PetscObject)(eps))->type_name,(msg))); \
      if ((mask) & EPS_FEATURE_REGION) { \
        PetscBool __istrivial; \
        PetscCall(RGIsTrivial((eps)->rg,&__istrivial)); \
        if (!__istrivial) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the specified region\n",((PetscObject)(eps))->type_name,(msg))); \
      } \
      if (((mask) & EPS_FEATURE_EXTRACTION) && (eps)->extraction!=EPS_RITZ) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the extraction settings\n",((PetscObject)(eps))->type_name,(msg))); \
      if (((mask) & EPS_FEATURE_CONVERGENCE) && (eps)->converged!=EPSConvergedRelative) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the convergence test settings\n",((PetscObject)(eps))->type_name,(msg))); \
      if (((mask) & EPS_FEATURE_STOPPING) && (eps)->stopping!=EPSStoppingBasic) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the stopping test settings\n",((PetscObject)(eps))->type_name,(msg))); \
      if (((mask) & EPS_FEATURE_TWOSIDED) && (eps)->twosided) PetscCall(PetscInfo((eps),"The solver '%s'%s ignores the two-sided flag\n",((PetscObject)(eps))->type_name,(msg))); \
    } \
  } while (0)
#define EPSCheckIgnored(eps,mask) EPSCheckIgnoredCondition(eps,mask,PETSC_TRUE,"")

/*
    EPSSetCtxThreshold - Fills EPSStoppingCtx with data needed for the threshold stopping test
*/
#define EPSSetCtxThreshold(eps,eigr,eigi,k) \
  do { \
    if (eps->stop==EPS_STOP_THRESHOLD && k) { \
      PetscScalar __kr=eigr[k-1],__ki=eigi[k-1],__kr0=eigr[0],__ki0=eigi[0]; \
      PetscCall(STBackTransform(eps->st,1,&__kr,&__ki)); \
      PetscCall(STBackTransform(eps->st,1,&__kr0,&__ki0)); \
      if (eps->which==EPS_LARGEST_MAGNITUDE || eps->which==EPS_SMALLEST_MAGNITUDE) { \
        ((EPSStoppingCtx)eps->stoppingctx)->firstev = SlepcAbsEigenvalue(__kr0,__ki0); \
        ((EPSStoppingCtx)eps->stoppingctx)->lastev  = SlepcAbsEigenvalue(__kr,__ki); \
      } else { \
        ((EPSStoppingCtx)eps->stoppingctx)->firstev = PetscRealPart(__kr0); \
        ((EPSStoppingCtx)eps->stoppingctx)->lastev  = PetscRealPart(__kr); \
      } \
    } \
  } while (0)

/*
  EPS_SetInnerProduct - set B matrix for inner product if appropriate.
*/
static inline PetscErrorCode EPS_SetInnerProduct(EPS eps)
{
  Mat            B;

  PetscFunctionBegin;
  if (!eps->V) PetscCall(EPSGetBV(eps,&eps->V));
  if (eps->ispositive || (eps->isgeneralized && eps->ishermitian)) {
    PetscCall(STGetBilinearForm(eps->st,&B));
    PetscCall(BVSetMatrix(eps->V,B,PetscNot(eps->ispositive)));
    PetscCall(MatDestroy(&B));
  } else PetscCall(BVSetMatrix(eps->V,NULL,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPS_Purify - purify the first k vectors in the V basis
*/
static inline PetscErrorCode EPS_Purify(EPS eps,PetscInt k)
{
  PetscInt       i;
  Vec            v,z;

  PetscFunctionBegin;
  PetscCall(BVCreateVec(eps->V,&v));
  for (i=0;i<k;i++) {
    PetscCall(BVCopyVec(eps->V,i,v));
    PetscCall(BVGetColumn(eps->V,i,&z));
    PetscCall(STApply(eps->st,v,z));
    PetscCall(BVRestoreColumn(eps->V,i,&z));
  }
  PetscCall(VecDestroy(&v));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPS_KSPSetOperators - Sets the KSP matrices, see also ST_KSPSetOperators()
*/
static inline PetscErrorCode EPS_KSPSetOperators(KSP ksp,Mat A,Mat B)
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
  EPS_GetActualConverged - Gets the actual value of nconv; in special cases the
  number of available eigenvalues is larger than the computed ones
*/
static inline PetscErrorCode EPS_GetActualConverged(EPS eps,PetscInt *nconv)
{
  PetscFunctionBegin;
  *nconv = eps->nconv;
  if (eps->isstructured) {
    if (eps->problem_type == EPS_BSE && (eps->which == EPS_SMALLEST_MAGNITUDE || eps->which == EPS_LARGEST_MAGNITUDE || eps->which == EPS_TARGET_MAGNITUDE)) *nconv *= 2;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  EPS_GetEigenvector - Gets the i-th eigenvector taking into account the case
  where i exceeds the number of computed vectors (structure-preserving solver).
  The argument V should be eps->V for right eigenvectors, eps->W for left ones.
*/
static inline PetscErrorCode EPS_GetEigenvector(EPS eps,BV V,PetscInt i,Vec Vr,Vec Vi)
{
  PetscInt k;
  Vec      v0,v1,w,w0,w1;
  Mat      H;
  IS       is[2];

  PetscFunctionBegin;
  if (!eps->isstructured) {
    k = eps->perm[i];
    PetscCall(BV_GetEigenvector(V,k,eps->eigi[k],Vr,Vi));
  } else {
    if (eps->problem_type == EPS_BSE && (eps->which == EPS_SMALLEST_MAGNITUDE || eps->which == EPS_LARGEST_MAGNITUDE || eps->which == EPS_TARGET_MAGNITUDE)) {
      /* BSE problem, even index is +lambda, odd index is -lambda */
      k = eps->perm[i/2];
      if (i%2) {
        /* eigenvector of -lambda is J*conj(X) where J=[0 I; I 0] and x is eigenvector of lambda */
        PetscCall(VecDuplicate(Vr?Vr:Vi,&w));
        PetscCall(STGetMatrix(eps->st,0,&H));
        PetscCall(MatNestGetISs(H,is,NULL));
        if (Vr) {
          PetscCall(BV_GetEigenvector(V,k,eps->eigi[k],w,NULL));
          PetscCall(VecConjugate(w));
          PetscCall(VecGetSubVector(w,is[0],&w0));
          PetscCall(VecGetSubVector(w,is[1],&w1));
          PetscCall(VecGetSubVector(Vr,is[0],&v0));
          PetscCall(VecGetSubVector(Vr,is[1],&v1));
          PetscCall(VecCopy(w1,v0));
          PetscCall(VecCopy(w0,v1));
          PetscCall(VecRestoreSubVector(w,is[0],&w0));
          PetscCall(VecRestoreSubVector(w,is[1],&w1));
          PetscCall(VecRestoreSubVector(Vr,is[0],&v0));
          PetscCall(VecRestoreSubVector(Vr,is[1],&v1));
        }
#if !defined(PETSC_USE_COMPLEX)
        if (Vi) {
          PetscCall(BV_GetEigenvector(V,k,eps->eigi[k],NULL,w));
          PetscCall(VecScale(w,-1.0));
          PetscCall(VecGetSubVector(w,is[0],&w0));
          PetscCall(VecGetSubVector(w,is[1],&w1));
          PetscCall(VecGetSubVector(Vi,is[0],&v0));
          PetscCall(VecGetSubVector(Vi,is[1],&v1));
          PetscCall(VecCopy(w1,v0));
          PetscCall(VecCopy(w0,v1));
          PetscCall(VecRestoreSubVector(w,is[0],&w0));
          PetscCall(VecRestoreSubVector(w,is[1],&w1));
          PetscCall(VecRestoreSubVector(Vi,is[0],&v0));
          PetscCall(VecRestoreSubVector(Vi,is[1],&v1));
        }
#endif
        PetscCall(VecDestroy(&w));
      } else {
        PetscCall(BV_GetEigenvector(V,k,eps->eigi[k],Vr,Vi));
      }
    } else SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Inconsistent state");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_INTERN PetscErrorCode EPSSetWhichEigenpairs_Default(EPS);
SLEPC_INTERN PetscErrorCode EPSSetDimensions_Default(EPS,PetscInt*,PetscInt*,PetscInt*);
SLEPC_INTERN PetscErrorCode EPSBackTransform_Default(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_Hermitian(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_Schur(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_Indefinite(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_Twosided(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_Slice(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeResidualNorm_Private(EPS,PetscBool,PetscScalar,PetscScalar,Vec,Vec,Vec*,PetscReal*);
SLEPC_INTERN PetscErrorCode EPSComputeRitzVector(EPS,PetscScalar*,PetscScalar*,BV,Vec,Vec);
SLEPC_INTERN PetscErrorCode EPSGetStartVector(EPS,PetscInt,PetscBool*);
SLEPC_INTERN PetscErrorCode EPSGetLeftStartVector(EPS,PetscInt,PetscBool*);
SLEPC_INTERN PetscErrorCode MatEstimateSpectralRange_EPS(Mat,PetscReal*,PetscReal*);

/* Private functions of the solver implementations */

SLEPC_INTERN PetscErrorCode EPSDelayedArnoldi(EPS,PetscScalar*,PetscInt,PetscInt,PetscInt*,PetscReal*,PetscBool*);
SLEPC_INTERN PetscErrorCode EPSDelayedArnoldi1(EPS,PetscScalar*,PetscInt,PetscInt,PetscInt*,PetscReal*,PetscBool*);
SLEPC_INTERN PetscErrorCode EPSKrylovConvergence(EPS,PetscBool,PetscInt,PetscInt,PetscReal,PetscReal,PetscReal,PetscInt*);
SLEPC_INTERN PetscErrorCode EPSPseudoLanczos(EPS,PetscReal*,PetscReal*,PetscReal*,PetscInt,PetscInt*,PetscBool*,PetscBool*,PetscReal*,Vec);
SLEPC_INTERN PetscErrorCode EPSBuildBalance_Krylov(EPS);
SLEPC_INTERN PetscErrorCode EPSSetDefaultST(EPS);
SLEPC_INTERN PetscErrorCode EPSSetDefaultST_Precond(EPS);
SLEPC_INTERN PetscErrorCode EPSSetDefaultST_GMRES(EPS);
SLEPC_INTERN PetscErrorCode EPSSetDefaultST_NoFactor(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUpSort_Basic(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUpSort_Default(EPS);
