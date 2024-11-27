/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to eigensolvers in ELPA.
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include <slepc/private/slepcscalapack.h>
#include <elpa/elpa.h>

#define PetscCallELPA(func, ...) do {                                                   \
    PetscErrorCode elpa_ierr_; \
    PetscStackPushExternal(PetscStringize(func));                                   \
    func(__VA_ARGS__,&elpa_ierr_);                                              \
    PetscStackPop;                                                                             \
    PetscCheck(!elpa_ierr_,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error calling %s: error code %d",PetscStringize(func(__VA_ARGS__,&elpa_ierr)),elpa_ierr_); \
  } while (0)

#define PetscCallELPARET(func, ...) do {                                                   \
    PetscStackPushExternal(PetscStringize(func));                                                      \
    PetscErrorCode elpa_ierr_ = func(__VA_ARGS__);                                              \
    PetscStackPop;                                                                             \
    PetscCheck(!elpa_ierr_,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error calling %s: error code %d",PetscStringize(func(__VA_ARGS__)),elpa_ierr_); \
  } while (0)

#define PetscCallELPANOARG(func) do {                                                   \
    PetscErrorCode elpa_ierr_; \
    PetscStackPushExternal(PetscStringize(func));                                   \
    func(&elpa_ierr_);                                              \
    PetscStackPop;                                                                             \
    PetscCheck(!elpa_ierr_,PETSC_COMM_SELF,PETSC_ERR_LIB,"Error calling %s: error code %d",PetscStringize(func(&elpa_ierr)),elpa_ierr_); \
  } while (0)

typedef struct {
  Mat As,Bs;        /* converted matrices */
} EPS_ELPA;

static PetscErrorCode EPSSetUp_ELPA(EPS eps)
{
  EPS_ELPA       *ctx = (EPS_ELPA*)eps->data;
  Mat            A,B;
  PetscInt       nmat;
  PetscBool      isshift;
  PetscScalar    shift;

  PetscFunctionBegin;
  EPSCheckHermitianDefinite(eps);
  EPSCheckNotStructured(eps);
  PetscCall(PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isshift));
  PetscCheck(isshift,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support spectral transformations");
  if (eps->nev==0) eps->nev = 1;
  eps->ncv = eps->n;
  if (eps->mpd!=PETSC_DETERMINE) PetscCall(PetscInfo(eps,"Warning: parameter mpd ignored\n"));
  if (eps->max_it==PETSC_DETERMINE) eps->max_it = 1;
  if (!eps->which) PetscCall(EPSSetWhichEigenpairs_Default(eps));
  PetscCheck(eps->which!=EPS_ALL || eps->inta==eps->intb,PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support interval computation");
  EPSCheckUnsupported(eps,EPS_FEATURE_BALANCE | EPS_FEATURE_ARBITRARY | EPS_FEATURE_REGION);
  EPSCheckIgnored(eps,EPS_FEATURE_EXTRACTION | EPS_FEATURE_CONVERGENCE | EPS_FEATURE_STOPPING);
  PetscCall(EPSAllocateSolution(eps,0));

  /* convert matrices */
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatDestroy(&ctx->Bs));
  PetscCall(STGetNumMatrices(eps->st,&nmat));
  PetscCall(STGetMatrix(eps->st,0,&A));
  PetscCall(MatConvert(A,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->As));
  if (nmat>1) {
    PetscCall(STGetMatrix(eps->st,1,&B));
    PetscCall(MatConvert(B,MATSCALAPACK,MAT_INITIAL_MATRIX,&ctx->Bs));
  }
  PetscCall(STGetShift(eps->st,&shift));
  if (shift != 0.0) {
    if (nmat>1) PetscCall(MatAXPY(ctx->As,-shift,ctx->Bs,SAME_NONZERO_PATTERN));
    else PetscCall(MatShift(ctx->As,-shift));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSSolve_ELPA(EPS eps)
{
  EPS_ELPA       *ctx = (EPS_ELPA*)eps->data;
  Mat            A = ctx->As,B = ctx->Bs,Q,V;
  Mat_ScaLAPACK  *a = (Mat_ScaLAPACK*)A->data,*b,*q;
  PetscReal      *w = eps->errest;  /* used to store real eigenvalues */
  PetscInt       i;
  elpa_t         handle;

  PetscFunctionBegin;
  PetscCall(MatDuplicate(A,MAT_DO_NOT_COPY_VALUES,&Q));
  q = (Mat_ScaLAPACK*)Q->data;

  PetscCallELPARET(elpa_init,20200417);    /* 20171201 */
  PetscCallELPANOARG(handle = elpa_allocate);

  /* set parameters of the matrix and its MPI distribution */
  PetscCallELPA(elpa_set,handle,"na",a->N);                         /* matrix size */
  PetscCallELPA(elpa_set,handle,"nev",a->N);                        /* number of eigenvectors to be computed (1<=nev<=na) */
  PetscCallELPA(elpa_set,handle,"local_nrows",a->locr);             /* number of local rows of the distributed matrix on this MPI task  */
  PetscCallELPA(elpa_set,handle,"local_ncols",a->locc);             /* number of local columns of the distributed matrix on this MPI task */
  PetscCallELPA(elpa_set,handle,"nblk",a->mb);                      /* size of the BLACS block cyclic distribution */
  PetscCallELPA(elpa_set,handle,"mpi_comm_parent",MPI_Comm_c2f(PetscObjectComm((PetscObject)eps)));
  PetscCallELPA(elpa_set,handle,"process_row",a->grid->myrow);      /* row coordinate of MPI process */
  PetscCallELPA(elpa_set,handle,"process_col",a->grid->mycol);      /* column coordinate of MPI process */
  if (B) PetscCallELPA(elpa_set,handle,"blacs_context",a->grid->ictxt);

  /* setup and set tunable run-time options */
  PetscCallELPARET(elpa_setup,handle);
  PetscCallELPA(elpa_set,handle,"solver",ELPA_SOLVER_2STAGE);
  /* PetscCallELPA(elpa_print_settings,handle); */

  /* solve the eigenvalue problem */
  if (B) {
    b = (Mat_ScaLAPACK*)B->data;
    PetscCallELPA(elpa_generalized_eigenvectors,handle,a->loc,b->loc,w,q->loc,0);
  } else PetscCallELPA(elpa_eigenvectors,handle,a->loc,w,q->loc);

  /* cleanup */
  PetscCallELPA(elpa_deallocate,handle);
  PetscCallELPANOARG(elpa_uninit);

  for (i=0;i<eps->ncv;i++) {
    eps->eigr[i]   = eps->errest[i];
    eps->errest[i] = PETSC_MACHINE_EPSILON;
  }

  PetscCall(BVGetMat(eps->V,&V));
  PetscCall(MatConvert(Q,MATDENSE,MAT_REUSE_MATRIX,&V));
  PetscCall(BVRestoreMat(eps->V,&V));
  PetscCall(MatDestroy(&Q));

  eps->nconv  = eps->ncv;
  eps->its    = 1;
  eps->reason = EPS_CONVERGED_TOL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSDestroy_ELPA(EPS eps)
{
  PetscFunctionBegin;
  PetscCall(PetscFree(eps->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode EPSReset_ELPA(EPS eps)
{
  EPS_ELPA       *ctx = (EPS_ELPA*)eps->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&ctx->As));
  PetscCall(MatDestroy(&ctx->Bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_ELPA(EPS eps)
{
  EPS_ELPA       *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->solve          = EPSSolve_ELPA;
  eps->ops->setup          = EPSSetUp_ELPA;
  eps->ops->setupsort      = EPSSetUpSort_Basic;
  eps->ops->destroy        = EPSDestroy_ELPA;
  eps->ops->reset          = EPSReset_ELPA;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_NoFactor;
  PetscFunctionReturn(PETSC_SUCCESS);
}
