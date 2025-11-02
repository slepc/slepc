/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "stoar"

   Method: S-TOAR with spectrum slicing for symmetric quadratic eigenproblems

   Algorithm:

       Symmetric Two-Level Orthogonal Arnoldi.

   References:

       [1] C. Campos and J.E. Roman, "Inertia-based spectrum slicing
           for symmetric quadratic eigenvalue problems", Numer. Linear
           Algebra Appl. 27(4):e2293, 2020.
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-slice-qep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Inertia-based spectrum slicing for symmetric quadratic eigenvalue problems\",\n"
  "   journal = \"Numer. Linear Algebra Appl.\",\n"
  "   volume = \"27\",\n"
  "   number = \"4\",\n"
  "   pages = \"e2293\",\n"
  "   year = \"2020,\"\n"
  "   doi = \"https://doi.org/10.1002/nla.2293\"\n"
  "}\n";

#define SLICE_PTOL PETSC_SQRT_MACHINE_EPSILON

static PetscErrorCode PEPQSliceResetSR(PEP pep)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;
  PEP_shift      s;
  PetscInt       i;

  PetscFunctionBegin;
  if (sr) {
    /* Reviewing list of shifts to free memory */
    s = sr->s0;
    if (s) {
      while (s->neighb[1]) {
        s = s->neighb[1];
        PetscCall(PetscFree(s->neighb[0]));
      }
      PetscCall(PetscFree(s));
    }
    PetscCall(PetscFree(sr->S));
    for (i=0;i<pep->nconv;i++) PetscCall(PetscFree(sr->qinfo[i].q));
    PetscCall(PetscFree(sr->qinfo));
    for (i=0;i<3;i++) PetscCall(VecDestroy(&sr->v[i]));
    PetscCall(EPSDestroy(&sr->eps));
    PetscCall(PetscFree(sr));
  }
  ctx->sr = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPReset_STOAR_QSlice(PEP pep)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  PetscCall(PEPQSliceResetSR(pep));
  PetscCall(PetscFree(ctx->inertias));
  PetscCall(PetscFree(ctx->shifts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  PEPQSliceAllocateSolution - Allocate memory storage for common variables such
  as eigenvalues and eigenvectors.
*/
static PetscErrorCode PEPQSliceAllocateSolution(PEP pep)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PetscInt       k;
  BVType         type;
  Vec            t;
  PEP_SR         sr = ctx->sr;

  PetscFunctionBegin;
  /* allocate space for eigenvalues and friends */
  k = PetscMax(1,sr->numEigs);
  PetscCall(PetscFree4(sr->eigr,sr->eigi,sr->errest,sr->perm));
  PetscCall(PetscCalloc4(k,&sr->eigr,k,&sr->eigi,k,&sr->errest,k,&sr->perm));
  PetscCall(PetscFree(sr->qinfo));
  PetscCall(PetscCalloc1(k,&sr->qinfo));

  /* allocate sr->V and transfer options from pep->V */
  PetscCall(BVDestroy(&sr->V));
  PetscCall(BVCreate(PetscObjectComm((PetscObject)pep),&sr->V));
  if (!pep->V) PetscCall(PEPGetBV(pep,&pep->V));
  if (!((PetscObject)pep->V)->type_name) PetscCall(BVSetType(sr->V,BVMAT));
  else {
    PetscCall(BVGetType(pep->V,&type));
    PetscCall(BVSetType(sr->V,type));
  }
  PetscCall(STMatCreateVecsEmpty(pep->st,&t,NULL));
  PetscCall(BVSetSizesFromVec(sr->V,t,k+1));
  PetscCall(VecDestroy(&t));
  sr->ld = k;
  PetscCall(PetscFree(sr->S));
  PetscCall(PetscMalloc1((k+1)*sr->ld*(pep->nmat-1),&sr->S));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Convergence test to compute positive Ritz values */
static PetscErrorCode ConvergedPositive(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = (PetscRealPart(eigr)>0.0)?0.0:res;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQSliceMatGetInertia(PEP pep,PetscReal shift,PetscInt *inertia,PetscInt *zeros)
{
  KSP            ksp,kspr;
  PC             pc;
  Mat            F;
  PetscBool      flg;

  PetscFunctionBegin;
  if (!pep->solvematcoeffs) PetscCall(PetscMalloc1(pep->nmat,&pep->solvematcoeffs));
  if (shift==PETSC_MAX_REAL) { /* Inertia of matrix A[2] */
    pep->solvematcoeffs[0] = 0.0; pep->solvematcoeffs[1] = 0.0; pep->solvematcoeffs[2] = 1.0;
  } else PetscCall(PEPEvaluateBasis(pep,shift,0,pep->solvematcoeffs,NULL));
  PetscCall(STMatSetUp(pep->st,pep->sfactor,pep->solvematcoeffs));
  PetscCall(STGetKSP(pep->st,&ksp));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCREDUNDANT,&flg));
  if (flg) {
    PetscCall(PCRedundantGetKSP(pc,&kspr));
    PetscCall(KSPGetPC(kspr,&pc));
  }
  PetscCall(PCFactorGetMatrix(pc,&F));
  PetscCall(MatGetInertia(F,inertia,zeros,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQSliceGetInertia(PEP pep,PetscReal shift,PetscInt *inertia,PetscInt *zeros,PetscInt correction)
{
  KSP            ksp;
  Mat            P;
  PetscReal      nzshift=0.0,dot;
  PetscRandom    rand;
  PetscInt       nconv;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;

  PetscFunctionBegin;
  if (shift >= PETSC_MAX_REAL) { /* Right-open interval */
    *inertia = 0;
  } else if (shift <= PETSC_MIN_REAL) {
    *inertia = 0;
    if (zeros) *zeros = 0;
  } else {
    /* If the shift is zero, perturb it to a very small positive value.
       The goal is that the nonzero pattern is the same in all cases and reuse
       the symbolic factorizations */
    nzshift = (shift==0.0)? 10.0/PETSC_MAX_REAL: shift;
    PetscCall(PEPQSliceMatGetInertia(pep,nzshift,inertia,zeros));
    PetscCall(STSetShift(pep->st,nzshift));
  }
  if (!correction) {
    if (shift >= PETSC_MAX_REAL) *inertia = 2*pep->n;
    else if (shift>PETSC_MIN_REAL) {
      PetscCall(STGetKSP(pep->st,&ksp));
      PetscCall(KSPGetOperators(ksp,&P,NULL));
      if (*inertia!=pep->n && !sr->v[0]) {
        PetscCall(MatCreateVecs(P,&sr->v[0],NULL));
        PetscCall(VecDuplicate(sr->v[0],&sr->v[1]));
        PetscCall(VecDuplicate(sr->v[0],&sr->v[2]));
        PetscCall(BVGetRandomContext(pep->V,&rand));
        PetscCall(VecSetRandom(sr->v[0],rand));
      }
      if (*inertia<pep->n && *inertia>0) {
        if (!sr->eps) {
          PetscCall(EPSCreate(PetscObjectComm((PetscObject)pep),&sr->eps));
          PetscCall(EPSSetProblemType(sr->eps,EPS_HEP));
          PetscCall(EPSSetWhichEigenpairs(sr->eps,EPS_LARGEST_REAL));
        }
        PetscCall(EPSSetConvergenceTestFunction(sr->eps,ConvergedPositive,NULL,NULL));
        PetscCall(EPSSetOperators(sr->eps,P,NULL));
        PetscCall(EPSSolve(sr->eps));
        PetscCall(EPSGetConverged(sr->eps,&nconv));
        PetscCheck(nconv,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",(double)nzshift);
        PetscCall(EPSGetEigenpair(sr->eps,0,NULL,NULL,sr->v[0],sr->v[1]));
      }
      if (*inertia!=pep->n) {
        PetscCall(MatMult(pep->A[1],sr->v[0],sr->v[1]));
        PetscCall(MatMult(pep->A[2],sr->v[0],sr->v[2]));
        PetscCall(VecAXPY(sr->v[1],2*nzshift,sr->v[2]));
        PetscCall(VecDotRealPart(sr->v[1],sr->v[0],&dot));
        if (dot>0.0) *inertia = 2*pep->n-*inertia;
      }
    }
  } else if (correction<0) *inertia = 2*pep->n-*inertia;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Check eigenvalue type - used only in non-hyperbolic problems.
   All computed eigenvalues must have the same definite type (positive or negative).
   If ini=TRUE the type is available in omega, otherwise we compute an eigenvalue
   closest to shift and determine its type.
 */
static PetscErrorCode PEPQSliceCheckEigenvalueType(PEP pep,PetscReal shift,PetscReal omega,PetscBool ini)
{
  PEP            pep2;
  ST             st;
  PetscInt       nconv;
  PetscScalar    lambda;
  PetscReal      dot;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;

  PetscFunctionBegin;
  if (!ini) {
    PetscCheck(-(omega/(shift*ctx->alpha+ctx->beta))*sr->type>=0,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected in eigenvalue %g",(double)shift);
  } else {
    PetscCall(PEPCreate(PetscObjectComm((PetscObject)pep),&pep2));
    PetscCall(PEPSetOptionsPrefix(pep2,((PetscObject)pep)->prefix));
    PetscCall(PEPAppendOptionsPrefix(pep2,"pep_eigenvalue_type_"));
    PetscCall(PEPSetTolerances(pep2,PETSC_CURRENT,pep->max_it/4));
    PetscCall(PEPSetType(pep2,PEPTOAR));
    PetscCall(PEPSetOperators(pep2,pep->nmat,pep->A));
    PetscCall(PEPSetWhichEigenpairs(pep2,PEP_TARGET_MAGNITUDE));
    PetscCall(PEPGetRG(pep2,&pep2->rg));
    PetscCall(RGSetType(pep2->rg,RGINTERVAL));
#if defined(PETSC_USE_COMPLEX)
    PetscCall(RGIntervalSetEndpoints(pep2->rg,pep->inta,pep->intb,-PETSC_SQRT_MACHINE_EPSILON,PETSC_SQRT_MACHINE_EPSILON));
#else
    PetscCall(RGIntervalSetEndpoints(pep2->rg,pep->inta,pep->intb,0.0,0.0));
#endif
    pep2->target = shift;
    st = pep2->st;
    pep2->st = pep->st;
    PetscCall(PEPSolve(pep2));
    PetscCall(PEPGetConverged(pep2,&nconv));
    if (nconv) {
      PetscCall(PEPGetEigenpair(pep2,0,&lambda,NULL,pep2->work[0],NULL));
      PetscCall(MatMult(pep->A[1],pep2->work[0],pep2->work[1]));
      PetscCall(MatMult(pep->A[2],pep2->work[0],pep2->work[2]));
      PetscCall(VecAXPY(pep2->work[1],2.0*lambda*pep->sfactor,pep2->work[2]));
      PetscCall(VecDotRealPart(pep2->work[1],pep2->work[0],&dot));
      PetscCall(PetscInfo(pep,"lambda=%g, %s type\n",(double)PetscRealPart(lambda),(dot>0.0)?"positive":"negative"));
      if (!sr->type) sr->type = (dot>0.0)?1:-1;
      else PetscCheck(sr->type*dot>=0.0,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected in eigenvalue %g",(double)PetscRealPart(lambda));
    }
    pep2->st = st;
    PetscCall(PEPDestroy(&pep2));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode PEPQSliceDiscriminant(PEP pep,Vec u,Vec w,PetscReal *d,PetscReal *smas,PetscReal *smenos)
{
  PetscReal      ap,bp,cp,dis;

  PetscFunctionBegin;
  PetscCall(MatMult(pep->A[0],u,w));
  PetscCall(VecDotRealPart(w,u,&cp));
  PetscCall(MatMult(pep->A[1],u,w));
  PetscCall(VecDotRealPart(w,u,&bp));
  PetscCall(MatMult(pep->A[2],u,w));
  PetscCall(VecDotRealPart(w,u,&ap));
  dis = bp*bp-4*ap*cp;
  if (dis>=0.0 && smas) {
    if (ap>0) *smas = (-bp+PetscSqrtReal(dis))/(2*ap);
    else if (ap<0) *smas = (-bp-PetscSqrtReal(dis))/(2*ap);
    else {
      if (bp >0) *smas = -cp/bp;
      else *smas = PETSC_MAX_REAL;
    }
  }
  if (dis>=0.0 && smenos) {
    if (ap>0) *smenos = (-bp-PetscSqrtReal(dis))/(2*ap);
    else if (ap<0) *smenos = (-bp+PetscSqrtReal(dis))/(2*ap);
    else {
      if (bp<0) *smenos = -cp/bp;
      else *smenos = PETSC_MAX_REAL;
    }
  }
  if (d) *d = dis;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static inline PetscErrorCode PEPQSliceEvaluateQEP(PEP pep,PetscScalar x,Mat M,MatStructure str)
{
  PetscFunctionBegin;
  PetscCall(MatCopy(pep->A[0],M,SAME_NONZERO_PATTERN));
  PetscCall(MatAXPY(M,x,pep->A[1],str));
  PetscCall(MatAXPY(M,x*x,pep->A[2],str));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   PEPCheckDefiniteQEP - Determines if a symmetric/Hermitian quadratic eigenvalue problem
   is definite or not.

   Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  xi - first computed parameter
.  mu - second computed parameter
.  definite - flag indicating that the problem is definite
-  hyperbolic - flag indicating that the problem is hyperbolic

   Notes:
   This function is intended for quadratic eigenvalue problems, $Q(\lambda)=K+\lambda C+\lambda^2M$,
   with symmetric (or Hermitian) coefficient matrices $K$, $C$, $M$.

   On output, the flag `definite` may have the values -1 (meaning that the QEP is not
   definite), 1 (if the problem is definite), or 0 if the algorithm was not able to
   determine whether the problem is definite or not.

   If `definite`=1, the output flag `hyperbolic` informs in a similar way about whether the
   problem is hyperbolic or not.

   If `definite`=1, the computed values `xi` and `mu` satisfy $Q(\xi)<0$ and $Q(\mu)>0$,
   as obtained via the method proposed by {cite:t}`Nie10`. Furthermore, if
   `hyperbolic`=1 then only `xi` is computed.

   Level: advanced

.seealso: [](ch:pep), `PEPSetProblemType()`
@*/
PetscErrorCode PEPCheckDefiniteQEP(PEP pep,PetscReal *xi,PetscReal *mu,PetscInt *definite,PetscInt *hyperbolic)
{
  PetscRandom    rand;
  Vec            u,w;
  PetscReal      d=0.0,s=0.0,sp,mut=0.0,omg=0.0,omgp;
  PetscInt       k,its=10,hyp=0,check=0,nconv,inertia,n;
  Mat            M=NULL;
  MatStructure   str;
  EPS            eps;
  PetscBool      transform,ptypehyp;

  PetscFunctionBegin;
  PetscCheck(pep->problem_type==PEP_HERMITIAN || pep->problem_type==PEP_HYPERBOLIC,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Only available for Hermitian (or hyperbolic) problems");
  ptypehyp = (pep->problem_type==PEP_HYPERBOLIC)? PETSC_TRUE: PETSC_FALSE;
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(PEPSetDefaultST(pep));
  PetscCall(STSetMatrices(pep->st,pep->nmat,pep->A));
  PetscCall(MatGetSize(pep->A[0],&n,NULL));
  PetscCall(STGetTransform(pep->st,&transform));
  PetscCall(STSetTransform(pep->st,PETSC_FALSE));
  PetscCall(STSetUp(pep->st));
  PetscCall(MatCreateVecs(pep->A[0],&u,&w));
  PetscCall(PEPGetBV(pep,&pep->V));
  PetscCall(BVGetRandomContext(pep->V,&rand));
  PetscCall(VecSetRandom(u,rand));
  PetscCall(VecNormalize(u,NULL));
  PetscCall(PEPQSliceDiscriminant(pep,u,w,&d,&s,NULL));
  if (d<0.0) check = -1;
  if (!check) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)pep),&eps));
    PetscCall(EPSSetProblemType(eps,EPS_HEP));
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
    PetscCall(EPSSetTolerances(eps,PetscSqrtReal(PETSC_SQRT_MACHINE_EPSILON),PETSC_DECIDE));
    PetscCall(MatDuplicate(pep->A[0],MAT_DO_NOT_COPY_VALUES,&M));
    PetscCall(STGetMatStructure(pep->st,&str));
  }
  for (k=0;k<its&&!check;k++) {
    PetscCall(PEPQSliceEvaluateQEP(pep,s,M,str));
    PetscCall(EPSSetOperators(eps,M,NULL));
    PetscCall(EPSSolve(eps));
    PetscCall(EPSGetConverged(eps,&nconv));
    if (!nconv) break;
    PetscCall(EPSGetEigenpair(eps,0,NULL,NULL,u,w));
    sp = s;
    PetscCall(PEPQSliceDiscriminant(pep,u,w,&d,&s,&omg));
    if (d<0.0) {check = -1; break;}
    if (PetscAbsReal((s-sp)/s)<100*PETSC_MACHINE_EPSILON) break;
    if (s>sp) hyp = -1;
    mut = 2*s-sp;
    PetscCall(PEPQSliceMatGetInertia(pep,mut,&inertia,NULL));
    if (inertia == n) {check = 1; break;}
  }
  for (;k<its&&!check;k++) {
    mut = (s-omg)/2;
    PetscCall(PEPQSliceMatGetInertia(pep,mut,&inertia,NULL));
    if (inertia == n) {check = 1; break;}
    if (PetscAbsReal((s-omg)/omg)<100*PETSC_MACHINE_EPSILON) break;
    PetscCall(PEPQSliceEvaluateQEP(pep,omg,M,str));
    PetscCall(EPSSetOperators(eps,M,NULL));
    PetscCall(EPSSolve(eps));
    PetscCall(EPSGetConverged(eps,&nconv));
    if (!nconv) break;
    PetscCall(EPSGetEigenpair(eps,0,NULL,NULL,u,w));
    omgp = omg;
    PetscCall(PEPQSliceDiscriminant(pep,u,w,&d,NULL,&omg));
    if (d<0.0) {check = -1; break;}
    if (omg<omgp) hyp = -1;
  }
  if (check==1) *xi = mut;
  PetscCheck(hyp!=-1 || !ptypehyp,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Problem does not satisfy hyperbolic test; consider removing the hyperbolicity flag");
  if (check==1 && hyp==0) {
    PetscCall(PEPQSliceMatGetInertia(pep,PETSC_MAX_REAL,&inertia,NULL));
    if (inertia == 0) hyp = 1;
    else hyp = -1;
  }
  if (check==1 && hyp!=1) {
    check = 0;
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
    for (;k<its&&!check;k++) {
      PetscCall(PEPQSliceEvaluateQEP(pep,s,M,str));
      PetscCall(EPSSetOperators(eps,M,NULL));
      PetscCall(EPSSolve(eps));
      PetscCall(EPSGetConverged(eps,&nconv));
      if (!nconv) break;
      PetscCall(EPSGetEigenpair(eps,0,NULL,NULL,u,w));
      sp = s;
      PetscCall(PEPQSliceDiscriminant(pep,u,w,&d,&s,&omg));
      if (d<0.0) {check = -1; break;}
      if (PetscAbsReal((s-sp)/s)<100*PETSC_MACHINE_EPSILON) break;
      mut = 2*s-sp;
      PetscCall(PEPQSliceMatGetInertia(pep,mut,&inertia,NULL));
      if (inertia == 0) {check = 1; break;}
    }
    for (;k<its&&!check;k++) {
      mut = (s-omg)/2;
      PetscCall(PEPQSliceMatGetInertia(pep,mut,&inertia,NULL));
      if (inertia == 0) {check = 1; break;}
      if (PetscAbsReal((s-omg)/omg)<100*PETSC_MACHINE_EPSILON) break;
      PetscCall(PEPQSliceEvaluateQEP(pep,omg,M,str));
      PetscCall(EPSSetOperators(eps,M,NULL));
      PetscCall(EPSSolve(eps));
      PetscCall(EPSGetConverged(eps,&nconv));
      if (!nconv) break;
      PetscCall(EPSGetEigenpair(eps,0,NULL,NULL,u,w));
      PetscCall(PEPQSliceDiscriminant(pep,u,w,&d,NULL,&omg));
      if (d<0.0) {check = -1; break;}
    }
  }
  if (check==1) *mu = mut;
  *definite = check;
  *hyperbolic = hyp;
  if (M) PetscCall(MatDestroy(&M));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&w));
  PetscCall(EPSDestroy(&eps));
  PetscCall(STSetTransform(pep->st,transform));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Dummy backtransform operation
 */
static PetscErrorCode PEPBackTransform_Skip(PEP pep)
{
  PetscFunctionBegin;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPSetUp_STOAR_QSlice(PEP pep)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr;
  PetscInt       ld,i,zeros=0;
  SlepcSC        sc;
  PetscReal      r;

  PetscFunctionBegin;
  PEPCheckSinvertCayley(pep);
  PetscCheck(pep->inta<pep->intb,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"This solver does not support computing all eigenvalues unless you provide a computational interval with PEPSetInterval()");
  PetscCheck(pep->intb<PETSC_MAX_REAL || pep->inta>PETSC_MIN_REAL,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"The defined computational interval should have at least one of their sides bounded");
  PEPCheckUnsupportedCondition(pep,PEP_FEATURE_STOPPING,PETSC_TRUE," (with spectrum slicing)");
  if (pep->tol==(PetscReal)PETSC_DETERMINE) {
#if defined(PETSC_USE_REAL_SINGLE)
    pep->tol = SLEPC_DEFAULT_TOL;
#else
    /* use tighter tolerance */
    pep->tol = SLEPC_DEFAULT_TOL*1e-2;
#endif
  }
  if (ctx->nev==1) ctx->nev = PetscMin(20,pep->n);  /* nev not set, use default value */
  PetscCheck(pep->n<=10 || ctx->nev>=10,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"nev cannot be less than 10 in spectrum slicing runs");
  pep->ops->backtransform = PEPBackTransform_Skip;
  if (pep->max_it==PETSC_DETERMINE) pep->max_it = 100;

  /* create spectrum slicing context and initialize it */
  PetscCall(PEPQSliceResetSR(pep));
  PetscCall(PetscNew(&sr));
  ctx->sr   = sr;
  sr->itsKs = 0;
  sr->nleap = 0;
  sr->sPres = NULL;

  if (pep->solvematcoeffs) PetscCall(PetscFree(pep->solvematcoeffs));
  PetscCall(PetscMalloc1(pep->nmat,&pep->solvematcoeffs));
  if (!pep->st) PetscCall(PEPGetST(pep,&pep->st));
  PetscCall(STSetTransform(pep->st,PETSC_FALSE));
  PetscCall(STSetUp(pep->st));

  ctx->hyperbolic = (pep->problem_type==PEP_HYPERBOLIC)? PETSC_TRUE: PETSC_FALSE;

  /* check presence of ends and finding direction */
  if (pep->inta > PETSC_MIN_REAL || pep->intb >= PETSC_MAX_REAL) {
    sr->int0 = pep->inta;
    sr->int1 = pep->intb;
    sr->dir = 1;
    if (pep->intb >= PETSC_MAX_REAL) { /* Right-open interval */
      sr->hasEnd = PETSC_FALSE;
    } else sr->hasEnd = PETSC_TRUE;
  } else {
    sr->int0 = pep->intb;
    sr->int1 = pep->inta;
    sr->dir = -1;
    sr->hasEnd = PetscNot(pep->inta <= PETSC_MIN_REAL);
  }

  /* compute inertia0 */
  PetscCall(PEPQSliceGetInertia(pep,sr->int0,&sr->inertia0,ctx->detect?&zeros:NULL,ctx->hyperbolic?0:1));
  PetscCheck(!zeros || (sr->int0!=pep->inta && sr->int0!=pep->intb),((PetscObject)pep)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in the interval endpoint");
  if (!ctx->hyperbolic && ctx->checket) PetscCall(PEPQSliceCheckEigenvalueType(pep,sr->int0,0.0,PETSC_TRUE));

  /* compute inertia1 */
  PetscCall(PEPQSliceGetInertia(pep,sr->int1,&sr->inertia1,ctx->detect?&zeros:NULL,ctx->hyperbolic?0:1));
  PetscCheck(!zeros,((PetscObject)pep)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in an interval endpoint defined by user");
  if (!ctx->hyperbolic && ctx->checket && sr->hasEnd) {
    PetscCall(PEPQSliceCheckEigenvalueType(pep,sr->int1,0.0,PETSC_TRUE));
    PetscCheck(sr->type || sr->inertia1==sr->inertia0,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"No information of eigenvalue type in interval");
    PetscCheck(!sr->type || sr->inertia1!=sr->inertia0,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected");
    if (sr->dir*(sr->inertia1-sr->inertia0)<0) {
      sr->intcorr = -1;
      sr->inertia0 = 2*pep->n-sr->inertia0;
      sr->inertia1 = 2*pep->n-sr->inertia1;
    } else sr->intcorr = 1;
  } else {
    if (sr->inertia0<=pep->n && sr->inertia1<=pep->n) sr->intcorr = 1;
    else if (sr->inertia0>=pep->n && sr->inertia1>=pep->n) sr->intcorr = -1;
  }

  if (sr->hasEnd) {
    sr->dir = -sr->dir; r = sr->int0; sr->int0 = sr->int1; sr->int1 = r;
    i = sr->inertia0; sr->inertia0 = sr->inertia1; sr->inertia1 = i;
  }

  /* number of eigenvalues in interval */
  sr->numEigs = (sr->dir)*(sr->inertia1 - sr->inertia0);
  PetscCall(PetscInfo(pep,"QSlice setup: allocating for %" PetscInt_FMT " eigenvalues in [%g,%g]\n",sr->numEigs,(double)pep->inta,(double)pep->intb));
  if (sr->numEigs) {
    PetscCall(PEPQSliceAllocateSolution(pep));
    PetscCall(PEPSetDimensions_Default(pep,ctx->nev,&ctx->ncv,&ctx->mpd));
    pep->nev = ctx->nev; pep->ncv = ctx->ncv; pep->mpd = ctx->mpd;
    ld   = ctx->ncv+2;
    PetscCall(DSSetType(pep->ds,DSGHIEP));
    PetscCall(DSSetCompact(pep->ds,PETSC_TRUE));
    PetscCall(DSSetExtraRow(pep->ds,PETSC_TRUE));
    PetscCall(DSAllocate(pep->ds,ld));
    PetscCall(DSGetSlepcSC(pep->ds,&sc));
    sc->rg            = NULL;
    sc->comparison    = SlepcCompareLargestMagnitude;
    sc->comparisonctx = NULL;
    sc->map           = NULL;
    sc->mapobj        = NULL;
  } else {pep->ncv = 0; pep->nev = 0; pep->mpd = 0;}
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Fills the fields of a shift structure
*/
static PetscErrorCode PEPCreateShift(PEP pep,PetscReal val,PEP_shift neighb0,PEP_shift neighb1)
{
  PEP_shift      s,*pending2;
  PetscInt       i;
  PEP_SR         sr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  PetscCall(PetscNew(&s));
  s->value = val;
  s->neighb[0] = neighb0;
  if (neighb0) neighb0->neighb[1] = s;
  s->neighb[1] = neighb1;
  if (neighb1) neighb1->neighb[0] = s;
  s->comp[0] = PETSC_FALSE;
  s->comp[1] = PETSC_FALSE;
  s->index = -1;
  s->neigs = 0;
  s->nconv[0] = s->nconv[1] = 0;
  s->nsch[0] = s->nsch[1]=0;
  /* Inserts in the stack of pending shifts */
  /* If needed, the array is resized */
  if (sr->nPend >= sr->maxPend) {
    sr->maxPend *= 2;
    PetscCall(PetscMalloc1(sr->maxPend,&pending2));
    for (i=0;i<sr->nPend;i++) pending2[i] = sr->pending[i];
    PetscCall(PetscFree(sr->pending));
    sr->pending = pending2;
  }
  sr->pending[sr->nPend++]=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Provides next shift to be computed */
static PetscErrorCode PEPExtractShift(PEP pep)
{
  PetscInt       iner,zeros=0;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr;
  PetscReal      newShift,aux;
  PEP_shift      sPres;

  PetscFunctionBegin;
  sr = ctx->sr;
  if (sr->nPend > 0) {
    if (sr->dirch) {
      aux = sr->int1; sr->int1 = sr->int0; sr->int0 = aux;
      iner = sr->inertia1; sr->inertia1 = sr->inertia0; sr->inertia0 = iner;
      sr->dir *= -1;
      PetscCall(PetscFree(sr->s0->neighb[1]));
      PetscCall(PetscFree(sr->s0));
      sr->nPend--;
      PetscCall(PEPCreateShift(pep,sr->int0,NULL,NULL));
      sr->sPrev = NULL;
      sr->sPres = sr->pending[--sr->nPend];
      pep->target = sr->sPres->value;
      sr->s0 = sr->sPres;
      pep->reason = PEP_CONVERGED_ITERATING;
    } else {
      sr->sPrev = sr->sPres;
      sr->sPres = sr->pending[--sr->nPend];
    }
    sPres = sr->sPres;
    PetscCall(PEPQSliceGetInertia(pep,sPres->value,&iner,ctx->detect?&zeros:NULL,sr->intcorr));
    if (zeros) {
      newShift = sPres->value*(1.0+SLICE_PTOL);
      if (sr->dir*(sPres->neighb[0] && newShift-sPres->neighb[0]->value) < 0) newShift = (sPres->value+sPres->neighb[0]->value)/2;
      else if (sPres->neighb[1] && sr->dir*(sPres->neighb[1]->value-newShift) < 0) newShift = (sPres->value+sPres->neighb[1]->value)/2;
      PetscCall(PEPQSliceGetInertia(pep,newShift,&iner,&zeros,sr->intcorr));
      PetscCheck(!zeros,((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",(double)newShift);
      sPres->value = newShift;
    }
    sr->sPres->inertia = iner;
    pep->target = sr->sPres->value;
    pep->reason = PEP_CONVERGED_ITERATING;
    pep->its = 0;
  } else sr->sPres = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Obtains value of subsequent shift
*/
static PetscErrorCode PEPGetNewShiftValue(PEP pep,PetscInt side,PetscReal *newS)
{
  PetscReal lambda,d_prev;
  PetscInt  i,idxP;
  PEP_SR    sr;
  PEP_shift sPres,s;
  PEP_STOAR *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;
  if (sPres->neighb[side]) {
  /* Completing a previous interval */
    if (!sPres->neighb[side]->neighb[side] && sPres->neighb[side]->nconv[side]==0) { /* One of the ends might be too far from eigenvalues */
      if (side) *newS = (sPres->value + PetscRealPart(sr->eigr[sr->perm[sr->indexEig-1]]))/2;
      else *newS = (sPres->value + PetscRealPart(sr->eigr[sr->perm[0]]))/2;
    } else *newS=(sPres->value + sPres->neighb[side]->value)/2;
  } else { /* (Only for side=1). Creating a new interval. */
    if (sPres->neigs==0) {/* No value has been accepted*/
      if (sPres->neighb[0]) {
        /* Multiplying by 10 the previous distance */
        *newS = sPres->value + 10*sr->dir*PetscAbsReal(sPres->value - sPres->neighb[0]->value);
        sr->nleap++;
        /* Stops when the interval is open and no values are found in the last 5 shifts (there might be infinite eigenvalues) */
        PetscCheck(sr->hasEnd || sr->nleap<=5,PetscObjectComm((PetscObject)pep),PETSC_ERR_CONV_FAILED,"Unable to compute the wanted eigenvalues with open interval");
      } else { /* First shift */
        if (pep->nconv != 0) {
          /* Unaccepted values give information for next shift */
          idxP=0;/* Number of values left from shift */
          for (i=0;i<pep->nconv;i++) {
            lambda = PetscRealPart(pep->eigr[i]);
            if (sr->dir*(lambda - sPres->value) <0) idxP++;
            else break;
          }
          /* Avoiding subtraction of eigenvalues (might be the same).*/
          if (idxP>0) {
            d_prev = PetscAbsReal(sPres->value - PetscRealPart(pep->eigr[0]))/(idxP+0.3);
          } else {
            d_prev = PetscAbsReal(sPres->value - PetscRealPart(pep->eigr[pep->nconv-1]))/(pep->nconv+0.3);
          }
          *newS = sPres->value + (sr->dir*d_prev*pep->nev)/2;
          sr->dirch = PETSC_FALSE;
        } else { /* No values found, no information for next shift */
          PetscCheck(!sr->dirch,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"First shift renders no information");
          sr->dirch = PETSC_TRUE;
          *newS = sr->int1;
        }
      }
    } else { /* Accepted values found */
      sr->dirch = PETSC_FALSE;
      sr->nleap = 0;
      /* Average distance of values in previous subinterval */
      s = sPres->neighb[0];
      while (s && PetscAbs(s->inertia - sPres->inertia)==0) {
        s = s->neighb[0];/* Looking for previous shifts with eigenvalues within */
      }
      if (s) {
        d_prev = PetscAbsReal((sPres->value - s->value)/(sPres->inertia - s->inertia));
      } else { /* First shift. Average distance obtained with values in this shift */
        /* first shift might be too far from first wanted eigenvalue (no values found outside the interval)*/
        if (sr->dir*(PetscRealPart(sr->eigr[0])-sPres->value)>0 && PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0]))/PetscRealPart(sr->eigr[0])) > PetscSqrtReal(pep->tol)) {
          d_prev =  PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0])))/(sPres->neigs+0.3);
        } else {
          d_prev = PetscAbsReal(PetscRealPart(sr->eigr[sr->indexEig-1]) - sPres->value)/(sPres->neigs+0.3);
        }
      }
      /* Average distance is used for next shift by adding it to value on the right or to shift */
      if (sr->dir*(PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1]) - sPres->value)>0) {
        *newS = PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1])+ (sr->dir*d_prev*pep->nev)/2;
      } else { /* Last accepted value is on the left of shift. Adding to shift */
        *newS = sPres->value + (sr->dir*d_prev*pep->nev)/2;
      }
    }
    /* End of interval can not be surpassed */
    if (sr->dir*(sr->int1 - *newS) < 0) *newS = sr->int1;
  }/* of neighb[side]==null */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Function for sorting an array of real values
*/
static PetscErrorCode sortRealEigenvalues(PetscScalar *r,PetscInt *perm,PetscInt nr,PetscBool prev,PetscInt dir)
{
  PetscReal re;
  PetscInt  i,j,tmp;

  PetscFunctionBegin;
  if (!prev) for (i=0;i<nr;i++) perm[i] = i;
  /* Insertion sort */
  for (i=1;i<nr;i++) {
    re = PetscRealPart(r[perm[i]]);
    j = i-1;
    while (j>=0 && dir*(re - PetscRealPart(r[perm[j]])) <= 0) {
      tmp = perm[j]; perm[j] = perm[j+1]; perm[j+1] = tmp; j--;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Stores the pairs obtained since the last shift in the global arrays */
static PetscErrorCode PEPStoreEigenpairs(PEP pep)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PetscReal      lambda,err,*errest;
  PetscInt       i,*aux,count=0,ndef,ld,nconv=pep->nconv,d=pep->nmat-1,idx;
  PetscBool      iscayley,divide=PETSC_FALSE;
  PEP_SR         sr = ctx->sr;
  PEP_shift      sPres;
  Vec            w,vomega;
  Mat            MS;
  BV             tV;
  PetscScalar    *S,*eigr,*tS,*omega;

  PetscFunctionBegin;
  sPres = sr->sPres;
  sPres->index = sr->indexEig;

  if (nconv>sr->ndef0+sr->ndef1) {
    /* Back-transform */
    PetscCall(STBackTransform(pep->st,nconv,pep->eigr,pep->eigi));
    for (i=0;i<nconv;i++) {
#if defined(PETSC_USE_COMPLEX)
      if (PetscImaginaryPart(pep->eigr[i])) pep->eigr[i] = sr->int0-sr->dir;
#else
      if (pep->eigi[i]) pep->eigr[i] = sr->int0-sr->dir;
#endif
    }
    PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STCAYLEY,&iscayley));
    /* Sort eigenvalues */
    PetscCall(sortRealEigenvalues(pep->eigr,pep->perm,nconv,PETSC_FALSE,sr->dir));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,nconv,&vomega));
    PetscCall(BVGetSignature(ctx->V,vomega));
    PetscCall(VecGetArray(vomega,&omega));
    PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArray(MS,&S));
    /* Values stored in global array */
    PetscCall(PetscCalloc4(nconv,&eigr,nconv,&errest,nconv*nconv*d,&tS,nconv,&aux));
    ndef = sr->ndef0+sr->ndef1;
    for (i=0;i<nconv;i++) {
      lambda = PetscRealPart(pep->eigr[pep->perm[i]]);
      err = pep->errest[pep->perm[i]];
      if (sr->dir*(lambda - sPres->ext[0]) > 0 && (sr->dir)*(sPres->ext[1] - lambda) > 0) {/* Valid value */
        PetscCheck(sr->indexEig+count-ndef<sr->numEigs,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Unexpected error in Spectrum Slicing");
        PetscCall(PEPQSliceCheckEigenvalueType(pep,lambda,PetscRealPart(omega[pep->perm[i]]),PETSC_FALSE));
        eigr[count] = lambda;
        errest[count] = err;
        if ((sr->dir*(sPres->value - lambda) > 0) && (sr->dir*(lambda - sPres->ext[0]) > 0)) sPres->nconv[0]++;
        if ((sr->dir*(lambda - sPres->value) > 0) && (sr->dir*(sPres->ext[1] - lambda) > 0)) sPres->nconv[1]++;
        PetscCall(PetscArraycpy(tS+count*(d*nconv),S+pep->perm[i]*(d*ld),nconv));
        PetscCall(PetscArraycpy(tS+count*(d*nconv)+nconv,S+pep->perm[i]*(d*ld)+ld,nconv));
        count++;
      }
    }
    PetscCall(VecRestoreArray(vomega,&omega));
    PetscCall(VecDestroy(&vomega));
    for (i=0;i<count;i++) {
      PetscCall(PetscArraycpy(S+i*(d*ld),tS+i*nconv*d,nconv));
      PetscCall(PetscArraycpy(S+i*(d*ld)+ld,tS+i*nconv*d+nconv,nconv));
    }
    PetscCall(MatDenseRestoreArray(MS,&S));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    PetscCall(BVSetActiveColumns(ctx->V,0,count));
    PetscCall(BVTensorCompress(ctx->V,count));
    if (sr->sPres->nconv[0] && sr->sPres->nconv[1]) {
      divide = PETSC_TRUE;
      PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
      PetscCall(MatDenseGetArray(MS,&S));
      PetscCall(PetscArrayzero(tS,nconv*nconv*d));
      for (i=0;i<count;i++) {
        PetscCall(PetscArraycpy(tS+i*nconv*d,S+i*(d*ld),count));
        PetscCall(PetscArraycpy(tS+i*nconv*d+nconv,S+i*(d*ld)+ld,count));
      }
      PetscCall(MatDenseRestoreArray(MS,&S));
      PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
      PetscCall(BVSetActiveColumns(pep->V,0,count));
      PetscCall(BVDuplicateResize(pep->V,count,&tV));
      PetscCall(BVCopy(pep->V,tV));
    }
    if (sr->sPres->nconv[0]) {
      if (divide) {
        PetscCall(BVSetActiveColumns(ctx->V,0,sr->sPres->nconv[0]));
        PetscCall(BVTensorCompress(ctx->V,sr->sPres->nconv[0]));
      }
      for (i=0;i<sr->ndef0;i++) aux[i] = sr->idxDef0[i];
      for (i=sr->ndef0;i<sr->sPres->nconv[0];i++) aux[i] = sr->indexEig+i-sr->ndef0;
      PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
      PetscCall(MatDenseGetArray(MS,&S));
      for (i=0;i<sr->sPres->nconv[0];i++) {
        sr->eigr[aux[i]] = eigr[i];
        sr->errest[aux[i]] = errest[i];
        PetscCall(BVGetColumn(pep->V,i,&w));
        PetscCall(BVInsertVec(sr->V,aux[i],w));
        PetscCall(BVRestoreColumn(pep->V,i,&w));
        idx = sr->ld*d*aux[i];
        PetscCall(PetscArrayzero(sr->S+idx,sr->ld*d));
        PetscCall(PetscArraycpy(sr->S+idx,S+i*(ld*d),sr->sPres->nconv[0]));
        PetscCall(PetscArraycpy(sr->S+idx+sr->ld,S+i*(ld*d)+ld,sr->sPres->nconv[0]));
        PetscCall(PetscFree(sr->qinfo[aux[i]].q));
        PetscCall(PetscMalloc1(sr->sPres->nconv[0],&sr->qinfo[aux[i]].q));
        PetscCall(PetscArraycpy(sr->qinfo[aux[i]].q,aux,sr->sPres->nconv[0]));
        sr->qinfo[aux[i]].nq = sr->sPres->nconv[0];
      }
      PetscCall(MatDenseRestoreArray(MS,&S));
      PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    }

    if (sr->sPres->nconv[1]) {
      if (divide) {
        PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
        PetscCall(MatDenseGetArray(MS,&S));
        for (i=0;i<sr->sPres->nconv[1];i++) {
          PetscCall(PetscArraycpy(S+i*(d*ld),tS+(sr->sPres->nconv[0]+i)*nconv*d,count));
          PetscCall(PetscArraycpy(S+i*(d*ld)+ld,tS+(sr->sPres->nconv[0]+i)*nconv*d+nconv,count));
        }
        PetscCall(MatDenseRestoreArray(MS,&S));
        PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
        PetscCall(BVSetActiveColumns(pep->V,0,count));
        PetscCall(BVCopy(tV,pep->V));
        PetscCall(BVSetActiveColumns(ctx->V,0,sr->sPres->nconv[1]));
        PetscCall(BVTensorCompress(ctx->V,sr->sPres->nconv[1]));
      }
      for (i=0;i<sr->ndef1;i++) aux[i] = sr->idxDef1[i];
      for (i=sr->ndef1;i<sr->sPres->nconv[1];i++) aux[i] = sr->indexEig+sr->sPres->nconv[0]-sr->ndef0+i-sr->ndef1;
      PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
      PetscCall(MatDenseGetArray(MS,&S));
      for (i=0;i<sr->sPres->nconv[1];i++) {
        sr->eigr[aux[i]] = eigr[sr->sPres->nconv[0]+i];
        sr->errest[aux[i]] = errest[sr->sPres->nconv[0]+i];
        PetscCall(BVGetColumn(pep->V,i,&w));
        PetscCall(BVInsertVec(sr->V,aux[i],w));
        PetscCall(BVRestoreColumn(pep->V,i,&w));
        idx = sr->ld*d*aux[i];
        PetscCall(PetscArrayzero(sr->S+idx,sr->ld*d));
        PetscCall(PetscArraycpy(sr->S+idx,S+i*(ld*d),sr->sPres->nconv[1]));
        PetscCall(PetscArraycpy(sr->S+idx+sr->ld,S+i*(ld*d)+ld,sr->sPres->nconv[1]));
        PetscCall(PetscFree(sr->qinfo[aux[i]].q));
        PetscCall(PetscMalloc1(sr->sPres->nconv[1],&sr->qinfo[aux[i]].q));
        PetscCall(PetscArraycpy(sr->qinfo[aux[i]].q,aux,sr->sPres->nconv[1]));
        sr->qinfo[aux[i]].nq = sr->sPres->nconv[1];
      }
      PetscCall(MatDenseRestoreArray(MS,&S));
      PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    }
    sPres->neigs = count-sr->ndef0-sr->ndef1;
    sr->indexEig += sPres->neigs;
    sPres->nconv[0]-= sr->ndef0;
    sPres->nconv[1]-= sr->ndef1;
    PetscCall(PetscFree4(eigr,errest,tS,aux));
  } else {
    sPres->neigs = 0;
    sPres->nconv[0]= 0;
    sPres->nconv[1]= 0;
  }
  /* Global ordering array updating */
  PetscCall(sortRealEigenvalues(sr->eigr,sr->perm,sr->indexEig,PETSC_FALSE,sr->dir));
  /* Check for completion */
  sPres->comp[0] = PetscNot(sPres->nconv[0] < sPres->nsch[0]);
  sPres->comp[1] = PetscNot(sPres->nconv[1] < sPres->nsch[1]);
  PetscCheck(sPres->nconv[0]<=sPres->nsch[0] && sPres->nconv[1]<=sPres->nsch[1],PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia");
  if (divide) PetscCall(BVDestroy(&tV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPLookForDeflation(PEP pep)
{
  PetscReal val;
  PetscInt  i,count0=0,count1=0;
  PEP_shift sPres;
  PetscInt  ini,fin;
  PEP_SR    sr;
  PEP_STOAR *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;

  if (sPres->neighb[0]) ini = (sr->dir)*(sPres->neighb[0]->inertia - sr->inertia0);
  else ini = 0;
  fin = sr->indexEig;
  /* Selection of ends for searching new values */
  if (!sPres->neighb[0]) sPres->ext[0] = sr->int0;/* First shift */
  else sPres->ext[0] = sPres->neighb[0]->value;
  if (!sPres->neighb[1]) {
    if (sr->hasEnd) sPres->ext[1] = sr->int1;
    else sPres->ext[1] = (sr->dir > 0)?PETSC_MAX_REAL:PETSC_MIN_REAL;
  } else sPres->ext[1] = sPres->neighb[1]->value;
  /* Selection of values between right and left ends */
  for (i=ini;i<fin;i++) {
    val=PetscRealPart(sr->eigr[sr->perm[i]]);
    /* Values to the right of left shift */
    if (sr->dir*(val - sPres->ext[1]) < 0) {
      if (sr->dir*(val - sPres->value) < 0) count0++;
      else count1++;
    } else break;
  }
  /* The number of values on each side are found */
  if (sPres->neighb[0]) {
    sPres->nsch[0] = (sr->dir)*(sPres->inertia - sPres->neighb[0]->inertia)-count0;
    PetscCheck(sPres->nsch[0]>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia");
  } else sPres->nsch[0] = 0;

  if (sPres->neighb[1]) {
    sPres->nsch[1] = (sr->dir)*(sPres->neighb[1]->inertia - sPres->inertia) - count1;
    PetscCheck(sPres->nsch[1]>=0,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Mismatch between number of values found and information from inertia");
  } else sPres->nsch[1] = (sr->dir)*(sr->inertia1 - sPres->inertia);

  /* Completing vector of indexes for deflation */
  for (i=0;i<count0;i++) sr->idxDef0[i] = sr->perm[ini+i];
  sr->ndef0 = count0;
  for (i=0;i<count1;i++) sr->idxDef1[i] = sr->perm[ini+count0+i];
  sr->ndef1 = count1;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute a run of Lanczos iterations
*/
static PetscErrorCode PEPSTOARrun_QSlice(PEP pep,PetscReal *a,PetscReal *b,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscBool *symmlost,Vec *t_)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       i,j,m=*M,l,lock;
  PetscInt       lds,d,ld,offq,nqt,ldds;
  Vec            v=t_[0],t=t_[1],q=t_[2];
  PetscReal      norm,sym=0.0,fro=0.0,*f;
  PetscScalar    *y,*S,sigma;
  PetscBLASInt   j_,one=1;
  PetscBool      lindep;
  Mat            MS;

  PetscFunctionBegin;
  PetscCall(PetscMalloc1(*M,&y));
  PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
  PetscCall(BVTensorGetDegree(ctx->V,&d));
  PetscCall(BVGetActiveColumns(pep->V,&lock,&nqt));
  lds = d*ld;
  offq = ld;
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));

  *breakdown = PETSC_FALSE; /* ----- */
  PetscCall(STGetShift(pep->st,&sigma));
  PetscCall(DSGetDimensions(pep->ds,NULL,&l,NULL,NULL));
  PetscCall(BVSetActiveColumns(ctx->V,0,m));
  PetscCall(BVSetActiveColumns(pep->V,0,nqt));
  for (j=k;j<m;j++) {
    /* apply operator */
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArray(MS,&S));
    PetscCall(BVGetColumn(pep->V,nqt,&t));
    PetscCall(BVMultVec(pep->V,1.0,0.0,v,S+j*lds));
    PetscCall(MatMult(pep->A[1],v,q));
    PetscCall(MatMult(pep->A[2],v,t));
    PetscCall(VecAXPY(q,sigma*pep->sfactor,t));
    PetscCall(VecScale(q,pep->sfactor));
    PetscCall(BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds));
    PetscCall(MatMult(pep->A[2],v,t));
    PetscCall(VecAXPY(q,pep->sfactor*pep->sfactor,t));
    PetscCall(STMatSolve(pep->st,q,t));
    PetscCall(VecScale(t,-1.0));
    PetscCall(BVRestoreColumn(pep->V,nqt,&t));

    /* orthogonalize */
    PetscCall(BVOrthogonalizeColumn(pep->V,nqt,S+(j+1)*lds,&norm,&lindep));
    if (!lindep) {
      *(S+(j+1)*lds+nqt) = norm;
      PetscCall(BVScaleColumn(pep->V,nqt,1.0/norm));
      nqt++;
    }
    for (i=0;i<nqt;i++) *(S+(j+1)*lds+offq+i) = *(S+j*lds+i)+sigma*(*(S+(j+1)*lds+i));
    PetscCall(BVSetActiveColumns(pep->V,0,nqt));
    PetscCall(MatDenseRestoreArray(MS,&S));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));

    /* level-2 orthogonalization */
    PetscCall(BVOrthogonalizeColumn(ctx->V,j+1,y,&norm,&lindep));
    a[j] = PetscRealPart(y[j]);
    omega[j+1] = (norm > 0)?1.0:-1.0;
    PetscCall(BVScaleColumn(ctx->V,j+1,1.0/norm));
    b[j] = PetscAbsReal(norm);

    /* check symmetry */
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_T,&f));
    if (j==k) {
      for (i=l;i<j-1;i++) y[i] = PetscAbsScalar(y[i])-PetscAbsReal(f[2*ldds+i]);
      for (i=0;i<l;i++) y[i] = 0.0;
    }
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_T,&f));
    if (j>0) y[j-1] = PetscAbsScalar(y[j-1])-PetscAbsReal(b[j-1]);
    PetscCall(PetscBLASIntCast(j,&j_));
    sym = SlepcAbs(BLASnrm2_(&j_,y,&one),sym);
    fro = SlepcAbs(fro,SlepcAbs(a[j],b[j]));
    if (j>0) fro = SlepcAbs(fro,b[j-1]);
    if (sym/fro>PetscMax(PETSC_SQRT_MACHINE_EPSILON,10*pep->tol)) {
      *symmlost = PETSC_TRUE;
      *M=j;
      break;
    }
  }
  PetscCall(BVSetActiveColumns(pep->V,lock,nqt));
  PetscCall(BVSetActiveColumns(ctx->V,0,*M));
  PetscCall(PetscFree(y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPSTOAR_QSlice(PEP pep,Mat B)
{
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       j,k,l,nv=0,ld,ldds,t,nq=0,idx;
  PetscInt       nconv=0,deg=pep->nmat-1,count0=0,count1=0;
  PetscScalar    *om,sigma,*back,*S,*pQ;
  PetscReal      beta,norm=1.0,*omega,*a,*b,eta,lambda;
  PetscBool      breakdown,symmlost=PETSC_FALSE,sinv,falselock=PETSC_TRUE;
  Mat            MS,MQ,D;
  Vec            v,vomega;
  PEP_SR         sr;
  BVOrthogType   otype;
  BVOrthogBlockType obtype;

  PetscFunctionBegin;
  /* Resize if needed for deflating vectors  */
  sr = ctx->sr;
  sigma = sr->sPres->value;
  k = sr->ndef0+sr->ndef1;
  pep->ncv = ctx->ncv+k;
  pep->nev = ctx->nev+k;
  PetscCall(PEPAllocateSolution(pep,3));
  PetscCall(BVDestroy(&ctx->V));
  PetscCall(BVCreateTensor(pep->V,pep->nmat-1,&ctx->V));
  PetscCall(BVGetOrthogonalization(pep->V,&otype,NULL,&eta,&obtype));
  PetscCall(BVSetOrthogonalization(ctx->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype));
  PetscCall(DSAllocate(pep->ds,pep->ncv+2));
  PetscCall(PetscMalloc1(pep->ncv,&back));
  PetscCall(DSGetLeadingDimension(pep->ds,&ldds));
  PetscCall(BVSetMatrix(ctx->V,B,PETSC_TRUE));
  PetscCheck(ctx->lock,PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"A locking variant is needed for spectrum slicing");
  /* undocumented option to use a cheaper locking instead of the true locking */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-pep_stoar_falselocking",&falselock,NULL));
  PetscCall(PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv));
  PetscCall(RGPushScale(pep->rg,sinv?pep->sfactor:1.0/pep->sfactor));
  PetscCall(STScaleShift(pep->st,sinv?pep->sfactor:1.0/pep->sfactor));

  /* Get the starting Arnoldi vector */
  PetscCall(BVSetActiveColumns(pep->V,0,1));
  PetscCall(BVTensorBuildFirstColumn(ctx->V,pep->nini));
  PetscCall(BVSetActiveColumns(ctx->V,0,1));
  if (k) {
    /* Insert deflated vectors */
    PetscCall(BVSetActiveColumns(pep->V,0,0));
    idx = sr->ndef0?sr->idxDef0[0]:sr->idxDef1[0];
    for (j=0;j<k;j++) {
      PetscCall(BVGetColumn(pep->V,j,&v));
      PetscCall(BVCopyVec(sr->V,sr->qinfo[idx].q[j],v));
      PetscCall(BVRestoreColumn(pep->V,j,&v));
    }
    /* Update innerproduct matrix */
    PetscCall(BVSetActiveColumns(ctx->V,0,0));
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(BVSetActiveColumns(pep->V,0,k));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));

    PetscCall(BVGetSizes(pep->V,NULL,NULL,&ld));
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArray(MS,&S));
    for (j=0;j<sr->ndef0;j++) {
      PetscCall(PetscArrayzero(S+j*ld*deg,ld*deg));
      PetscCall(PetscArraycpy(S+j*ld*deg,sr->S+sr->idxDef0[j]*sr->ld*deg,k));
      PetscCall(PetscArraycpy(S+j*ld*deg+ld,sr->S+sr->idxDef0[j]*sr->ld*deg+sr->ld,k));
      pep->eigr[j] = sr->eigr[sr->idxDef0[j]];
      pep->errest[j] = sr->errest[sr->idxDef0[j]];
    }
    for (j=0;j<sr->ndef1;j++) {
      PetscCall(PetscArrayzero(S+(j+sr->ndef0)*ld*deg,ld*deg));
      PetscCall(PetscArraycpy(S+(j+sr->ndef0)*ld*deg,sr->S+sr->idxDef1[j]*sr->ld*deg,k));
      PetscCall(PetscArraycpy(S+(j+sr->ndef0)*ld*deg+ld,sr->S+sr->idxDef1[j]*sr->ld*deg+sr->ld,k));
      pep->eigr[j+sr->ndef0] = sr->eigr[sr->idxDef1[j]];
      pep->errest[j+sr->ndef0] = sr->errest[sr->idxDef1[j]];
    }
    PetscCall(MatDenseRestoreArray(MS,&S));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    PetscCall(BVSetActiveColumns(ctx->V,0,k+1));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,k+1,&vomega));
    PetscCall(VecGetArray(vomega,&om));
    for (j=0;j<k;j++) {
      PetscCall(BVOrthogonalizeColumn(ctx->V,j,NULL,&norm,NULL));
      PetscCall(BVScaleColumn(ctx->V,j,1/norm));
      om[j] = (norm>=0.0)?1.0:-1.0;
    }
    PetscCall(BVTensorGetFactors(ctx->V,NULL,&MS));
    PetscCall(MatDenseGetArray(MS,&S));
    for (j=0;j<deg;j++) {
      PetscCall(BVSetRandomColumn(pep->V,k+j));
      PetscCall(BVOrthogonalizeColumn(pep->V,k+j,S+k*ld*deg+j*ld,&norm,NULL));
      PetscCall(BVScaleColumn(pep->V,k+j,1.0/norm));
      S[k*ld*deg+j*ld+k+j] = norm;
    }
    PetscCall(MatDenseRestoreArray(MS,&S));
    PetscCall(BVSetActiveColumns(pep->V,0,k+deg));
    PetscCall(BVTensorRestoreFactors(ctx->V,NULL,&MS));
    PetscCall(BVOrthogonalizeColumn(ctx->V,k,NULL,&norm,NULL));
    PetscCall(BVScaleColumn(ctx->V,k,1.0/norm));
    om[k] = (norm>=0.0)?1.0:-1.0;
    PetscCall(VecRestoreArray(vomega,&om));
    PetscCall(BVSetSignature(ctx->V,vomega));
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_T,&a));
    PetscCall(VecGetArray(vomega,&om));
    for (j=0;j<k;j++) a[j] = PetscRealPart(om[j]/(pep->eigr[j]-sigma));
    PetscCall(VecRestoreArray(vomega,&om));
    PetscCall(VecDestroy(&vomega));
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_T,&a));
    PetscCall(DSGetArray(pep->ds,DS_MAT_Q,&pQ));
    PetscCall(PetscArrayzero(pQ,ldds*k));
    for (j=0;j<k;j++) pQ[j+j*ldds] = 1.0;
    PetscCall(DSRestoreArray(pep->ds,DS_MAT_Q,&pQ));
  }
  PetscCall(BVSetActiveColumns(ctx->V,0,k+1));
  PetscCall(DSSetDimensions(pep->ds,k+1,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(DSGetMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));
  PetscCall(BVGetSignature(ctx->V,vomega));
  PetscCall(DSRestoreMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));

  PetscCall(PetscInfo(pep,"Start STOAR: sigma=%g in [%g,%g], for deflation: left=%" PetscInt_FMT " right=%" PetscInt_FMT ", searching: left=%" PetscInt_FMT " right=%" PetscInt_FMT "\n",(double)sr->sPres->value,(double)(sr->sPres->neighb[0]?sr->sPres->neighb[0]->value:sr->int0),(double)(sr->sPres->neighb[1]?sr->sPres->neighb[1]->value:sr->int1),sr->ndef0,sr->ndef1,sr->sPres->nsch[0],sr->sPres->nsch[1]));

  /* Restart loop */
  l = 0;
  pep->nconv = k;
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_T,&a));
    b = a+ldds;
    PetscCall(DSGetArrayReal(pep->ds,DS_MAT_D,&omega));

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(pep->nconv+pep->mpd,pep->ncv);
    PetscCall(PEPSTOARrun_QSlice(pep,a,b,omega,pep->nconv+l,&nv,&breakdown,&symmlost,pep->work));
    beta = b[nv-1];
    if (symmlost && nv==pep->nconv+l) {
      pep->reason = PEP_DIVERGED_SYMMETRY_LOST;
      pep->nconv = nconv;
      PetscCall(PetscInfo(pep,"Symmetry lost in STOAR sigma=%g nconv=%" PetscInt_FMT "\n",(double)sr->sPres->value,nconv));
      if (falselock || !ctx->lock) {
        PetscCall(BVSetActiveColumns(ctx->V,0,pep->nconv));
        PetscCall(BVTensorCompress(ctx->V,0));
      }
      break;
    }
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_T,&a));
    PetscCall(DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega));
    PetscCall(DSSetDimensions(pep->ds,nv,pep->nconv,pep->nconv+l));
    if (l==0) PetscCall(DSSetState(pep->ds,DS_STATE_INTERMEDIATE));
    else PetscCall(DSSetState(pep->ds,DS_STATE_RAW));

    /* Solve projected problem */
    PetscCall(DSSolve(pep->ds,pep->eigr,pep->eigi));
    PetscCall(DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL));
    PetscCall(DSUpdateExtraRow(pep->ds));
    PetscCall(DSSynchronize(pep->ds,pep->eigr,pep->eigi));

    /* Check convergence */
    /* PetscCall(PEPSTOARpreKConvergence(pep,nv,&norm,pep->work));*/
    norm = 1.0;
    PetscCall(DSGetDimensions(pep->ds,NULL,NULL,NULL,&t));
    PetscCall(PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,t-pep->nconv,PetscAbsReal(beta)*norm,&k));
    PetscCall((*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx));
    for (j=0;j<k;j++) back[j] = pep->eigr[j];
    PetscCall(STBackTransform(pep->st,k,back,pep->eigi));
    count0=count1=0;
    for (j=0;j<k;j++) {
      lambda = PetscRealPart(back[j]);
      if ((sr->dir*(sr->sPres->value - lambda) > 0) && (sr->dir*(lambda - sr->sPres->ext[0]) > 0)) count0++;
      if ((sr->dir*(lambda - sr->sPres->value) > 0) && (sr->dir*(sr->sPres->ext[1] - lambda) > 0)) count1++;
    }
    if ((count0-sr->ndef0 >= sr->sPres->nsch[0]) && (count1-sr->ndef1 >= sr->sPres->nsch[1])) pep->reason = PEP_CONVERGED_TOL;
    /* Update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)/2));
      l = PetscMin(l,t);
      PetscCall(DSGetTruncateSize(pep->ds,k,t,&l));
      if (!breakdown) {
        /* Prepare the Rayleigh quotient for restart */
        PetscCall(DSTruncate(pep->ds,k+l,PETSC_FALSE));
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */
    if (l) PetscCall(PetscInfo(pep,"Preparing to restart keeping l=%" PetscInt_FMT " vectors\n",l));

    /* Update S */
    PetscCall(DSGetMat(pep->ds,DS_MAT_Q,&MQ));
    PetscCall(BVMultInPlace(ctx->V,MQ,pep->nconv,k+l));
    PetscCall(DSRestoreMat(pep->ds,DS_MAT_Q,&MQ));

    /* Copy last column of S */
    PetscCall(BVCopyColumn(ctx->V,nv,k+l));
    PetscCall(BVSetActiveColumns(ctx->V,0,k+l));
    if (k+l) {
      PetscCall(DSSetDimensions(pep->ds,k+l,PETSC_DETERMINE,PETSC_DETERMINE));
      PetscCall(DSGetMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));
      PetscCall(BVSetSignature(ctx->V,vomega));
      PetscCall(DSRestoreMatAndColumn(pep->ds,DS_MAT_D,0,&D,&vomega));
    }

    if (breakdown && pep->reason == PEP_CONVERGED_ITERATING) {
      /* stop if breakdown */
      PetscCall(PetscInfo(pep,"Breakdown TOAR method (it=%" PetscInt_FMT " norm=%g)\n",pep->its,(double)beta));
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
    PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
    if (k+l+deg<=nq) {
      PetscCall(BVSetActiveColumns(ctx->V,pep->nconv,k+l+1));
      if (!falselock && ctx->lock) PetscCall(BVTensorCompress(ctx->V,k-pep->nconv));
      else PetscCall(BVTensorCompress(ctx->V,0));
    }
    pep->nconv = k;
    PetscCall(PEPMonitor(pep,pep->its,nconv,pep->eigr,pep->eigi,pep->errest,nv));
  }
  sr->itsKs += pep->its;
  if (pep->nconv>0) {
    PetscCall(BVSetActiveColumns(ctx->V,0,pep->nconv));
    PetscCall(BVGetActiveColumns(pep->V,NULL,&nq));
    PetscCall(BVSetActiveColumns(pep->V,0,nq));
    if (nq>pep->nconv) {
      PetscCall(BVTensorCompress(ctx->V,pep->nconv));
      PetscCall(BVSetActiveColumns(pep->V,0,pep->nconv));
    }
    for (j=0;j<pep->nconv;j++) {
      pep->eigr[j] *= pep->sfactor;
      pep->eigi[j] *= pep->sfactor;
    }
  }
  PetscCall(PetscInfo(pep,"Finished STOAR: nconv=%" PetscInt_FMT " (deflated=%" PetscInt_FMT ", left=%" PetscInt_FMT ", right=%" PetscInt_FMT ")\n",pep->nconv,sr->ndef0+sr->ndef1,count0-sr->ndef0,count1-sr->ndef1));
  PetscCall(STScaleShift(pep->st,sinv?1.0/pep->sfactor:pep->sfactor));
  PetscCall(RGPopScale(pep->rg));

  PetscCheck(pep->reason!=PEP_DIVERGED_SYMMETRY_LOST || nconv>=sr->ndef0+sr->ndef1,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Symmetry lost at sigma=%g",(double)sr->sPres->value);
  if (pep->reason == PEP_DIVERGED_SYMMETRY_LOST && nconv==sr->ndef0+sr->ndef1) {
    PetscCheck(++sr->symmlost<=10,PetscObjectComm((PetscObject)pep),PETSC_ERR_PLIB,"Symmetry lost at sigma=%g",(double)sr->sPres->value);
  } else sr->symmlost = 0;

  PetscCall(DSTruncate(pep->ds,pep->nconv,PETSC_TRUE));
  PetscCall(PetscFree(back));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PEPQSliceGetInertias(PEP pep,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR          sr=ctx->sr;
  PetscInt        i=0,j,tmpi;
  PetscReal       v,tmpr;
  PEP_shift       s;

  PetscFunctionBegin;
  PetscCheck(pep->state,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSetUp() first");
  PetscCheck(sr,PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see PEPSetInterval()");
  if (!sr->s0) {  /* PEPSolve not called yet */
    *n = 2;
  } else {
    *n = 1;
    s = sr->s0;
    while (s) {
      (*n)++;
      s = s->neighb[1];
    }
  }
  PetscCall(PetscMalloc1(*n,shifts));
  PetscCall(PetscMalloc1(*n,inertias));
  if (!sr->s0) {  /* PEPSolve not called yet */
    (*shifts)[0]   = sr->int0;
    (*shifts)[1]   = sr->int1;
    (*inertias)[0] = sr->inertia0;
    (*inertias)[1] = sr->inertia1;
  } else {
    s = sr->s0;
    while (s) {
      (*shifts)[i]     = s->value;
      (*inertias)[i++] = s->inertia;
      s = s->neighb[1];
    }
    (*shifts)[i]   = sr->int1;
    (*inertias)[i] = sr->inertia1;
  }
  /* remove possible duplicate in last position */
  if ((*shifts)[(*n)-1]==(*shifts)[(*n)-2]) (*n)--;
  /* sort result */
  for (i=0;i<*n;i++) {
    v = (*shifts)[i];
    for (j=i+1;j<*n;j++) {
      if (v > (*shifts)[j]) {
        SlepcSwap((*shifts)[i],(*shifts)[j],tmpr);
        SlepcSwap((*inertias)[i],(*inertias)[j],tmpi);
        v = (*shifts)[i];
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PEPSolve_STOAR_QSlice(PEP pep)
{
  PetscInt       i,j,ti,deg=pep->nmat-1;
  PetscReal      newS;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;
  Mat            S,B;
  PetscScalar    *pS;

  PetscFunctionBegin;
  PetscCall(PetscCitationsRegister(citation,&cited));

  /* Only with eigenvalues present in the interval ...*/
  if (sr->numEigs==0) {
    pep->reason = PEP_CONVERGED_TOL;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* Inner product matrix */
  PetscCall(PEPSTOARSetUpInnerMatrix(pep,&B));

  /* Array of pending shifts */
  sr->maxPend = 100; /* Initial size */
  sr->nPend = 0;
  PetscCall(PetscMalloc1(sr->maxPend,&sr->pending));
  PetscCall(PEPCreateShift(pep,sr->int0,NULL,NULL));
  /* extract first shift */
  sr->sPrev = NULL;
  sr->sPres = sr->pending[--sr->nPend];
  sr->sPres->inertia = sr->inertia0;
  pep->target = sr->sPres->value;
  sr->s0 = sr->sPres;
  sr->indexEig = 0;

  for (i=0;i<sr->numEigs;i++) {
    sr->eigr[i]   = 0.0;
    sr->eigi[i]   = 0.0;
    sr->errest[i] = 0.0;
    sr->perm[i]   = i;
  }
  /* Vectors for deflation */
  PetscCall(PetscMalloc2(sr->numEigs,&sr->idxDef0,sr->numEigs,&sr->idxDef1));
  sr->indexEig = 0;
  while (sr->sPres) {
    /* Search for deflation */
    PetscCall(PEPLookForDeflation(pep));
    /* KrylovSchur */
    PetscCall(PEPSTOAR_QSlice(pep,B));

    PetscCall(PEPStoreEigenpairs(pep));
    /* Select new shift */
    if (!sr->sPres->comp[1]) {
      PetscCall(PEPGetNewShiftValue(pep,1,&newS));
      PetscCall(PEPCreateShift(pep,newS,sr->sPres,sr->sPres->neighb[1]));
    }
    if (!sr->sPres->comp[0]) {
      /* Completing earlier interval */
      PetscCall(PEPGetNewShiftValue(pep,0,&newS));
      PetscCall(PEPCreateShift(pep,newS,sr->sPres->neighb[0],sr->sPres));
    }
    /* Preparing for a new search of values */
    PetscCall(PEPExtractShift(pep));
  }

  /* Updating pep values prior to exit */
  PetscCall(PetscFree2(sr->idxDef0,sr->idxDef1));
  PetscCall(PetscFree(sr->pending));
  pep->nconv  = sr->indexEig;
  pep->reason = PEP_CONVERGED_TOL;
  pep->its    = sr->itsKs;
  pep->nev    = sr->indexEig;
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,pep->nconv,pep->nconv,NULL,&S));
  PetscCall(MatDenseGetArray(S,&pS));
  for (i=0;i<pep->nconv;i++) {
    for (j=0;j<sr->qinfo[i].nq;j++) pS[i*pep->nconv+sr->qinfo[i].q[j]] = *(sr->S+i*sr->ld*deg+j);
  }
  PetscCall(MatDenseRestoreArray(S,&pS));
  PetscCall(BVSetActiveColumns(sr->V,0,pep->nconv));
  PetscCall(BVMultInPlace(sr->V,S,0,pep->nconv));
  PetscCall(MatDestroy(&S));
  PetscCall(BVDestroy(&pep->V));
  pep->V = sr->V;
  PetscCall(PetscFree4(pep->eigr,pep->eigi,pep->errest,pep->perm));
  pep->eigr   = sr->eigr;
  pep->eigi   = sr->eigi;
  pep->perm   = sr->perm;
  pep->errest = sr->errest;
  if (sr->dir<0) {
    for (i=0;i<pep->nconv/2;i++) {
      ti = sr->perm[i]; sr->perm[i] = sr->perm[pep->nconv-1-i]; sr->perm[pep->nconv-1-i] = ti;
    }
  }
  PetscCall(PetscFree(ctx->inertias));
  PetscCall(PetscFree(ctx->shifts));
  PetscCall(MatDestroy(&B));
  PetscCall(PEPQSliceGetInertias(pep,&ctx->nshifts,&ctx->shifts,&ctx->inertias));
  PetscFunctionReturn(PETSC_SUCCESS);
}
