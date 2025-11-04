/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Filter spectral transformation, to encapsulate polynomial filters
*/

#include <slepc/private/stimpl.h>         /*I "slepcst.h" I*/
#include "filter.h"

const char *STFilterTypes[] = {"","FILTLAN","CHEBYSHEV","STFilterType","ST_FILTER_",NULL};
const char *STFilterDampings[] = {"NONE","JACKSON","LANCZOS","FEJER","STFilterDamping","ST_FILTER_DAMPING_",NULL};

static PetscErrorCode STFilterSetType_Private(ST st,STFilterType type)
{
  ST_FILTER    *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  ctx->type = type;
  switch(type) {
    case ST_FILTER_FILTLAN:
      PetscCall(STCreate_Filter_FILTLAN(st));
      break;
    case ST_FILTER_CHEBYSHEV:
      PetscCall(STCreate_Filter_Chebyshev(st));
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Invalid type");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Operator (filter):
               Op               P         M
   if nmat=1:  p(A)             NULL      p(A)
*/
static PetscErrorCode STComputeOperator_Filter(ST st)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  PetscCheck(st->nmat==1,PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Only implemented for standard eigenvalue problem");
  PetscCheck(ctx->intb<PETSC_MAX_REAL || ctx->inta>PETSC_MIN_REAL,PetscObjectComm((PetscObject)st),PETSC_ERR_ORDER,"Must pass an interval with STFilterSetInterval()");
  PetscCheck(ctx->right!=0.0 || ctx->left!=0.0,PetscObjectComm((PetscObject)st),PETSC_ERR_ORDER,"Must pass an approximate numerical range with STFilterSetRange()");
  PetscCheck(ctx->left<=ctx->inta && ctx->right>=ctx->intb,PetscObjectComm((PetscObject)st),PETSC_ERR_USER_INPUT,"The requested interval [%g,%g] must be contained in the numerical range [%g,%g]",(double)ctx->inta,(double)ctx->intb,(double)ctx->left,(double)ctx->right);

  if (!ctx->type) PetscCall(STFilterSetType_Private(st,ST_FILTER_FILTLAN)); /* default type */
  if (!ctx->polyDegree) ctx->polyDegree = 100;
  PetscCall(ctx->computeoperator(st,&st->T[0]));
  st->M = st->T[0];
  PetscCall(MatDestroy(&st->P));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STSetFromOptions_Filter(ST st,PetscOptionItems PetscOptionsObject)
{
  ST_FILTER        *ctx = (ST_FILTER*)st->data;
  PetscReal        array[2]={0,0};
  PetscInt         k;
  PetscBool        flg;
  STFilterType     type;
  STFilterDamping  damping;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"ST Filter Options");

    PetscCall(PetscOptionsEnum("-st_filter_type","How to construct the filter","STFilterSetType",STFilterTypes,(PetscEnum)ctx->type,(PetscEnum*)&type,&flg));
    if (flg) PetscCall(STFilterSetType(st,type));

    k = 2;
    PetscCall(PetscOptionsRealArray("-st_filter_interval","Interval containing the desired eigenvalues (two real values separated with a comma without spaces)","STFilterSetInterval",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Must pass two values in -st_filter_interval (comma-separated without spaces)");
      PetscCall(STFilterSetInterval(st,array[0],array[1]));
    }
    k = 2;
    PetscCall(PetscOptionsRealArray("-st_filter_range","Interval containing all eigenvalues (two real values separated with a comma without spaces)","STFilterSetRange",array,&k,&flg));
    if (flg) {
      PetscCheck(k>1,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Must pass two values in -st_filter_range (comma-separated without spaces)");
      PetscCall(STFilterSetRange(st,array[0],array[1]));
    }
    PetscCall(PetscOptionsInt("-st_filter_degree","Degree of filter polynomial","STFilterSetDegree",100,&k,&flg));
    if (flg) PetscCall(STFilterSetDegree(st,k));

    PetscCall(PetscOptionsEnum("-st_filter_damping","Type of damping","STFilterSetDamping",STFilterDampings,(PetscEnum)ctx->damping,(PetscEnum*)&damping,&flg));
    if (flg) PetscCall(STFilterSetDamping(st,damping));

  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STReset_Filter(ST st)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  ctx->left  = 0.0;
  ctx->right = 0.0;
  PetscCall(MatDestroy(&ctx->T));
  PetscCall(MatDestroyMatrices(ctx->nW,&ctx->W));
  ctx->nW = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterSetType_Filter(ST st,STFilterType type)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (ctx->type != type) {
    PetscCall(STReset_Filter(st));
    PetscCall(STFilterSetType_Private(st,type));
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
    ctx->filtch = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterSetType - Sets the method to be used to build the polynomial filter.

   Logically Collective

   Input Parameters:
+  st   - the spectral transformation context
-  type - the type of filter

   Options Database Key:
.  -st_filter_type \<type\> - set the type of filter

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterGetType()`
@*/
PetscErrorCode STFilterSetType(ST st,STFilterType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveEnum(st,type,2);
  PetscTryMethod(st,"STFilterSetType_C",(ST,STFilterType),(st,type));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetType_Filter(ST st,STFilterType *type)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  *type = ctx->type ;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetType - Gets the method to be used to build the polynomial filter.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  type - the type of filter

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterSetType()`
@*/
PetscErrorCode STFilterGetType(ST st,STFilterType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(type,2);
  PetscUseMethod(st,"STFilterGetType_C",(ST,STFilterType*),(st,type));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterSetInterval_Filter(ST st,PetscReal inta,PetscReal intb)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  PetscCheck(inta<intb,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be inta<intb");
  if (ctx->inta != inta || ctx->intb != intb) {
    ctx->inta   = inta;
    ctx->intb   = intb;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
    ctx->filtch = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterSetInterval - Defines the interval containing the desired eigenvalues.

   Logically Collective

   Input Parameters:
+  st   - the spectral transformation context
.  inta - left end of the interval
-  intb - right end of the interval

   Options Database Key:
.  -st_filter_interval <a,b> - set $[a,b]$ as the interval of interest

   Notes:
   The filter will be configured to emphasize eigenvalues contained in the given
   interval, and damp out eigenvalues outside it. If the interval is open, then
   the filter is low- or high-pass, otherwise it is mid-pass.

   Common usage is to set the interval in `EPS` with `EPSSetInterval()`.

   The interval must be contained within the numerical range of the matrix, see
   `STFilterSetRange()`.

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterGetInterval()`, `STFilterSetRange()`, `EPSSetInterval()`
@*/
PetscErrorCode STFilterSetInterval(ST st,PetscReal inta,PetscReal intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveReal(st,inta,2);
  PetscValidLogicalCollectiveReal(st,intb,3);
  PetscTryMethod(st,"STFilterSetInterval_C",(ST,PetscReal,PetscReal),(st,inta,intb));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetInterval_Filter(ST st,PetscReal *inta,PetscReal *intb)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (inta) *inta = ctx->inta;
  if (intb) *intb = ctx->intb;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetInterval - Gets the interval containing the desired eigenvalues.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameters:
+  inta - left end of the interval
-  intb - right end of the interval

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterSetInterval()`
@*/
PetscErrorCode STFilterGetInterval(ST st,PetscReal *inta,PetscReal *intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscUseMethod(st,"STFilterGetInterval_C",(ST,PetscReal*,PetscReal*),(st,inta,intb));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterSetRange_Filter(ST st,PetscReal left,PetscReal right)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  PetscCheck(left<right,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be left<right");
  if (ctx->left != left || ctx->right != right) {
    ctx->left   = left;
    ctx->right  = right;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
    ctx->filtch = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterSetRange - Defines the numerical range (or field of values) of the matrix, that is,
   the interval containing all eigenvalues.

   Logically Collective

   Input Parameters:
+  st    - the spectral transformation context
.  left  - left end of the spectral range
-  right - right end of the spectral range

   Options Database Key:
.  -st_filter_range <lmin,lmax> - set $[\lambda_\mathrm{min},\lambda_\mathrm{max}]$ as the numerical range

   Notes:
   The filter will be most effective if the numerical range is tight, that is,
   `left` and `right` are good approximations to the leftmost and rightmost
   eigenvalues, respectively.

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterGetRange()`, `STFilterSetInterval()`
@*/
PetscErrorCode STFilterSetRange(ST st,PetscReal left,PetscReal right)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveReal(st,left,2);
  PetscValidLogicalCollectiveReal(st,right,3);
  PetscTryMethod(st,"STFilterSetRange_C",(ST,PetscReal,PetscReal),(st,left,right));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetRange_Filter(ST st,PetscReal *left,PetscReal *right)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (left)  *left  = ctx->left;
  if (right) *right = ctx->right;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetRange - Gets the interval containing all eigenvalues.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameters:
+  left  - left end of the spectral range
-  right - right end of the spectral range

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterSetRange()`
@*/
PetscErrorCode STFilterGetRange(ST st,PetscReal *left,PetscReal *right)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscUseMethod(st,"STFilterGetRange_C",(ST,PetscReal*,PetscReal*),(st,left,right));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterSetDegree_Filter(ST st,PetscInt deg)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (deg == PETSC_DEFAULT || deg == PETSC_DECIDE) {
    ctx->polyDegree = 0;
    st->state       = ST_STATE_INITIAL;
    st->opready     = PETSC_FALSE;
    ctx->filtch     = PETSC_TRUE;
  } else {
    PetscCheck(deg>0,PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->polyDegree != deg) {
      ctx->polyDegree = deg;
      st->state       = ST_STATE_INITIAL;
      st->opready     = PETSC_FALSE;
      ctx->filtch     = PETSC_TRUE;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterSetDegree - Sets the degree of the filter polynomial.

   Logically Collective

   Input Parameters:
+  st  - the spectral transformation context
-  deg - polynomial degree

   Options Database Key:
.  -st_filter_degree \<deg\> - sets the degree of the filter polynomial

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterGetDegree()`
@*/
PetscErrorCode STFilterSetDegree(ST st,PetscInt deg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,deg,2);
  PetscTryMethod(st,"STFilterSetDegree_C",(ST,PetscInt),(st,deg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetDegree_Filter(ST st,PetscInt *deg)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  *deg = ctx->polyDegree;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetDegree - Gets the degree of the filter polynomial.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  deg - polynomial degree

   Level: intermediate

.seealso: [](ch:st), `STFILTER`, `STFilterSetDegree()`
@*/
PetscErrorCode STFilterGetDegree(ST st,PetscInt *deg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(deg,2);
  PetscUseMethod(st,"STFilterGetDegree_C",(ST,PetscInt*),(st,deg));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetThreshold_Filter(ST st,PetscReal *gamma)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (ctx->getthreshold) PetscCall(ctx->getthreshold(st,gamma));
  else *gamma = 0.5;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetThreshold - Gets the threshold value $\gamma$ used to decide
   whether eigenvalue approximations are inside or outside the wanted interval.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  gamma - the threshold value

   Note:
   An eigenvalue $\lambda$ is considered to be inside the wanted interval
   if $|p(\lambda)|\geq\gamma$, and outside otherwise.

   Level: developer

.seealso: [](ch:st), `STFILTER`, `STFilterGetRange()`
@*/
PetscErrorCode STFilterGetThreshold(ST st,PetscReal *gamma)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(gamma,2);
  PetscUseMethod(st,"STFilterGetThreshold_C",(ST,PetscReal*),(st,gamma));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterSetDamping_Filter(ST st,STFilterDamping damping)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (ctx->damping != damping) {
    ctx->damping = damping;
    st->state    = ST_STATE_INITIAL;
    st->opready  = PETSC_FALSE;
    ctx->filtch  = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterSetDamping - Sets the type of damping to be used in the polynomial filter.

   Logically Collective

   Input Parameters:
+  st      - the spectral transformation context
-  damping - the type of damping

   Options Database Key:
.  -st_filter_damping \<damping\> - sets the type of damping

   Note:
   Only used in `ST_FILTER_CHEBYSHEV` filters.

   Level: advanced

.seealso: [](ch:st), `STFILTER`, `STFilterGetDamping()`
@*/
PetscErrorCode STFilterSetDamping(ST st,STFilterDamping damping)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveEnum(st,damping,2);
  PetscTryMethod(st,"STFilterSetDamping_C",(ST,STFilterDamping),(st,damping));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STFilterGetDamping_Filter(ST st,STFilterDamping *damping)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  *damping = ctx->damping;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   STFilterGetDamping - Gets the type of damping used in the polynomial filter.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  damping - the type of damping

   Level: advanced

.seealso: [](ch:st), `STFILTER`, `STFilterSetDamping()`
@*/
PetscErrorCode STFilterGetDamping(ST st,STFilterDamping *damping)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscAssertPointer(damping,2);
  PetscUseMethod(st,"STFilterGetDamping_C",(ST,STFilterDamping*),(st,damping));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STView_Filter(ST st,PetscViewer viewer)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;
  PetscReal gamma;
  PetscBool isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"  filter type: %s\n",STFilterTypes[ctx->type]));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  interval of desired eigenvalues: [%g,%g]\n",(double)ctx->inta,(double)ctx->intb));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  numerical range: [%g,%g]\n",(double)ctx->left,(double)ctx->right));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  degree of filter polynomial: %" PetscInt_FMT "\n",ctx->polyDegree));
    if (ctx->damping && ctx->type==ST_FILTER_CHEBYSHEV) PetscCall(PetscViewerASCIIPrintf(viewer,"  type of damping = %s\n",STFilterDampings[ctx->damping]));
    if (st->state>=ST_STATE_SETUP && ctx->getthreshold) {
      PetscCall(STFilterGetThreshold(st,&gamma));
      PetscCall(PetscViewerASCIIPrintf(viewer,"  limit to accept eigenvalues: gamma=%g\n",(double)gamma));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode STDestroy_Filter(ST st)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (ctx->destroy) PetscCall(ctx->destroy(st));
  PetscCall(PetscFree(st->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetType_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetInterval_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetInterval_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetRange_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetRange_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetDegree_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetThreshold_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetDamping_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetDamping_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   STFILTER - STFILTER = "filter" - A special type of `ST` that represents
   a polynomial filter.

   Level: beginner

   Notes:
   In standard eigenvalue problems, when the eigenvalues of interest are interior
   to the spectrum and we want to avoid the high cost associated with the matrix
   factorization of the shift-and-invert spectral transformation, an alternative
   is to build a high-order polynomial such that $p(A)$ enhances the wanted
   eigenvalues and filters out the unwanted ones.

   The definition of the filter is done with functions such as `STFilterSetDegree()`,
   `STFilterSetInterval()` or `STFilterSetType()`.

.seealso: [](ch:st), `ST`, `STType`, `STSetType()`, `STSetMatrices()`, `STFilterSetDegree()`, `STFilterSetInterval()`, `STFilterSetType()`
M*/

SLEPC_EXTERN PetscErrorCode STCreate_Filter(ST st)
{
  ST_FILTER   *ctx;

  PetscFunctionBegin;
  PetscCall(PetscNew(&ctx));
  st->data = (void*)ctx;

  st->usesksp = PETSC_FALSE;

  ctx->type               = (STFilterType)0;
  ctx->inta               = PETSC_MIN_REAL;
  ctx->intb               = PETSC_MAX_REAL;
  ctx->left               = 0.0;
  ctx->right              = 0.0;
  ctx->polyDegree         = 0;

  st->ops->apply           = STApply_Generic;
  st->ops->computeoperator = STComputeOperator_Filter;
  st->ops->setfromoptions  = STSetFromOptions_Filter;
  st->ops->destroy         = STDestroy_Filter;
  st->ops->reset           = STReset_Filter;
  st->ops->view            = STView_Filter;

  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetType_C",STFilterSetType_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetType_C",STFilterGetType_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetInterval_C",STFilterSetInterval_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetInterval_C",STFilterGetInterval_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetRange_C",STFilterSetRange_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetRange_C",STFilterGetRange_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetDegree_C",STFilterSetDegree_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetDegree_C",STFilterGetDegree_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetThreshold_C",STFilterGetThreshold_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterSetDamping_C",STFilterSetDamping_Filter));
  PetscCall(PetscObjectComposeFunction((PetscObject)st,"STFilterGetDamping_C",STFilterGetDamping_Filter));
  PetscFunctionReturn(PETSC_SUCCESS);
}
