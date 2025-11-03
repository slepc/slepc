/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV orthogonalization routines
*/

#include <slepc/private/bvimpl.h>          /*I   "slepcbv.h"   I*/

/*
   BV_NormVecOrColumn - Compute the 2-norm of the working vector, irrespective of
   whether it is in a column or not
*/
static inline PetscErrorCode BV_NormVecOrColumn(BV bv,PetscInt j,Vec v,PetscReal *nrm)
{
  PetscFunctionBegin;
  if (v) PetscCall(BVNormVec(bv,v,NORM_2,nrm));
  else PetscCall(BVNormColumn(bv,j,NORM_2,nrm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BVDotColumnInc - Same as BVDotColumn() but also including column j, which
   is multiplied by itself
*/
static inline PetscErrorCode BVDotColumnInc(BV X,PetscInt j,PetscScalar *q)
{
  PetscInt       ksave;
  Vec            y;

  PetscFunctionBegin;
  PetscCall(PetscLogEventBegin(BV_DotVec,X,0,0,0));
  ksave = X->k;
  X->k = j+1;
  PetscCall(BVGetColumn(X,j,&y));
  PetscUseTypeMethod(X,dotvec,y,q);
  PetscCall(BVRestoreColumn(X,j,&y));
  X->k = ksave;
  PetscCall(PetscLogEventEnd(BV_DotVec,X,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BVOrthogonalizeMGS1 - Compute one step of Modified Gram-Schmidt
*/
static PetscErrorCode BVOrthogonalizeMGS1(BV bv,PetscInt j,Vec v,PetscBool *which,PetscScalar *h,PetscScalar *c,PetscReal *onrm,PetscReal *nrm)
{
  PetscInt          i;
  PetscScalar       dot;
  PetscBool         indef=bv->indef;
  Vec               vi,z,w=v;
  const PetscScalar *omega;

  PetscFunctionBegin;
  if (!v) PetscCall(BVGetColumn(bv,j,&w));
  if (onrm) PetscCall(BVNormVec(bv,w,NORM_2,onrm));
  z = w;
  if (indef) PetscCall(VecGetArrayRead(bv->omega,&omega));
  for (i=-bv->nc;i<j;i++) {
    if (which && i>=0 && !which[i]) continue;
    PetscCall(BVGetColumn(bv,i,&vi));
    /* h_i = (v, v_i) */
    if (bv->matrix) {
      PetscCall(BV_IPMatMult(bv,w));
      z = bv->Bx;
    }
    PetscCall(VecDot(z,vi,&dot));
    /* v <- v - h_i v_i */
    PetscCall(BV_SetValue(bv,i,0,c,dot));
    if (indef) dot /= PetscRealPart(omega[bv->nc+i]);
    PetscCall(VecAXPY(w,-dot,vi));
    PetscCall(BVRestoreColumn(bv,i,&vi));
  }
  if (nrm) PetscCall(BVNormVec(bv,w,NORM_2,nrm));
  if (!v) PetscCall(BVRestoreColumn(bv,j,&w));
  PetscCall(BV_AddCoefficients(bv,j,h,c));
  if (indef) PetscCall(VecRestoreArrayRead(bv->omega,&omega));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BVOrthogonalizeCGS1 - Compute |v'| (estimated), |v| and one step of CGS with
   only one global synchronization
*/
static PetscErrorCode BVOrthogonalizeCGS1(BV bv,PetscInt j,Vec v,PetscBool *which,PetscScalar *h,PetscScalar *c,PetscReal *onorm,PetscReal *norm)
{
  PetscReal      sum,beta;

  PetscFunctionBegin;
  /* h = W^* v ; alpha = (v, v) */
  bv->k = j;
  if (onorm || norm) {
    if (!v) {
      PetscCall(BVDotColumnInc(bv,j,c));
      PetscCall(BV_SquareRoot(bv,j,c,&beta));
    } else {
      PetscCall(BVDotVec(bv,v,c));
      PetscCall(BVNormVec(bv,v,NORM_2,&beta));
    }
  } else {
    if (!v) PetscCall(BVDotColumn(bv,j,c));
    else PetscCall(BVDotVec(bv,v,c));
  }

  /* q = v - V h */
  if (PetscUnlikely(bv->indef)) PetscCall(BV_ApplySignature(bv,j,c,PETSC_TRUE));
  if (!v) PetscCall(BVMultColumn(bv,-1.0,1.0,j,c));
  else PetscCall(BVMultVec(bv,-1.0,1.0,v,c));
  if (PetscUnlikely(bv->indef)) PetscCall(BV_ApplySignature(bv,j,c,PETSC_FALSE));

  /* compute |v| */
  if (onorm) *onorm = beta;

  if (norm) {
    if (PetscUnlikely(bv->indef)) PetscCall(BV_NormVecOrColumn(bv,j,v,norm));
    else {
      /* estimate |v'| from |v| */
      PetscCall(BV_SquareSum(bv,j,c,&sum));
      *norm = beta*beta-sum;
      if (PetscUnlikely(*norm <= 0.0)) PetscCall(BV_NormVecOrColumn(bv,j,v,norm));
      else *norm = PetscSqrtReal(*norm);
    }
  }
  PetscCall(BV_AddCoefficients(bv,j,h,c));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#define BVOrthogonalizeGS1(a,b,c,d,e,f,g,h) (bv->ops->gramschmidt?(*bv->ops->gramschmidt):(mgs?BVOrthogonalizeMGS1:BVOrthogonalizeCGS1))(a,b,c,d,e,f,g,h)

/*
   BVOrthogonalizeGS - Orthogonalize with (classical or modified) Gram-Schmidt

   j      - the index of the column to orthogonalize (cannot use both j and v)
   v      - the vector to orthogonalize (cannot use both j and v)
   which  - logical array indicating selected columns (only used in MGS)
   norm   - (optional) norm of the vector after being orthogonalized
   lindep - (optional) flag indicating possible linear dependence
*/
static PetscErrorCode BVOrthogonalizeGS(BV bv,PetscInt j,Vec v,PetscBool *which,PetscReal *norm,PetscBool *lindep)
{
  PetscScalar    *h,*c,*omega;
  PetscReal      onrm,nrm;
  PetscInt       k,l;
  PetscBool      mgs,dolindep,signature;

  PetscFunctionBegin;
  if (v) {
    k = bv->k;
    h = bv->h;
    c = bv->c;
  } else {
    k = j;
    h = NULL;
    c = NULL;
  }

  mgs = (bv->orthog_type==BV_ORTHOG_MGS)? PETSC_TRUE: PETSC_FALSE;

  /* if indefinite inner product, skip the computation of lindep */
  if (bv->indef && lindep) *lindep = PETSC_FALSE;
  dolindep = (!bv->indef && lindep)? PETSC_TRUE: PETSC_FALSE;

  /* if indefinite and we are orthogonalizing a column, the norm must always be computed */
  signature = (bv->indef && !v)? PETSC_TRUE: PETSC_FALSE;

  PetscCall(BV_CleanCoefficients(bv,k,h));

  switch (bv->orthog_ref) {

  case BV_ORTHOG_REFINE_IFNEEDED:
    PetscCall(BVOrthogonalizeGS1(bv,k,v,which,h,c,&onrm,&nrm));
    /* repeat if ||q|| < eta ||h|| */
    l = 1;
    while (l<3 && nrm && PetscAbsReal(nrm) < bv->orthog_eta*PetscAbsReal(onrm)) {
      l++;
      if (mgs||bv->indef) onrm = nrm;
      PetscCall(BVOrthogonalizeGS1(bv,k,v,which,h,c,(mgs||bv->indef)?NULL:&onrm,&nrm));
    }
    /* linear dependence check: criterion not satisfied in the last iteration */
    if (dolindep) *lindep = PetscNot(nrm && PetscAbsReal(nrm) >= bv->orthog_eta*PetscAbsReal(onrm));
    break;

  case BV_ORTHOG_REFINE_NEVER:
    PetscCall(BVOrthogonalizeGS1(bv,k,v,which,h,c,NULL,NULL));
    /* compute ||v|| */
    if (norm || dolindep || signature) PetscCall(BV_NormVecOrColumn(bv,k,v,&nrm));
    /* linear dependence check: just test for exactly zero norm */
    if (dolindep) *lindep = PetscNot(nrm);
    break;

  case BV_ORTHOG_REFINE_ALWAYS:
    PetscCall(BVOrthogonalizeGS1(bv,k,v,which,h,c,NULL,NULL));
    PetscCall(BVOrthogonalizeGS1(bv,k,v,which,h,c,dolindep?&onrm:NULL,(norm||dolindep||signature)?&nrm:NULL));
    /* linear dependence check: criterion not satisfied in the second iteration */
    if (dolindep) *lindep = PetscNot(nrm && PetscAbsReal(nrm) >= bv->orthog_eta*PetscAbsReal(onrm));
    break;
  }
  if (signature) {
    PetscCall(VecGetArray(bv->omega,&omega));
    omega[bv->nc+k] = (nrm<0.0)? -1.0: 1.0;
    PetscCall(VecRestoreArray(bv->omega,&omega));
  }
  if (norm) {
    *norm = nrm;
    if (!v) { /* store norm value next to the orthogonalization coefficients */
      if (dolindep && *lindep) PetscCall(BV_SetValue(bv,k,k,h,0.0));
      else PetscCall(BV_SetValue(bv,k,k,h,nrm));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVOrthogonalizeVec - Orthogonalize a given vector with respect to all
   active columns.

   Collective

   Input Parameters:
+  bv     - the basis vectors context
-  v      - the vector

   Output Parameters:
+  H      - (optional) coefficients computed during orthogonalization
.  norm   - (optional) norm of the vector after being orthogonalized
-  lindep - (optional) flag indicating that refinement did not improve the quality
            of orthogonalization

   Notes:
   This function is equivalent to `BVOrthogonalizeColumn()` but orthogonalizes
   a vector as an argument rather than taking one of the `BV` columns. The
   vector is orthogonalized against all active columns (`k`) and the constraints.
   If `H` is given, it must have enough space to store `k-l` coefficients, where `l`
   is the number of leading columns.

   In the case of an indefinite inner product, the `lindep` parameter is not
   computed (set to false).

   Level: advanced

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVSetOrthogonalization()`, `BVSetActiveColumns()`, `BVGetNumConstraints()`
@*/
PetscErrorCode BVOrthogonalizeVec(BV bv,Vec v,PetscScalar H[],PetscReal *norm,PetscBool *lindep)
{
  PetscInt       ksave,lsave;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidType(v,2);
  PetscCheckSameComm(bv,1,v,2);

  PetscCall(PetscLogEventBegin(BV_OrthogonalizeVec,bv,0,0,0));
  ksave = bv->k;
  lsave = bv->l;
  bv->l = -bv->nc;  /* must also orthogonalize against constraints and leading columns */
  PetscCall(BV_AllocateCoeffs(bv));
  PetscCall(BV_AllocateSignature(bv));
  PetscCall(BVOrthogonalizeGS(bv,0,v,NULL,norm,lindep));
  bv->k = ksave;
  bv->l = lsave;
  if (H) PetscCall(BV_StoreCoefficients(bv,bv->k,bv->h,H));
  PetscCall(PetscLogEventEnd(BV_OrthogonalizeVec,bv,0,0,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVOrthogonalizeColumn - Orthogonalize one of the column vectors with respect to
   the previous ones.

   Collective

   Input Parameters:
+  bv     - the basis vectors context
-  j      - index of column to be orthogonalized

   Output Parameters:
+  H      - (optional) coefficients computed during orthogonalization
.  norm   - (optional) norm of the vector after being orthogonalized
-  lindep - (optional) flag indicating that refinement did not improve the quality
            of orthogonalization

   Notes:
   This function applies an orthogonal projector to project vector $v_j$ onto
   the orthogonal complement of the span of the columns $v_{0:j-1}$.
   The columns $v_{0:j-1}$ are assumed to be mutually orthonormal.

   Leading columns $v_{0:\ell-1}$ also participate in the orthogonalization, as well
   as the constraints. If `H` is given, it must have enough space to store
   $j-\ell+1$ coefficients (the last coefficient will contain the value `norm`, unless
   the `norm` argument is `NULL`).

   If a non-standard inner product has been specified with `BVSetMatrix()`,
   then the vector is B-orthogonalized, using the non-standard inner product
   defined by matrix $B$. The output vector satisfies $v_j^*Bv_{0:j-1} = 0$.

   This routine does not normalize the resulting vector, see `BVOrthonormalizeColumn()`.

   In the case of an indefinite inner product, the `lindep` parameter is not
   computed (set to false).

   Level: advanced

.seealso: [](sec:bv), `BVSetOrthogonalization()`, `BVSetMatrix()`, `BVSetActiveColumns()`, `BVOrthogonalize()`, `BVOrthogonalizeVec()`, `BVGetNumConstraints()`, `BVOrthonormalizeColumn()`
@*/
PetscErrorCode BVOrthogonalizeColumn(BV bv,PetscInt j,PetscScalar H[],PetscReal *norm,PetscBool *lindep)
{
  PetscInt       ksave,lsave;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,bv->m);

  PetscCall(PetscLogEventBegin(BV_OrthogonalizeVec,bv,0,0,0));
  ksave = bv->k;
  lsave = bv->l;
  bv->l = -bv->nc;  /* must also orthogonalize against constraints and leading columns */
  if (!bv->buffer) PetscCall(BVGetBufferVec(bv,&bv->buffer));
  PetscCall(BV_AllocateSignature(bv));
  PetscCall(BVOrthogonalizeGS(bv,j,NULL,NULL,norm,lindep));
  bv->k = ksave;
  bv->l = lsave;
  if (H) PetscCall(BV_StoreCoefficients(bv,j,NULL,H));
  PetscCall(PetscLogEventEnd(BV_OrthogonalizeVec,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVOrthonormalizeColumn - Orthonormalize one of the column vectors with respect to
   the previous ones.

   Collective

   Input Parameters:
+  bv      - the basis vectors context
.  j       - index of column to be orthonormalized
-  replace - whether it is allowed to set the vector randomly

   Output Parameters:
+  norm    - (optional) norm of the vector after orthogonalization and before normalization
-  lindep  - (optional) flag indicating that linear dependence was determined during
             orthogonalization

   Notes:
   This is equivalent to a call to `BVOrthogonalizeColumn()` followed by a
   call to `BVScaleColumn()` with the reciprocal of the norm.

   This function first orthogonalizes vector $v_j$ with respect to $v_{0:j-1}$.
   A byproduct of this computation is `norm`, the norm of the vector after
   orthogonalization. Secondly, it scales the vector with `1/norm`, so that the
   resulting vector has unit norm.

   If after orthogonalization the vector $v_j$ is exactly zero, it cannot be normalized
   because `norm=0`. In that case, it could be left as zero or replaced by a random
   vector that is then orthonormalized. The latter is achieved by setting the
   argument `replace` to `PETSC_TRUE`. The vector will be replaced by a random vector also
   if `lindep` was set to `PETSC_TRUE`, even if the norm is not exactly zero.

   If the vector has been replaced by a random vector, the output arguments `norm` and
   `lindep` will be set according to the orthogonalization of this new vector.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVScaleColumn()`
@*/
PetscErrorCode BVOrthonormalizeColumn(BV bv,PetscInt j,PetscBool replace,PetscReal *norm,PetscBool *lindep)
{
  PetscScalar    alpha;
  PetscReal      nrm;
  PetscInt       ksave,lsave;
  PetscBool      lndep;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,bv->m);

  /* orthogonalize */
  PetscCall(PetscLogEventBegin(BV_OrthogonalizeVec,bv,0,0,0));
  ksave = bv->k;
  lsave = bv->l;
  bv->l = -bv->nc;  /* must also orthogonalize against constraints and leading columns */
  if (!bv->buffer) PetscCall(BVGetBufferVec(bv,&bv->buffer));
  PetscCall(BV_AllocateSignature(bv));
  PetscCall(BVOrthogonalizeGS(bv,j,NULL,NULL,&nrm,&lndep));
  if (replace && (nrm==0.0 || lndep)) {
    PetscCall(PetscInfo(bv,"Vector was linearly dependent, generating a new random vector\n"));
    PetscCall(BVSetRandomColumn(bv,j));
    PetscCall(BVOrthogonalizeGS(bv,j,NULL,NULL,&nrm,&lndep));
    if (nrm==0.0 || lndep) {  /* yet another attempt */
      PetscCall(BVSetRandomColumn(bv,j));
      PetscCall(BVOrthogonalizeGS(bv,j,NULL,NULL,&nrm,&lndep));
    }
  }
  bv->k = ksave;
  bv->l = lsave;
  PetscCall(PetscLogEventEnd(BV_OrthogonalizeVec,bv,0,0,0));

  /* scale */
  if (nrm!=1.0 && nrm!=0.0) {
    alpha = 1.0/nrm;
    PetscCall(PetscLogEventBegin(BV_Scale,bv,0,0,0));
    PetscUseTypeMethod(bv,scale,j,alpha);
    PetscCall(PetscLogEventEnd(BV_Scale,bv,0,0,0));
  }
  if (norm) *norm = nrm;
  if (lindep) *lindep = lndep;
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVOrthogonalizeSomeColumn - Orthogonalize one of the column vectors with
   respect to some of the previous ones.

   Collective

   Input Parameters:
+  bv     - the basis vectors context
.  j      - index of column to be orthogonalized
-  which  - logical array indicating selected columns

   Output Parameters:
+  H      - (optional) coefficients computed during orthogonalization
.  norm   - (optional) norm of the vector after being orthogonalized
-  lindep - (optional) flag indicating that refinement did not improve the quality
            of orthogonalization

   Notes:
   This function is similar to `BVOrthogonalizeColumn()`, but $v_j$ is
   orthogonalized only against columns $v_i$ having `which[i]=PETSC_TRUE`.
   The length of array `which` must be `j` at least.

   The use of this operation is restricted to MGS orthogonalization type,
   `BV_ORTHOG_MGS`.

   In the case of an indefinite inner product, the `lindep` parameter is not
   computed (set to false).

   Level: advanced

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVSetOrthogonalization()`
@*/
PetscErrorCode BVOrthogonalizeSomeColumn(BV bv,PetscInt j,PetscBool *which,PetscScalar H[],PetscReal *norm,PetscBool *lindep)
{
  PetscInt       ksave,lsave;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscAssertPointer(which,3);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscCheck(j>=0,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  PetscCheck(j<bv->m,PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%" PetscInt_FMT " but BV only has %" PetscInt_FMT " columns",j,bv->m);
  PetscCheck(bv->orthog_type==BV_ORTHOG_MGS,PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Operation only available for MGS orthogonalization");

  PetscCall(PetscLogEventBegin(BV_OrthogonalizeVec,bv,0,0,0));
  ksave = bv->k;
  lsave = bv->l;
  bv->l = -bv->nc;  /* must also orthogonalize against constraints and leading columns */
  if (!bv->buffer) PetscCall(BVGetBufferVec(bv,&bv->buffer));
  PetscCall(BV_AllocateSignature(bv));
  PetscCall(BVOrthogonalizeGS(bv,j,NULL,which,norm,lindep));
  bv->k = ksave;
  bv->l = lsave;
  if (H) PetscCall(BV_StoreCoefficients(bv,j,NULL,H));
  PetscCall(PetscLogEventEnd(BV_OrthogonalizeVec,bv,0,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)bv));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   Block Gram-Schmidt: V2 = V2 - V1*R12, where R12 = V1'*V2
 */
static PetscErrorCode BVOrthogonalize_BlockGS(BV V,Mat R)
{
  BV             V1;

  PetscFunctionBegin;
  PetscCall(BVGetSplit(V,&V1,NULL));
  PetscCall(BVDot(V,V1,R));
  PetscCall(BVMult(V,-1.0,1.0,V1,R));
  PetscCall(BVRestoreSplit(V,&V1,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BV_ORTHOG_BLOCK_GS - Orthogonalize a set of vectors with Gram-Schmidt, column by column.

   Note:
   In this case, the block orthogonalization is done in the same way as the
   orthogonalization of individual vectors, e.g., `BVOrthogonalizeColumn`, and
   hence the settings provided via `BVSetOrthogonalization()` apply.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogBlockType`, `BVSetOrthogonalization()`, `BVOrthogonalize()`, `BV_ORTHOG_BLOCK_CHOL`, `BV_ORTHOG_BLOCK_TSQR`, `BV_ORTHOG_BLOCK_TSQRCHOL`, `BV_ORTHOG_BLOCK_SVQB`
M*/
static PetscErrorCode BVOrthogonalize_GS(BV V,Mat R)
{
  PetscScalar    *r=NULL;
  PetscReal      norm;
  PetscInt       j,ldr,lsave;
  Vec            v,w;

  PetscFunctionBegin;
  if (R) {
    PetscCall(MatDenseGetLDA(R,&ldr));
    PetscCall(MatDenseGetArray(R,&r));
  }
  if (V->matrix) {
    PetscCall(BVGetCachedBV(V,&V->cached));
    PetscCall(BVSetActiveColumns(V->cached,V->l,V->k));
  }
  for (j=V->l;j<V->k;j++) {
    if (V->matrix && V->orthog_type==BV_ORTHOG_MGS) {  /* fill cached BV */
      PetscCall(BVGetColumn(V->cached,j,&v));
      PetscCall(BVGetColumn(V,j,&w));
      PetscCall(MatMult(V->matrix,w,v));
      PetscCall(BVRestoreColumn(V,j,&w));
      PetscCall(BVRestoreColumn(V->cached,j,&v));
    }
    if (R) {
      PetscCall(BVOrthogonalizeColumn(V,j,NULL,&norm,NULL));
      lsave = V->l;
      V->l = -V->nc;
      PetscCall(BV_StoreCoefficients(V,j,NULL,r+j*ldr));
      V->l = lsave;
      r[j+j*ldr] = norm;
    } else PetscCall(BVOrthogonalizeColumn(V,j,NULL,&norm,NULL));
    PetscCheck(norm,PetscObjectComm((PetscObject)V),PETSC_ERR_CONV_FAILED,"Breakdown in BVOrthogonalize due to a linearly dependent column");
    if (V->matrix && V->orthog_type==BV_ORTHOG_CGS) {  /* fill cached BV */
      PetscCall(BVGetColumn(V->cached,j,&v));
      PetscCall(VecCopy(V->Bx,v));
      PetscCall(BVRestoreColumn(V->cached,j,&v));
    }
    PetscCall(BVScaleColumn(V,j,1.0/norm));
  }
  if (R) PetscCall(MatDenseRestoreArray(R,&r));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  BV_GetBufferMat - Create auxiliary seqdense matrix that wraps the bv->buffer.
*/
static inline PetscErrorCode BV_GetBufferMat(BV bv)
{
  PetscInt       ld;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (!bv->Abuffer) {
    if (!bv->buffer) PetscCall(BVGetBufferVec(bv,&bv->buffer));
    ld = bv->m+bv->nc;
    PetscCall(VecGetArray(bv->buffer,&array));
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF,ld,bv->m,array,&bv->Abuffer));
    PetscCall(VecRestoreArray(bv->buffer,&array));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   BV_StoreCoeffsBlock_Default - Copy the contents of the BV buffer to a dense Mat
   provided by the caller. Only columns l:k-1 are copied, restricting to the upper
   triangular part if tri=PETSC_TRUE.
*/
static inline PetscErrorCode BV_StoreCoeffsBlock_Default(BV bv,Mat R,PetscBool tri)
{
  const PetscScalar *bb;
  PetscScalar       *rr;
  PetscInt          j,ldr,ldb;

  PetscFunctionBegin;
  PetscCall(MatDenseGetLDA(R,&ldr));
  PetscCall(MatDenseGetArray(R,&rr));
  ldb  = bv->m+bv->nc;
  PetscCall(VecGetArrayRead(bv->buffer,&bb));
  for (j=bv->l;j<bv->k;j++) PetscCall(PetscArraycpy(rr+j*ldr,bb+j*ldb,(tri?(j+1):bv->k)+bv->nc));
  PetscCall(VecRestoreArrayRead(bv->buffer,&bb));
  PetscCall(MatDenseRestoreArray(R,&rr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BV_ORTHOG_BLOCK_CHOL - Orthogonalize a set of vectors $V$ with the Cholesky QR
   method, that is, compute $R=\mathrm{chol}(V^*V)$, then obtain the
   orthogonalized vectors as $Q=V R^{-1}$.

   Developer Note:
   This is plain Cholesky QR, which may fail in case of ill-conditioned $V$.
   An iterated version such as CholeskyQR2 should be preferred {cite:p}`Yam15`.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogBlockType`, `BVSetOrthogonalization()`, `BVOrthogonalize()`, `BV_ORTHOG_BLOCK_GS`, `BV_ORTHOG_BLOCK_TSQR`, `BV_ORTHOG_BLOCK_TSQRCHOL`, `BV_ORTHOG_BLOCK_SVQB`
M*/
static PetscErrorCode BVOrthogonalize_Chol(BV V,Mat Rin)
{
  Mat            R,S;

  PetscFunctionBegin;
  PetscCall(BV_GetBufferMat(V));
  R = V->Abuffer;
  if (Rin) S = Rin;   /* use Rin as a workspace for S */
  else S = R;
  if (V->l) PetscCall(BVOrthogonalize_BlockGS(V,R));
  PetscCall(BVDot(V,V,R));
  PetscCall(BVMatCholInv_LAPACK_Private(V,R,S));
  PetscCall(BVMultInPlace(V,S,V->l,V->k));
  if (Rin) PetscCall(BV_StoreCoeffsBlock_Default(V,Rin,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BV_ORTHOG_BLOCK_TSQR - Orthogonalize a set of vectors with the Tall-Skinny QR method {cite:p}`Dem12`.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogBlockType`, `BVSetOrthogonalization()`, `BVOrthogonalize()`, `BV_ORTHOG_BLOCK_GS`, `BV_ORTHOG_BLOCK_CHOL`, `BV_ORTHOG_BLOCK_TSQRCHOL`, `BV_ORTHOG_BLOCK_SVQB`
M*/
static PetscErrorCode BVOrthogonalize_TSQR(BV V,Mat Rin)
{
  PetscScalar    *pv,*r=NULL;
  PetscInt       ldr;
  Mat            R;

  PetscFunctionBegin;
  PetscCall(BV_GetBufferMat(V));
  R = V->Abuffer;
  if (V->l) PetscCall(BVOrthogonalize_BlockGS(V,R));
  PetscCall(MatDenseGetLDA(R,&ldr));
  PetscCall(MatDenseGetArray(R,&r));
  PetscCall(BVGetArray(V,&pv));
  PetscCall(BVOrthogonalize_LAPACK_TSQR(V,V->n,V->k-V->l,pv+(V->nc+V->l)*V->ld,V->ld,r+V->l*ldr+V->l,ldr));
  PetscCall(BVRestoreArray(V,&pv));
  PetscCall(MatDenseRestoreArray(R,&r));
  if (Rin) PetscCall(BV_StoreCoeffsBlock_Default(V,Rin,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BV_ORTHOG_BLOCK_TSQRCHOL - Orthogonalize a set of vectors $V$ with the TSQR
   method, but computing the triangular factor $R$ only, then retrieving the
   orthogonalized vectors as $Q=V R^{-1}$.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogBlockType`, `BVSetOrthogonalization()`, `BVOrthogonalize()`, `BV_ORTHOG_BLOCK_GS`, `BV_ORTHOG_BLOCK_CHOL`, `BV_ORTHOG_BLOCK_TSQR`, `BV_ORTHOG_BLOCK_SVQB`
M*/
static PetscErrorCode BVOrthogonalize_TSQRCHOL(BV V,Mat Rin)
{
  PetscScalar    *pv,*r=NULL;
  PetscInt       ldr;
  Mat            R,S;

  PetscFunctionBegin;
  PetscCall(BV_GetBufferMat(V));
  R = V->Abuffer;
  if (Rin) S = Rin;   /* use Rin as a workspace for S */
  else S = R;
  if (V->l) PetscCall(BVOrthogonalize_BlockGS(V,R));
  PetscCall(MatDenseGetLDA(R,&ldr));
  PetscCall(MatDenseGetArray(R,&r));
  PetscCall(BVGetArray(V,&pv));
  PetscCall(BVOrthogonalize_LAPACK_TSQR_OnlyR(V,V->n,V->k-V->l,pv+(V->nc+V->l)*V->ld,V->ld,r+V->l*ldr+V->l,ldr));
  PetscCall(BVRestoreArray(V,&pv));
  PetscCall(MatDenseRestoreArray(R,&r));
  PetscCall(BVMatTriInv_LAPACK_Private(V,R,S));
  PetscCall(BVMultInPlace(V,S,V->l,V->k));
  if (Rin) PetscCall(BV_StoreCoeffsBlock_Default(V,Rin,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*MC
   BV_ORTHOG_BLOCK_SVQB - Orthogonalize a set of vectors with the SVQB method {cite:p}`Sta02`.

   Note:
   As opposed to other block orthogonalization methods, SVQB does not compute a
   QR decomposition but a QB decomposition, i.e., the second factor is not
   triangular.

   Level: advanced

.seealso: [](sec:bv), `BVOrthogBlockType`, `BVSetOrthogonalization()`, `BVOrthogonalize()`, `BV_ORTHOG_BLOCK_GS`, `BV_ORTHOG_BLOCK_CHOL`, `BV_ORTHOG_BLOCK_TSQR`, `BV_ORTHOG_BLOCK_TSQRCHOL`
M*/
static PetscErrorCode BVOrthogonalize_SVQB(BV V,Mat Rin)
{
  Mat            R,S;

  PetscFunctionBegin;
  PetscCall(BV_GetBufferMat(V));
  R = V->Abuffer;
  if (Rin) S = Rin;   /* use Rin as a workspace for S */
  else S = R;
  if (V->l) PetscCall(BVOrthogonalize_BlockGS(V,R));
  PetscCall(BVDot(V,V,R));
  PetscCall(BVMatSVQB_LAPACK_Private(V,R,S));
  PetscCall(BVMultInPlace(V,S,V->l,V->k));
  if (Rin) PetscCall(BV_StoreCoeffsBlock_Default(V,Rin,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVOrthogonalize - Orthogonalize all columns (starting from the leading ones),
   that is, compute the QR decomposition.

   Collective

   Input Parameters:
+  V - basis vectors to be orthogonalized (or $B$-orthogonalized), modified on output
-  R - a sequential dense matrix (or `NULL`), on output the triangular factor of
       the QR decomposition

   Notes:
   On input, matrix `R` must be a square sequential dense `Mat`, with at least as many
   rows and columns as the number of active columns of `V`. The output satisfies
   $V = \tilde V R$ (where $\tilde V$ represents `V` on exit) and $\tilde V^*\tilde V = I$
   (or $\tilde V^*B\tilde V = I$ if an inner product matrix $B$ has been specified with
   `BVSetMatrix()`).

   If `V` has leading columns, then they are not modified (are assumed to be already
   orthonormal) and the leading columns of `R` are not referenced. Let the
   decomposition be
   $$
     \begin{bmatrix}
       V_1 & V_2
     \end{bmatrix}=
     \begin{bmatrix}
       \tilde V_1 & \tilde V_2
     \end{bmatrix}
     \begin{bmatrix}
       R_{11} & R_{12} \\
            0 & R_{22}
     \end{bmatrix}
   $$
   then $\tilde V_1$ is left unchanged (equal to $V_1$) as well as $R_{11}$ (they should
   satisfy $V_1 = \tilde V_1 R_{11}$).

   Can pass `NULL` if `R` is not required.

   The method to be used for block orthogonalization can be set with
   `BVSetOrthogonalization()`. If set to `BV_ORTHOG_BLOCK_GS`, the computation is done
   column by column with successive calls to `BVOrthogonalizeColumn()`. Note that in the
   `BV_ORTHOG_BLOCK_SVQB` method the `R` factor is not upper triangular.

   If `V` is rank-deficient or very ill-conditioned, that is, one or more columns are
   (almost) linearly dependent with respect to the rest, then the algorithm may
   break down or result in larger numerical error. Linearly dependent columns are
   essentially replaced by random directions, and the corresponding diagonal entry
   in `R` is set to (nearly) zero.

   Level: intermediate

.seealso: [](sec:bv), `BVOrthogonalizeColumn()`, `BVOrthogonalizeVec()`, `BVSetMatrix()`, `BVSetActiveColumns()`, `BVSetOrthogonalization()`, `BVOrthogBlockType`
@*/
PetscErrorCode BVOrthogonalize(BV V,Mat R)
{
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  if (R) {
    PetscValidHeaderSpecific(R,MAT_CLASSID,2);
    PetscValidType(R,2);
    PetscCheckTypeName(R,MATSEQDENSE);
    PetscCall(MatGetSize(R,&m,&n));
    PetscCheck(m==n,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument is not square, it has %" PetscInt_FMT " rows and %" PetscInt_FMT " columns",m,n);
    PetscCheck(n>=V->k,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat size %" PetscInt_FMT " is smaller than the number of BV active columns %" PetscInt_FMT,n,V->k);
  }
  PetscCheck(!V->nc,PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Not implemented for BV with constraints, use BVOrthogonalizeColumn() instead");

  PetscCall(PetscLogEventBegin(BV_Orthogonalize,V,R,0,0));
  switch (V->orthog_block) {
  case BV_ORTHOG_BLOCK_GS: /* proceed column by column with Gram-Schmidt */
    PetscCall(BVOrthogonalize_GS(V,R));
    break;
  case BV_ORTHOG_BLOCK_CHOL:
    PetscCall(BVOrthogonalize_Chol(V,R));
    break;
  case BV_ORTHOG_BLOCK_TSQR:
    PetscCheck(!V->matrix,PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Orthogonalization method not available for non-standard inner product");
    PetscCall(BVOrthogonalize_TSQR(V,R));
    break;
  case BV_ORTHOG_BLOCK_TSQRCHOL:
    PetscCheck(!V->matrix,PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Orthogonalization method not available for non-standard inner product");
    PetscCall(BVOrthogonalize_TSQRCHOL(V,R));
    break;
  case BV_ORTHOG_BLOCK_SVQB:
    PetscCall(BVOrthogonalize_SVQB(V,R));
    break;
  }
  PetscCall(PetscLogEventEnd(BV_Orthogonalize,V,R,0,0));
  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}
