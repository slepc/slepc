/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/vecimplslepc.h>       /*I "slepcvec.h" I*/

/*@
   VecNormalizeComplex - Normalizes a possibly complex vector by the 2-norm.

   Collective

   Input Parameters:
+  xr - the real part of the vector (overwritten on output)
.  xi - the imaginary part of the vector (not referenced if iscomplex is false)
-  iscomplex - a flag indicating if the vector is complex

   Output Parameter:
.  norm - the vector norm before normalization (can be set to NULL)

   Level: developer

.seealso: BVNormalize()
@*/
PetscErrorCode VecNormalizeComplex(Vec xr,Vec xi,PetscBool iscomplex,PetscReal *norm)
{
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      normr,normi,alpha;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(xr,VEC_CLASSID,1);
#if !defined(PETSC_USE_COMPLEX)
  if (iscomplex) {
    PetscValidHeaderSpecific(xi,VEC_CLASSID,2);
    PetscCall(VecNormBegin(xr,NORM_2,&normr));
    PetscCall(VecNormBegin(xi,NORM_2,&normi));
    PetscCall(VecNormEnd(xr,NORM_2,&normr));
    PetscCall(VecNormEnd(xi,NORM_2,&normi));
    alpha = SlepcAbsEigenvalue(normr,normi);
    if (norm) *norm = alpha;
    alpha = 1.0 / alpha;
    PetscCall(VecScale(xr,alpha));
    PetscCall(VecScale(xi,alpha));
  } else
#endif
    PetscCall(VecNormalize(xr,norm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode VecCheckOrthogonality_Private(Vec V[],PetscInt nv,Vec W[],PetscInt nw,Mat B,PetscViewer viewer,PetscReal *lev,PetscBool norm)
{
  PetscInt       i,j;
  PetscScalar    *vals;
  PetscBool      isascii;
  Vec            w;

  PetscFunctionBegin;
  if (!lev) {
    if (!viewer) PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)*V),&viewer));
    PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,6);
    PetscCheckSameComm(*V,1,viewer,6);
    PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
    if (!isascii) PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscCall(PetscMalloc1(nv,&vals));
  if (B) PetscCall(VecDuplicate(V[0],&w));
  if (lev) *lev = 0.0;
  for (i=0;i<nw;i++) {
    if (B) {
      if (W) PetscCall(MatMult(B,W[i],w));
      else PetscCall(MatMult(B,V[i],w));
    } else {
      if (W) w = W[i];
      else w = V[i];
    }
    PetscCall(VecMDot(w,nv,V,vals));
    for (j=0;j<nv;j++) {
      if (lev) {
        if (i!=j) *lev = PetscMax(*lev,PetscAbsScalar(vals[j]));
        else if (norm) {
          if (PetscRealPart(vals[j])<0.0) *lev = PetscMax(*lev,PetscAbsScalar(vals[j]+PetscRealConstant(1.0)));  /* indefinite case */
          else *lev = PetscMax(*lev,PetscAbsScalar(vals[j]-PetscRealConstant(1.0)));
        }
      } else {
#if !defined(PETSC_USE_COMPLEX)
        PetscCall(PetscViewerASCIIPrintf(viewer," %12g  ",(double)vals[j]));
#else
        PetscCall(PetscViewerASCIIPrintf(viewer," %12g%+12gi ",(double)PetscRealPart(vals[j]),(double)PetscImaginaryPart(vals[j])));
#endif
      }
    }
    if (!lev) PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(PetscFree(vals));
  if (B) PetscCall(VecDestroy(&w));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecCheckOrthogonality - Checks (or prints) the level of (bi-)orthogonality
   of a set of vectors.

   Collective

   Input Parameters:
+  V  - a set of vectors
.  nv - number of V vectors
.  W  - an alternative set of vectors (optional)
.  nw - number of W vectors
.  B  - Hermitian matrix defining the inner product (optional)
-  viewer - optional visualization context

   Output Parameter:
.  lev - level of orthogonality (optional)

   Notes:
   This function computes W'*V and prints the result. It is intended to check
   the level of bi-orthogonality of the vectors in the two sets. If W is equal
   to NULL then V is used, thus checking the orthogonality of the V vectors.

   If matrix B is provided then the check uses the B-inner product, W'*B*V,
   where B is assumed to be Hermitian.

   If V, W represent eigenvectors computed by SLEPc, this function will not work
   correctly if one of the eigenvalues is complex when running with real scalars.

   If lev is not NULL, it will contain the maximum entry of matrix
   W'*V - I (in absolute value) omitting the diagonal. Otherwise, the matrix W'*V
   is printed.

   Level: developer

.seealso: VecCheckOrthonormality()
@*/
PetscErrorCode VecCheckOrthogonality(Vec V[],PetscInt nv,Vec W[],PetscInt nw,Mat B,PetscViewer viewer,PetscReal *lev)
{
  PetscFunctionBegin;
  PetscAssertPointer(V,1);
  PetscValidHeaderSpecific(*V,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(*V,nv,2);
  PetscValidLogicalCollectiveInt(*V,nw,4);
  if (nv<=0 || nw<=0) PetscFunctionReturn(PETSC_SUCCESS);
  if (W) {
    PetscAssertPointer(W,3);
    PetscValidHeaderSpecific(*W,VEC_CLASSID,3);
    PetscCheckSameComm(*V,1,*W,3);
  }
  PetscCall(VecCheckOrthogonality_Private(V,nv,W,nw,B,viewer,lev,PETSC_FALSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecCheckOrthonormality - Checks (or prints) the level of (bi-)orthonormality
   of a set of vectors.

   Collective

   Input Parameters:
+  V  - a set of vectors
.  nv - number of V vectors
.  W  - an alternative set of vectors (optional)
.  nw - number of W vectors
.  B  - Hermitian matrix defining the inner product (optional)
-  viewer - optional visualization context

   Output Parameter:
.  lev - level of orthogonality (optional)

   Notes:
   This function is equivalent to VecCheckOrthonormality(), but in addition it checks
   that the diagonal of W'*V (or W'*B*V) is equal to all ones.

   Level: developer

.seealso: VecCheckOrthogonality()
@*/
PetscErrorCode VecCheckOrthonormality(Vec V[],PetscInt nv,Vec W[],PetscInt nw,Mat B,PetscViewer viewer,PetscReal *lev)
{
  PetscFunctionBegin;
  PetscAssertPointer(V,1);
  PetscValidHeaderSpecific(*V,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(*V,nv,2);
  PetscValidLogicalCollectiveInt(*V,nw,4);
  if (nv<=0 || nw<=0) PetscFunctionReturn(PETSC_SUCCESS);
  if (W) {
    PetscAssertPointer(W,3);
    PetscValidHeaderSpecific(*W,VEC_CLASSID,3);
    PetscCheckSameComm(*V,1,*W,3);
  }
  PetscCall(VecCheckOrthogonality_Private(V,nv,W,nw,B,viewer,lev,PETSC_TRUE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecDuplicateEmpty - Creates a new vector of the same type as an existing vector,
   but without internal array.

   Collective

   Input Parameters:
.  v - a vector to mimic

   Output Parameter:
.  newv - location to put new vector

   Note:
   This is similar to VecDuplicate(), but the new vector does not have an internal
   array, so the intended usage is with VecPlaceArray().

   Level: developer

.seealso: MatCreateVecsEmpty()
@*/
PetscErrorCode VecDuplicateEmpty(Vec v,Vec *newv)
{
  PetscBool      standard,cuda,hip,mpi;
  PetscInt       N,nloc,bs;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  PetscAssertPointer(newv,2);
  PetscValidType(v,1);

  PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&standard,VECSEQ,VECMPI,""));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&cuda,VECSEQCUDA,VECMPICUDA,""));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&hip,VECSEQHIP,VECMPIHIP,""));
  if (standard || cuda || hip) {
    PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&mpi,VECMPI,VECMPICUDA,VECMPIHIP,""));
    PetscCall(VecGetLocalSize(v,&nloc));
    PetscCall(VecGetSize(v,&N));
    PetscCall(VecGetBlockSize(v,&bs));
    if (cuda) {
#if defined(PETSC_HAVE_CUDA)
      if (mpi) PetscCall(VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)v),bs,nloc,N,NULL,newv));
      else PetscCall(VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)v),bs,N,NULL,newv));
#endif
    } else if (hip) {
#if defined(PETSC_HAVE_HIP)
      if (mpi) PetscCall(VecCreateMPIHIPWithArray(PetscObjectComm((PetscObject)v),bs,nloc,N,NULL,newv));
      else PetscCall(VecCreateSeqHIPWithArray(PetscObjectComm((PetscObject)v),bs,N,NULL,newv));
#endif
    } else {
      if (mpi) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)v),bs,nloc,N,NULL,newv));
      else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)v),bs,N,NULL,newv));
    }
  } else PetscCall(VecDuplicate(v,newv)); /* standard duplicate, with internal array */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   VecSetRandomNormal - Sets all components of a vector to normally distributed random values.

   Logically Collective

   Input Parameters:
+  v    - the vector to be filled with random values
.  rctx - the random number context (can be NULL)
.  w1   - first work vector (can be NULL)
-  w2   - second work vector (can be NULL)

   Notes:
   Fills the two work vectors with uniformly distributed random values (VecSetRandom)
   and then applies the Box-Muller transform to get normally distributed values on v.

   Level: developer

.seealso: VecSetRandom()
@*/
PetscErrorCode VecSetRandomNormal(Vec v,PetscRandom rctx,Vec w1,Vec w2)
{
  const PetscScalar *x,*y;
  PetscScalar       *z;
  PetscInt          n,i;
  PetscRandom       rand=NULL;
  Vec               v1=NULL,v2=NULL;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  PetscValidType(v,1);
  if (rctx) PetscValidHeaderSpecific(rctx,PETSC_RANDOM_CLASSID,2);
  if (w1) PetscValidHeaderSpecific(w1,VEC_CLASSID,3);
  if (w2) PetscValidHeaderSpecific(w2,VEC_CLASSID,4);

  if (!rctx) {
    PetscCall(PetscRandomCreate(PetscObjectComm((PetscObject)v),&rand));
    PetscCall(PetscRandomSetFromOptions(rand));
    rctx = rand;
  }
  if (!w1) {
    PetscCall(VecDuplicate(v,&v1));
    w1 = v1;
  }
  if (!w2) {
    PetscCall(VecDuplicate(v,&v2));
    w2 = v2;
  }
  PetscCheckSameTypeAndComm(v,1,w1,3);
  PetscCheckSameTypeAndComm(v,1,w2,4);

  PetscCall(VecSetRandom(w1,rctx));
  PetscCall(VecSetRandom(w2,rctx));
  PetscCall(VecGetLocalSize(v,&n));
  PetscCall(VecGetArrayWrite(v,&z));
  PetscCall(VecGetArrayRead(w1,&x));
  PetscCall(VecGetArrayRead(w2,&y));
  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    z[i] = PetscCMPLX(PetscSqrtReal(-2.0*PetscLogReal(PetscRealPart(x[i])))*PetscCosReal(2.0*PETSC_PI*PetscRealPart(y[i])),PetscSqrtReal(-2.0*PetscLogReal(PetscImaginaryPart(x[i])))*PetscCosReal(2.0*PETSC_PI*PetscImaginaryPart(y[i])));
#else
    z[i] = PetscSqrtReal(-2.0*PetscLogReal(x[i]))*PetscCosReal(2.0*PETSC_PI*y[i]);
#endif
  }
  PetscCall(VecRestoreArrayWrite(v,&z));
  PetscCall(VecRestoreArrayRead(w1,&x));
  PetscCall(VecRestoreArrayRead(w2,&y));

  PetscCall(VecDestroy(&v1));
  PetscCall(VecDestroy(&v2));
  if (!rctx) PetscCall(PetscRandomDestroy(&rand));
  PetscFunctionReturn(PETSC_SUCCESS);
}
