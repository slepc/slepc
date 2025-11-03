/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>            /*I "slepcsys.h" I*/

static PetscErrorCode MatCreateTile_Seq(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat G)
{
  PetscInt          i,j,M1,M2,N1,N2,ncols,*scols;
  PetscScalar       *svals,*buf;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,&M1,&N1));
  PetscCall(MatGetSize(D,&M2,&N2));

  PetscCall(PetscMalloc2(PetscMax(N1,N2),&buf,PetscMax(N1,N2),&scols));
  /* Transfer A */
  if (a!=0.0) {
    for (i=0;i<M1;i++) {
      PetscCall(MatGetRow(A,i,&ncols,&cols,&vals));
      if (a!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*a;
      } else svals=(PetscScalar*)vals;
      PetscCall(MatSetValues(G,1,&i,ncols,cols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(A,i,&ncols,&cols,&vals));
    }
  }
  /* Transfer B */
  if (b!=0.0) {
    for (i=0;i<M1;i++) {
      PetscCall(MatGetRow(B,i,&ncols,&cols,&vals));
      if (b!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*b;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = cols[j]+N1;
      PetscCall(MatSetValues(G,1,&i,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(B,i,&ncols,&cols,&vals));
    }
  }
  /* Transfer C */
  if (c!=0.0) {
    for (i=0;i<M2;i++) {
      PetscCall(MatGetRow(C,i,&ncols,&cols,&vals));
      if (c!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*c;
      } else svals=(PetscScalar*)vals;
      j = i+M1;
      PetscCall(MatSetValues(G,1,&j,ncols,cols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(C,i,&ncols,&cols,&vals));
    }
  }
  /* Transfer D */
  if (d!=0.0) {
    for (i=0;i<M2;i++) {
      PetscCall(MatGetRow(D,i,&ncols,&cols,&vals));
      if (d!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*d;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = cols[j]+N1;
      j = i+M1;
      PetscCall(MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(D,i,&ncols,&cols,&vals));
    }
  }
  PetscCall(PetscFree2(buf,scols));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MatCreateTile_MPI(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat G)
{
  PetscMPIInt       np;
  PetscInt          p,i,j,N1,N2,m1,m2,*map1,*map2;
  PetscInt          ncols,*scols,start,gstart;
  PetscScalar       *svals,*buf;
  const PetscInt    *cols,*mapptr1,*mapptr2;
  const PetscScalar *vals;

  PetscFunctionBegin;
  PetscCall(MatGetSize(A,NULL,&N1));
  PetscCall(MatGetLocalSize(A,&m1,NULL));
  PetscCall(MatGetSize(D,NULL,&N2));
  PetscCall(MatGetLocalSize(D,&m2,NULL));

  /* Create mappings */
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)G),&np));
  PetscCall(MatGetOwnershipRangesColumn(A,&mapptr1));
  PetscCall(MatGetOwnershipRangesColumn(B,&mapptr2));
  PetscCall(PetscMalloc4(PetscMax(N1,N2),&buf,PetscMax(N1,N2),&scols,N1,&map1,N2,&map2));
  for (p=0;p<np;p++) {
    for (i=mapptr1[p];i<mapptr1[p+1];i++) map1[i] = i+mapptr2[p];
  }
  for (p=0;p<np;p++) {
    for (i=mapptr2[p];i<mapptr2[p+1];i++) map2[i] = i+mapptr1[p+1];
  }
  PetscCall(MatGetOwnershipRange(G,&gstart,NULL));

  /* Transfer A */
  if (a!=0.0) {
    PetscCall(MatGetOwnershipRange(A,&start,NULL));
    for (i=0;i<m1;i++) {
      PetscCall(MatGetRow(A,i+start,&ncols,&cols,&vals));
      if (a!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*a;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      j = gstart+i;
      PetscCall(MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(A,i+start,&ncols,&cols,&vals));
    }
  }
  /* Transfer B */
  if (b!=0.0) {
    PetscCall(MatGetOwnershipRange(B,&start,NULL));
    for (i=0;i<m1;i++) {
      PetscCall(MatGetRow(B,i+start,&ncols,&cols,&vals));
      if (b!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*b;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      j = gstart+i;
      PetscCall(MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(B,i+start,&ncols,&cols,&vals));
    }
  }
  /* Transfer C */
  if (c!=0.0) {
    PetscCall(MatGetOwnershipRange(C,&start,NULL));
    for (i=0;i<m2;i++) {
      PetscCall(MatGetRow(C,i+start,&ncols,&cols,&vals));
      if (c!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*c;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      j = gstart+m1+i;
      PetscCall(MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(C,i+start,&ncols,&cols,&vals));
    }
  }
  /* Transfer D */
  if (d!=0.0) {
    PetscCall(MatGetOwnershipRange(D,&start,NULL));
    for (i=0;i<m2;i++) {
      PetscCall(MatGetRow(D,i+start,&ncols,&cols,&vals));
      if (d!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*d;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      j = gstart+m1+i;
      PetscCall(MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES));
      PetscCall(MatRestoreRow(D,i+start,&ncols,&cols,&vals));
    }
  }
  PetscCall(PetscFree4(buf,scols,map1,map2));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatCreateTile - Explicitly build a matrix from four blocks.

   Collective

   Input Parameters:
+  A - matrix for top-left block
.  a - scaling factor for block A
.  B - matrix for top-right block
.  b - scaling factor for block B
.  C - matrix for bottom-left block
.  c - scaling factor for block C
.  D - matrix for bottom-right block
-  d - scaling factor for block D

   Output Parameter:
.  G  - the resulting matrix

   Notes:
   The resulting matrix is built as
   $$G = \begin{bmatrix} aA & bB\\ cC & dD \end{bmatrix}.$$

   In the case of a parallel matrix, a permuted version of `G` is returned. The permutation
   is a perfect shuffle such that the local parts of `A`, `B`, `C`, `D` remain in the local part of
   `G` for the same process.

   Matrix `G` must be destroyed by the user.

   The blocks can be of different type. They can be either `MATCONSTANTDIAGONAL`, or a standard
   type such as `MATAIJ`, or any other type provided that it supports the `MatGetRow()` operation.
   The type of the output matrix will be the same as the first block that is not
   `MATCONSTANTDIAGONAL` (checked in the `A`,`B`,`C`,`D` order).

   Level: developer

.seealso: `MatCreateNest()`
@*/
PetscErrorCode MatCreateTile(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat *G)
{
  PetscInt       i,k,M1,M2,N1,N2,M,N,m1,m2,n1,n2,m,n,bs;
  PetscBool      diag[4];
  Mat            block[4] = {A,B,C,D};
  MatType        type[4];
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(B,MAT_CLASSID,4);
  PetscValidHeaderSpecific(C,MAT_CLASSID,6);
  PetscValidHeaderSpecific(D,MAT_CLASSID,8);
  PetscCheckSameTypeAndComm(A,2,B,4);
  PetscCheckSameTypeAndComm(A,2,C,6);
  PetscCheckSameTypeAndComm(A,2,D,8);
  PetscValidLogicalCollectiveScalar(A,a,1);
  PetscValidLogicalCollectiveScalar(A,b,3);
  PetscValidLogicalCollectiveScalar(A,c,5);
  PetscValidLogicalCollectiveScalar(A,d,7);
  PetscAssertPointer(G,9);

  /* check row 1 */
  PetscCall(MatGetSize(A,&M1,NULL));
  PetscCall(MatGetLocalSize(A,&m1,NULL));
  PetscCall(MatGetSize(B,&M,NULL));
  PetscCall(MatGetLocalSize(B,&m,NULL));
  PetscCheck(M==M1 && m==m1,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check row 2 */
  PetscCall(MatGetSize(C,&M2,NULL));
  PetscCall(MatGetLocalSize(C,&m2,NULL));
  PetscCall(MatGetSize(D,&M,NULL));
  PetscCall(MatGetLocalSize(D,&m,NULL));
  PetscCheck(M==M2 && m==m2,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check column 1 */
  PetscCall(MatGetSize(A,NULL,&N1));
  PetscCall(MatGetLocalSize(A,NULL,&n1));
  PetscCall(MatGetSize(C,NULL,&N));
  PetscCall(MatGetLocalSize(C,NULL,&n));
  PetscCheck(N==N1 && n==n1,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check column 2 */
  PetscCall(MatGetSize(B,NULL,&N2));
  PetscCall(MatGetLocalSize(B,NULL,&n2));
  PetscCall(MatGetSize(D,NULL,&N));
  PetscCall(MatGetLocalSize(D,NULL,&n));
  PetscCheck(N==N2 && n==n2,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");

  /* check matrix types */
  for (i=0;i<4;i++) {
    PetscCall(MatGetType(block[i],&type[i]));
    PetscCall(PetscStrcmp(type[i],MATCONSTANTDIAGONAL,&diag[i]));
  }
  for (k=0;k<4;k++) if (!diag[k]) break;
  PetscCheck(k<4,PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Not implemented for 4 diagonal blocks");

  PetscCall(MatGetBlockSize(block[k],&bs));
  PetscCall(MatCreate(PetscObjectComm((PetscObject)block[k]),G));
  PetscCall(MatSetSizes(*G,m1+m2,n1+n2,M1+M2,N1+N2));
  PetscCall(MatSetType(*G,type[k]));
  PetscCall(MatSetBlockSize(*G,bs));
  PetscCall(MatSetUp(*G));

  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)*G),&size));
  if (size>1) PetscCall(MatCreateTile_MPI(a,A,b,B,c,C,d,D,*G));
  else PetscCall(MatCreateTile_Seq(a,A,b,B,c,C,d,D,*G));
  PetscCall(MatAssemblyBegin(*G,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*G,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatCreateVecsEmpty - Get vector(s) compatible with the matrix, i.e., with the same
   parallel layout, but without internal array.

   Collective

   Input Parameter:
.  mat - the matrix

   Output Parameters:
+  right - (optional) vector that the matrix can be multiplied against
-  left - (optional) vector that the matrix vector product can be stored in

   Note:
   This is similar to `MatCreateVecs()`, but the new vectors do not have an internal
   array, so the intended usage is with `VecPlaceArray()`.

   Level: developer

.seealso: `MatCreateVecs()`, `VecDuplicateEmpty()`
@*/
PetscErrorCode MatCreateVecsEmpty(Mat mat,Vec *right,Vec *left)
{
  PetscBool      standard,cuda=PETSC_FALSE,hip=PETSC_FALSE,skip=PETSC_FALSE;
  PetscInt       M,N,mloc,nloc,rbs,cbs;
  PetscMPIInt    size;
  Vec            v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mat,MAT_CLASSID,1);
  PetscValidType(mat,1);

  PetscCall(PetscObjectTypeCompareAny((PetscObject)mat,&standard,MATSEQAIJ,MATMPIAIJ,MATSEQBAIJ,MATMPIBAIJ,MATSEQSBAIJ,MATMPISBAIJ,MATSEQDENSE,MATMPIDENSE,""));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)mat,&cuda,MATSEQAIJCUSPARSE,MATMPIAIJCUSPARSE,""));
  PetscCall(PetscObjectTypeCompareAny((PetscObject)mat,&hip,MATSEQAIJHIPSPARSE,MATMPIAIJHIPSPARSE,""));
  if (!standard && !cuda && !hip) {
    PetscCall(MatCreateVecs(mat,right,left));
    v = right? *right: *left;
    if (v) {
      PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&standard,VECSEQ,VECMPI,""));
      PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&cuda,VECSEQCUDA,VECMPICUDA,""));
      PetscCall(PetscObjectTypeCompareAny((PetscObject)v,&hip,VECSEQHIP,VECMPIHIP,""));
    }
    if (!standard && !cuda && !hip) skip = PETSC_TRUE;
    else {
      if (right) PetscCall(VecDestroy(right));
      if (left) PetscCall(VecDestroy(left));
    }
  }
  if (!skip) {
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)mat),&size));
    PetscCall(MatGetLocalSize(mat,&mloc,&nloc));
    PetscCall(MatGetSize(mat,&M,&N));
    PetscCall(MatGetBlockSizes(mat,&rbs,&cbs));
    if (right) {
      if (cuda) {
#if defined(PETSC_HAVE_CUDA)
        if (size>1) PetscCall(VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)mat),cbs,nloc,N,NULL,right));
        else PetscCall(VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)mat),cbs,N,NULL,right));
#endif
      } else if (hip) {
#if defined(PETSC_HAVE_HIP)
        if (size>1) PetscCall(VecCreateMPIHIPWithArray(PetscObjectComm((PetscObject)mat),cbs,nloc,N,NULL,right));
        else PetscCall(VecCreateSeqHIPWithArray(PetscObjectComm((PetscObject)mat),cbs,N,NULL,right));
#endif
      } else {
        if (size>1) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)mat),cbs,nloc,N,NULL,right));
        else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)mat),cbs,N,NULL,right));
      }
    }
    if (left) {
      if (cuda) {
#if defined(PETSC_HAVE_CUDA)
        if (size>1) PetscCall(VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)mat),rbs,mloc,M,NULL,left));
        else PetscCall(VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)mat),rbs,M,NULL,left));
#endif
      } else if (hip) {
#if defined(PETSC_HAVE_HIP)
        if (size>1) PetscCall(VecCreateMPIHIPWithArray(PetscObjectComm((PetscObject)mat),rbs,mloc,M,NULL,left));
        else PetscCall(VecCreateSeqHIPWithArray(PetscObjectComm((PetscObject)mat),rbs,M,NULL,left));
#endif
      } else {
        if (size>1) PetscCall(VecCreateMPIWithArray(PetscObjectComm((PetscObject)mat),rbs,mloc,M,NULL,left));
        else PetscCall(VecCreateSeqWithArray(PetscObjectComm((PetscObject)mat),rbs,M,NULL,left));
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatNormEstimate - Estimate the 2-norm of a matrix.

   Collective

   Input Parameters:
+  A   - the matrix
.  vrn - random vector with normally distributed entries (can be `NULL`)
-  w   - workspace vector (can be `NULL`)

   Output Parameter:
.  nrm - the norm estimate

   Notes:
   Does not need access to the matrix entries, just performs a matrix-vector product.
   Based on work by I. Ipsen and coworkers <https://ipsen.math.ncsu.edu/ps/slides_ima.pdf>.

   The input vector `vrn` must have unit 2-norm.
   If `vrn` is `NULL`, then it is created internally and filled with `VecSetRandomNormal()`.

   Level: developer

.seealso: `VecSetRandomNormal()`
@*/
PetscErrorCode MatNormEstimate(Mat A,Vec vrn,Vec w,PetscReal *nrm)
{
  PetscInt       n;
  Vec            vv=NULL,ww=NULL;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_CLASSID,1);
  PetscValidType(A,1);
  if (vrn) PetscValidHeaderSpecific(vrn,VEC_CLASSID,2);
  if (w) PetscValidHeaderSpecific(w,VEC_CLASSID,3);
  PetscAssertPointer(nrm,4);

  if (!vrn) {
    PetscCall(MatCreateVecs(A,&vv,NULL));
    vrn = vv;
    PetscCall(VecSetRandomNormal(vv,NULL,NULL,NULL));
    PetscCall(VecNormalize(vv,NULL));
  }
  if (!w) {
    PetscCall(MatCreateVecs(A,&ww,NULL));
    w = ww;
  }

  PetscCall(MatGetSize(A,&n,NULL));
  PetscCall(MatMult(A,vrn,w));
  PetscCall(VecNorm(w,NORM_2,nrm));
  *nrm *= PetscSqrtReal((PetscReal)n);

  PetscCall(VecDestroy(&vv));
  PetscCall(VecDestroy(&ww));
  PetscFunctionReturn(PETSC_SUCCESS);
}
