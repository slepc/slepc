/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV routines related to Krylov decompositions
*/

#include <slepc/private/bvimpl.h>          /*I   "slepcbv.h"   I*/

/*@
   BVMatArnoldi - Computes an Arnoldi factorization associated with a matrix.

   Collective

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
.  H - (optional) the upper Hessenberg matrix
.  k - number of locked columns
-  m - dimension of the Arnoldi basis, may be modified

   Output Parameters:
+  beta - (optional) norm of last vector before normalization
-  breakdown - (optional) flag indicating that breakdown occurred

   Notes:
   Computes an `m`-step Arnoldi factorization for matrix `A`. The first `k` columns
   are assumed to be locked and therefore they are not modified. On exit, the
   following relation is satisfied

                  $$ A V - V H = \beta v_m e_m^* $$

   where the columns of $V$ are the Arnoldi vectors (which are orthonormal), $H$ is
   an upper Hessenberg matrix, $e_m$ is the $m$-th vector of the canonical basis.
   On exit, `beta` contains the norm of $v_m$ before normalization.

   The breakdown flag indicates that orthogonalization failed, see
   `BVOrthonormalizeColumn()`. In that case, on exit `m` contains the index of
   the column that failed.

   The values of `k` and `m` are not restricted to the active columns of `V`.

   To create an Arnoldi factorization from scratch, set `k=0` and make sure the
   first column contains the normalized initial vector.

   Level: advanced

.seealso: [](sec:bv), `BVMatLanczos()`, `BVSetActiveColumns()`, `BVOrthonormalizeColumn()`
@*/
PetscErrorCode BVMatArnoldi(BV V,Mat A,Mat H,PetscInt k,PetscInt *m,PetscReal *beta,PetscBool *breakdown)
{
  PetscScalar       *h;
  const PetscScalar *a;
  PetscInt          j,ldh,rows,cols;
  PetscBool         lindep=PETSC_FALSE;
  Vec               buf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,k,4);
  PetscAssertPointer(m,5);
  PetscValidLogicalCollectiveInt(V,*m,5);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(A,2);
  PetscCheckSameComm(V,1,A,2);

  PetscCheck(k>=0 && k<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument k has wrong value %" PetscInt_FMT ", should be between 0 and %" PetscInt_FMT,k,V->m);
  PetscCheck(*m>0 && *m<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument m has wrong value %" PetscInt_FMT ", should be between 1 and %" PetscInt_FMT,*m,V->m);
  PetscCheck(*m>k,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument m should be at least equal to k+1");
  if (H) {
    PetscValidHeaderSpecific(H,MAT_CLASSID,3);
    PetscValidType(H,3);
    PetscCheckTypeName(H,MATSEQDENSE);
    PetscCall(MatGetSize(H,&rows,&cols));
    PetscCall(MatDenseGetLDA(H,&ldh));
    PetscCheck(rows>=*m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Matrix H has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,rows,*m);
    PetscCheck(cols>=*m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Matrix H has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,cols,*m);
  }

  for (j=k;j<*m;j++) {
    PetscCall(BVMatMultColumn(V,A,j));
    if (PetscUnlikely(j==V->N-1)) PetscCall(BV_OrthogonalizeColumn_Safe(V,j+1,NULL,beta,&lindep)); /* safeguard in case the full basis is requested */
    else PetscCall(BVOrthonormalizeColumn(V,j+1,PETSC_FALSE,beta,&lindep));
    if (PetscUnlikely(lindep)) {
      *m = j+1;
      break;
    }
  }
  if (breakdown) *breakdown = lindep;
  if (lindep) PetscCall(PetscInfo(V,"Arnoldi finished early at m=%" PetscInt_FMT "\n",*m));

  if (H) {
    PetscCall(MatDenseGetArray(H,&h));
    PetscCall(BVGetBufferVec(V,&buf));
    PetscCall(VecGetArrayRead(buf,&a));
    for (j=k;j<*m-1;j++) PetscCall(PetscArraycpy(h+j*ldh,a+V->nc+(j+1)*(V->nc+V->m),j+2));
    PetscCall(PetscArraycpy(h+(*m-1)*ldh,a+V->nc+(*m)*(V->nc+V->m),*m));
    if (ldh>*m) h[(*m)+(*m-1)*ldh] = a[V->nc+(*m)+(*m)*(V->nc+V->m)];
    PetscCall(VecRestoreArrayRead(buf,&a));
    PetscCall(MatDenseRestoreArray(H,&h));
  }

  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   BVMatLanczos - Computes a Lanczos factorization associated with a matrix.

   Collective

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
.  T - (optional) the tridiagonal matrix
.  k - number of locked columns
-  m - dimension of the Lanczos basis, may be modified

   Output Parameters:
+  beta - (optional) norm of last vector before normalization
-  breakdown - (optional) flag indicating that breakdown occurred

   Notes:
   Computes an `m`-step Lanczos factorization for matrix `A`, with full
   reorthogonalization. At each Lanczos step, the corresponding Lanczos
   vector is orthogonalized with respect to all previous Lanczos vectors.
   This is equivalent to computing an `m`-step Arnoldi factorization and
   exploting symmetry of the operator.

   The first `k` columns are assumed to be locked and therefore they are
   not modified. On exit, the following relation is satisfied

                  $$ A V - V T = \beta v_m e_m^* $$

   where the columns of $V$ are the Lanczos vectors (which are $B$-orthonormal),
   $T$ is a real symmetric tridiagonal matrix, and $e_m$ is the $m$-th vector of
   the canonical basis. On exit, `beta` contains the $B$-norm of $v_m$ before
   normalization. The $T$ matrix is stored in a special way, its first column
   contains the diagonal elements, and its second column the off-diagonal
   ones. In complex scalars, the elements are stored as `PetscReal` and thus
   occupy only the first column of the Mat object. This is the same storage
   scheme used in matrix `DS_MAT_T` obtained with `DSGetMat()`.

   The breakdown flag indicates that orthogonalization failed, see
   `BVOrthonormalizeColumn()`. In that case, on exit `m` contains the index of
   the column that failed.

   The values of `k` and `m` are not restricted to the active columns of `V`.

   To create a Lanczos factorization from scratch, set `k=0` and make sure the
   first column contains the normalized initial vector.

   Level: advanced

.seealso: [](sec:bv), `BVMatArnoldi()`, `BVSetActiveColumns()`, `BVOrthonormalizeColumn()`, `DSGetMat()`
@*/
PetscErrorCode BVMatLanczos(BV V,Mat A,Mat T,PetscInt k,PetscInt *m,PetscReal *beta,PetscBool *breakdown)
{
  PetscScalar       *t;
  const PetscScalar *a;
  PetscReal         *alpha,*betat;
  PetscInt          j,ldt,rows,cols,mincols=PetscDefined(USE_COMPLEX)?1:2;
  PetscBool         lindep=PETSC_FALSE;
  Vec               buf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,k,4);
  PetscAssertPointer(m,5);
  PetscValidLogicalCollectiveInt(V,*m,5);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(A,2);
  PetscCheckSameComm(V,1,A,2);

  PetscCheck(k>=0 || k<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument k has wrong value %" PetscInt_FMT ", should be between 0 and %" PetscInt_FMT,k,V->m);
  PetscCheck(*m>0 || *m<=V->m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument m has wrong value %" PetscInt_FMT ", should be between 1 and %" PetscInt_FMT,*m,V->m);
  PetscCheck(*m>k,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument m should be at least equal to k+1");
  if (T) {
    PetscValidHeaderSpecific(T,MAT_CLASSID,3);
    PetscValidType(T,3);
    PetscCheckTypeName(T,MATSEQDENSE);
    PetscCall(MatGetSize(T,&rows,&cols));
    PetscCall(MatDenseGetLDA(T,&ldt));
    PetscCheck(rows>=*m,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Matrix T has %" PetscInt_FMT " rows, should have at least %" PetscInt_FMT,rows,*m);
    PetscCheck(cols>=mincols,PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Matrix T has %" PetscInt_FMT " columns, should have at least %" PetscInt_FMT,cols,mincols);
  }

  for (j=k;j<*m;j++) {
    PetscCall(BVMatMultColumn(V,A,j));
    if (PetscUnlikely(j==V->N-1)) PetscCall(BV_OrthogonalizeColumn_Safe(V,j+1,NULL,beta,&lindep)); /* safeguard in case the full basis is requested */
    else PetscCall(BVOrthonormalizeColumn(V,j+1,PETSC_FALSE,beta,&lindep));
    if (PetscUnlikely(lindep)) {
      *m = j+1;
      break;
    }
  }
  if (breakdown) *breakdown = lindep;
  if (lindep) PetscCall(PetscInfo(V,"Lanczos finished early at m=%" PetscInt_FMT "\n",*m));

  if (T) {
    PetscCall(MatDenseGetArray(T,&t));
    alpha = (PetscReal*)t;
    betat = alpha+ldt;
    PetscCall(BVGetBufferVec(V,&buf));
    PetscCall(VecGetArrayRead(buf,&a));
    for (j=k;j<*m;j++) {
      alpha[j] = PetscRealPart(a[V->nc+j+(j+1)*(V->nc+V->m)]);
      betat[j] = PetscRealPart(a[V->nc+j+1+(j+1)*(V->nc+V->m)]);
    }
    PetscCall(VecRestoreArrayRead(buf,&a));
    PetscCall(MatDenseRestoreArray(T,&t));
  }

  PetscCall(PetscObjectStateIncrease((PetscObject)V));
  PetscFunctionReturn(PETSC_SUCCESS);
}
