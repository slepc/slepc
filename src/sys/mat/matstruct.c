/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>            /*I "slepcsys.h" I*/

/*@
   MatCreateBSE - Create a matrix that can be used to define a structured eigenvalue
   problem of type BSE (Bethe-Salpeter Equation).

   Collective

   Input Parameters:
+  R - matrix for the diagonal block (resonant)
-  C - matrix for the off-diagonal block (coupling)

   Output Parameter:
.  H  - the resulting matrix

   Notes:
   The resulting matrix has the block form $$H = \begin{bmatrix} R & C \\ -C^* & -R^T \end{bmatrix},$$
   where $R$ is assumed
   to be (complex) Hermitian and $C$ complex symmetric. Note that this function does
   not check these properties, so if the matrices provided by the user do not satisfy
   them, then the solver will not behave as expected.

   The obtained matrix can be used as an input matrix to `EPS` eigensolvers via
   `EPSSetOperators()` for the case that the problem type is `EPS_BSE`. Note that the user
   cannot just build a matrix with the required structure, it must be done via this
   function.

   In the current implementation, `H` is a `MATNEST` matrix, where `R` and `C` form the top
   block row, while the bottom block row is composed of matrices of type
   `MATTRANSPOSEVIRTUAL` and `MATHERMITIANTRANSPOSEVIRTUAL` scaled by -1.

   Level: intermediate

.seealso: [](sec:structured), `MatCreateNest()`, `EPSSetOperators()`, `EPSSetProblemType()`, `MatCreateHamiltonian()`, `MatCreateLREP()`
@*/
PetscErrorCode MatCreateBSE(Mat R,Mat C,Mat *H)
{
  PetscInt       Mr,Mc,Nr,Nc,mr,mc,nr,nc;
  Mat            block[4] = { R, C, NULL, NULL };
  SlepcMatStruct mctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(R,MAT_CLASSID,1);
  PetscValidHeaderSpecific(C,MAT_CLASSID,2);
  PetscCheckSameTypeAndComm(R,1,C,2);
  PetscAssertPointer(H,3);

  /* check sizes */
  PetscCall(MatGetSize(R,&Mr,&Nr));
  PetscCall(MatGetLocalSize(R,&mr,&nr));
  PetscCall(MatGetSize(C,&Mc,&Nc));
  PetscCall(MatGetLocalSize(C,&mc,&nc));
  PetscCheck(Mc==Mr && mc==mr,PetscObjectComm((PetscObject)R),PETSC_ERR_ARG_INCOMP,"Incompatible row dimensions");
  PetscCheck(Nc==Nr && nc==nr,PetscObjectComm((PetscObject)R),PETSC_ERR_ARG_INCOMP,"Incompatible column dimensions");

  /* bottom block row */
  PetscCall(MatCreateHermitianTranspose(C,&block[2]));
  PetscCall(MatScale(block[2],-1.0));
  PetscCall(MatCreateTranspose(R,&block[3]));
  PetscCall(MatScale(block[3],-1.0));

  /* create nest matrix and compose context */
  PetscCall(MatCreateNest(PetscObjectComm((PetscObject)R),2,NULL,2,NULL,block,H));
  PetscCall(MatDestroy(&block[2]));
  PetscCall(MatDestroy(&block[3]));

  PetscCall(PetscNew(&mctx));
  mctx->cookie = SLEPC_MAT_STRUCT_BSE;
  PetscCall(PetscObjectContainerCompose((PetscObject)*H,"SlepcMatStruct",mctx,PetscCtxDestroyDefault));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatCreateHamiltonian - Create a matrix that can be used to define a structured
   eigenvalue problem of Hamiltonian type.

   Collective

   Input Parameters:
+  A - matrix for the (0,0) block
.  B - matrix for the (0,1) block, must be real symmetric or complex Hermitian
-  C - matrix for the (1,0) block, must be real symmetric or complex Hermitian

   Output Parameter:
.  H  - the resulting matrix

   Notes:
   The resulting matrix has the block form $$H = \begin{bmatrix} A & B \\ C & -A^* \end{bmatrix},$$
   where $B$ and $C$ are
   assumed to be symmetric in the real case or Hermitian in the complex case. Note
   that this function does not check this property, so if the matrices provided by
   the user do not satisfy it, then the solver will not behave as expected.

   The obtained matrix can be used as an input matrix to `EPS` eigensolvers via
   `EPSSetOperators()` for the case that the problem type is `EPS_HAMILT`. Note that the
   user cannot just build a matrix with the required structure, it must be done via
   this function.

   In the current implementation, `H` is a `MATNEST` matrix, where the (1,1) block is
   a matrix of type `MATHERMITIANTRANSPOSEVIRTUAL` obtained from `A` and scaled by -1.

   Level: intermediate

.seealso: [](sec:structured), `MatCreateNest()`, `EPSSetOperators()`, `EPSSetProblemType()`, `MatCreateBSE()`, `MatCreateLREP()`
@*/
PetscErrorCode MatCreateHamiltonian(Mat A,Mat B,Mat C,Mat *H)
{
  PetscInt       Ma,Mb,Mc,Na,Nb,Nc,ma,mb,mc,na,nb,nc;
  Mat            block[4] = { A, B, C, NULL };
  SlepcMatStruct mctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_CLASSID,1);
  PetscValidHeaderSpecific(B,MAT_CLASSID,2);
  PetscValidHeaderSpecific(C,MAT_CLASSID,3);
  PetscCheckSameTypeAndComm(A,1,B,2);
  PetscCheckSameTypeAndComm(A,1,C,3);
  PetscAssertPointer(H,4);

  /* check sizes */
  PetscCall(MatGetSize(A,&Ma,&Na));
  PetscCall(MatGetLocalSize(A,&ma,&na));
  PetscCall(MatGetSize(B,&Mb,&Nb));
  PetscCall(MatGetLocalSize(B,&mb,&nb));
  PetscCall(MatGetSize(C,&Mc,&Nc));
  PetscCall(MatGetLocalSize(C,&mc,&nc));
  PetscCheck(Mb==Ma && mb==ma,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible row dimensions");
  PetscCheck(Nc==Na && nc==na,PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible column dimensions");

  /* (1,1) block */
  PetscCall(MatCreateHermitianTranspose(A,&block[3]));
  PetscCall(MatScale(block[3],-1.0));

  /* create nest matrix and compose context */
  PetscCall(MatCreateNest(PetscObjectComm((PetscObject)A),2,NULL,2,NULL,block,H));
  PetscCall(MatDestroy(&block[3]));

  PetscCall(PetscNew(&mctx));
  mctx->cookie = SLEPC_MAT_STRUCT_HAMILT;
  PetscCall(PetscObjectContainerCompose((PetscObject)*H,"SlepcMatStruct",mctx,PetscCtxDestroyDefault));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@
   MatCreateLREP - Create a matrix that can be used to define a structured Linear
   Response eigenvalue problem.

   Collective

   Input Parameters:
+  AK  - matrix for the diagonal block ($A$) or the top block ($K$)
.  BM  - matrix for the off-diagonal block ($B$) or the bottom block ($M$)
-  red - whether the reduced form should be built

   Output Parameter:
.  H  - the resulting matrix

   Notes:
   The resulting matrix has the block form $$H = \begin{bmatrix} A & B \\ -B & -A \end{bmatrix},$$
   where both $A$ and $B$ are assumed to be real symmetric. An alternative form,
   called reduced form (`red`=`PETSC_TRUE`) is $$H = \begin{bmatrix} 0 & K \\ M & 0 \end{bmatrix},$$
   where both $K$ and $M$ are also assumed to be real symmetric. Note that this function
   does not check these properties, so if the matrices provided by the user do not
   satisfy them, then the solver will not behave as expected.

   The reduced form can be obtained from the original formulation by setting $K=A-B$
   and $M=A+B$. However, $K$ and $M$ need not have such relation, in which case the
   computed eigenvalues will be $\pm\lambda_i$, where $\lambda_i^2$ are the eigenvalues
   of the _product eigenvalue problem_ $KM$ (or $MK$).

   The obtained matrix can be used as an input matrix to `EPS` eigensolvers via
   `EPSSetOperators()` for the case that the problem type is `EPS_LREP`. Note that the
   user cannot just build a matrix with the required structure, it must be done via
   this function.

   In the current implementation, `H` is a `MATNEST` matrix, where `A` and `B` form the
   top block row, while the bottom block row is composed of shell matrices scaled by -1.
   Similarly, in the reduced formulation, the `MATNEST` has two zero blocks.

   Internally, the eigensolver will work with the reduced form, and it will carry out
   the conversion transparently if the first form is provided by the user.

   Level: intermediate

.seealso: [](sec:structured), `MatCreateNest()`, `EPSSetOperators()`, `EPSSetProblemType()`, `MatCreateBSE()`, `MatCreateHamiltonian()`
@*/
PetscErrorCode MatCreateLREP(Mat AK,Mat BM,PetscBool red,Mat *H)
{
  PetscInt       Ma,Mb,Na,Nb,ma,mb,na,nb;
  SlepcMatStruct mbtx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(AK,MAT_CLASSID,1);
  PetscValidHeaderSpecific(BM,MAT_CLASSID,2);
  PetscCheckSameTypeAndComm(AK,1,BM,2);
  PetscValidLogicalCollectiveBool(AK,red,3);
  PetscAssertPointer(H,4);

  /* check sizes */
  PetscCall(MatGetSize(AK,&Ma,&Na));
  PetscCall(MatGetLocalSize(AK,&ma,&na));
  PetscCall(MatGetSize(BM,&Mb,&Nb));
  PetscCall(MatGetLocalSize(BM,&mb,&nb));
  PetscCheck(Mb==Ma && mb==ma,PetscObjectComm((PetscObject)AK),PETSC_ERR_ARG_INCOMP,"Incompatible row dimensions");
  PetscCheck(Nb==Na && nb==na,PetscObjectComm((PetscObject)AK),PETSC_ERR_ARG_INCOMP,"Incompatible column dimensions");

  if (red) {  /* reduced form */

    Mat block[4] = { NULL, AK, BM, NULL };

    /* create nest matrix */
    PetscCall(MatCreateNest(PetscObjectComm((PetscObject)AK),2,NULL,2,NULL,block,H));

  } else {  /* original LREP form */

    Mat block[4] = { AK, BM, NULL, NULL };
    const PetscScalar scal[] = { -1.0 };

    /* bottom block row, using MATCOMPOSITE */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)AK),&block[2]));
    PetscCall(MatSetSizes(block[2],mb,nb,Mb,Nb));
    PetscCall(MatSetType(block[2],MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(block[2],BM));
    PetscCall(MatAssemblyBegin(block[2],MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(block[2],MAT_FINAL_ASSEMBLY));
    PetscCall(MatCompositeSetScalings(block[2],scal));

    PetscCall(MatCreate(PetscObjectComm((PetscObject)AK),&block[3]));
    PetscCall(MatSetSizes(block[3],ma,na,Ma,Na));
    PetscCall(MatSetType(block[3],MATCOMPOSITE));
    PetscCall(MatCompositeAddMat(block[3],AK));
    PetscCall(MatAssemblyBegin(block[3],MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(block[3],MAT_FINAL_ASSEMBLY));
    PetscCall(MatCompositeSetScalings(block[3],scal));

    /* create nest matrix */
    PetscCall(MatCreateNest(PetscObjectComm((PetscObject)AK),2,NULL,2,NULL,block,H));
    PetscCall(MatDestroy(&block[2]));
    PetscCall(MatDestroy(&block[3]));

  }

  /* compose context */
  PetscCall(PetscNew(&mbtx));
  mbtx->cookie = SLEPC_MAT_STRUCT_LREP;
  PetscCall(PetscObjectContainerCompose((PetscObject)*H,"SlepcMatStruct",mbtx,PetscCtxDestroyDefault));
  PetscFunctionReturn(PETSC_SUCCESS);
}
