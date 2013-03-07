/*
     SVD routines for setting up the solver.

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2012, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
      
   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY 
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS 
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for 
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc-private/svdimpl.h>      /*I "slepcsvd.h" I*/
#include <slepc-private/ipimpl.h>

#undef __FUNCT__
#define __FUNCT__ "SVDSetOperator"
/*@
   SVDSetOperator - Set the matrix associated with the singular value problem.

   Collective on SVD and Mat

   Input Parameters:
+  svd - the singular value solver context
-  A  - the matrix associated with the singular value problem

   Level: beginner

.seealso: SVDSolve(), SVDGetOperator()
@*/
PetscErrorCode SVDSetOperator(SVD svd,Mat mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(mat,MAT_CLASSID,2);
  PetscCheckSameComm(svd,1,mat,2);
  if (svd->setupcalled) { ierr = SVDReset(svd);CHKERRQ(ierr); }
  ierr = PetscObjectReference((PetscObject)mat);CHKERRQ(ierr);
  ierr = MatDestroy(&svd->OP);CHKERRQ(ierr);
  svd->OP = mat;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SVDGetOperator"
/*@
   SVDGetOperator - Get the matrix associated with the singular value problem.

   Not collective, though parallel Mats are returned if the SVD is parallel

   Input Parameter:
.  svd - the singular value solver context

   Output Parameters:
.  A    - the matrix associated with the singular value problem

   Level: advanced

.seealso: SVDSolve(), SVDSetOperator()
@*/
PetscErrorCode SVDGetOperator(SVD svd,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(A,2);
  *A = svd->OP;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SVDSetUp"
/*@
   SVDSetUp - Sets up all the internal data structures necessary for the
   execution of the singular value solver.

   Collective on SVD

   Input Parameter:
.  svd   - singular value solver context

   Level: advanced

   Notes:
   This function need not be called explicitly in most cases, since SVDSolve()
   calls it. It can be useful when one wants to measure the set-up time 
   separately from the solve time.

.seealso: SVDCreate(), SVDSolve(), SVDDestroy()
@*/
PetscErrorCode SVDSetUp(SVD svd)
{
  PetscErrorCode ierr;
  PetscBool      flg;
  PetscInt       M,N,k;
  Vec            *T;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->setupcalled) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(SVD_SetUp,svd,0,0,0);CHKERRQ(ierr);

  /* reset the convergence flag from the previous solves */
  svd->reason = SVD_CONVERGED_ITERATING;

  /* Set default solver type (SVDSetFromOptions was not called) */
  if (!((PetscObject)svd)->type_name) {
    ierr = SVDSetType(svd,SVDCROSS);CHKERRQ(ierr);
  }
  if (!svd->ip) { ierr = SVDGetIP(svd,&svd->ip);CHKERRQ(ierr); }
  if (!((PetscObject)svd->ip)->type_name) {
    ierr = IPSetType_Default(svd->ip);CHKERRQ(ierr);
  }
  if (!svd->ds) { ierr = SVDGetDS(svd,&svd->ds);CHKERRQ(ierr); }
  ierr = DSReset(svd->ds);CHKERRQ(ierr);
  if (!((PetscObject)svd->rand)->type_name) {
    ierr = PetscRandomSetFromOptions(svd->rand);CHKERRQ(ierr);
  }

  /* check matrix */
  if (!svd->OP) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_WRONGSTATE,"SVDSetOperator must be called first"); 
  
  /* determine how to build the transpose */
  if (svd->transmode == PETSC_DECIDE) {
    ierr = MatHasOperation(svd->OP,MATOP_TRANSPOSE,&flg);CHKERRQ(ierr);    
    if (flg) svd->transmode = SVD_TRANSPOSE_EXPLICIT;
    else svd->transmode = SVD_TRANSPOSE_IMPLICIT;
  }
  
  /* build transpose matrix */
  ierr = MatDestroy(&svd->A);CHKERRQ(ierr);
  ierr = MatDestroy(&svd->AT);CHKERRQ(ierr);
  ierr = MatGetSize(svd->OP,&M,&N);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)svd->OP);CHKERRQ(ierr);
  switch (svd->transmode) {
    case SVD_TRANSPOSE_EXPLICIT:
      if (M>=N) {
        svd->A = svd->OP;
        ierr = MatTranspose(svd->OP,MAT_INITIAL_MATRIX,&svd->AT);CHKERRQ(ierr);
        ierr = MatConjugate(svd->AT);CHKERRQ(ierr);
      } else {
        ierr = MatTranspose(svd->OP,MAT_INITIAL_MATRIX,&svd->A);CHKERRQ(ierr);
        ierr = MatConjugate(svd->A);CHKERRQ(ierr);
        svd->AT = svd->OP;
      }
      break;
    case SVD_TRANSPOSE_IMPLICIT:
      if (M>=N) {
        svd->A = svd->OP;
        svd->AT = NULL;    
      } else {
        svd->A = NULL;
        svd->AT = svd->OP;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid transpose mode"); 
  }

  ierr = VecDestroy(&svd->tr);CHKERRQ(ierr);
  ierr = VecDestroy(&svd->tl);CHKERRQ(ierr);
  if (svd->A) {
    ierr = SlepcMatGetVecsTemplate(svd->A,&svd->tr,&svd->tl);CHKERRQ(ierr);
  } else {
    ierr = SlepcMatGetVecsTemplate(svd->AT,&svd->tl,&svd->tr);CHKERRQ(ierr);
  }

  /* swap initial vectors if necessary */
  if (M<N) {
    T=svd->ISL; svd->ISL=svd->IS; svd->IS=T;
    k=svd->ninil; svd->ninil=svd->nini; svd->nini=k;
  }

  /* call specific solver setup */
  ierr = (*svd->ops->setup)(svd);CHKERRQ(ierr);

  /* set tolerance if not yet set */
  if (svd->tol==PETSC_DEFAULT) svd->tol = SLEPC_DEFAULT_TOL;

  if (svd->ncv > M || svd->ncv > N) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"ncv bigger than matrix dimensions");
  if (svd->nsv > svd->ncv) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"nsv bigger than ncv");

  if (svd->ncv != svd->n) {
    /* free memory for previous solution  */
    if (svd->n) { 
      ierr = PetscFree(svd->sigma);CHKERRQ(ierr);
      ierr = PetscFree(svd->perm);CHKERRQ(ierr);
      ierr = PetscFree(svd->errest);CHKERRQ(ierr);
      ierr = VecDestroyVecs(svd->n,&svd->V);CHKERRQ(ierr);
    }
    /* allocate memory for next solution */
    ierr = PetscMalloc(svd->ncv*sizeof(PetscReal),&svd->sigma);CHKERRQ(ierr);
    ierr = PetscMalloc(svd->ncv*sizeof(PetscInt),&svd->perm);CHKERRQ(ierr);
    ierr = PetscMalloc(svd->ncv*sizeof(PetscReal),&svd->errest);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory(svd,PetscMax(0,svd->ncv-svd->n)*(2*sizeof(PetscReal)+sizeof(PetscInt)));CHKERRQ(ierr);
    ierr = VecDuplicateVecs(svd->tr,svd->ncv,&svd->V);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(svd,svd->ncv,svd->V);CHKERRQ(ierr);
    svd->n = svd->ncv;
  }

  /* process initial vectors */
  if (svd->nini<0) {
    svd->nini = -svd->nini;
    if (svd->nini>svd->ncv) SETERRQ(PetscObjectComm((PetscObject)svd),1,"The number of initial vectors is larger than ncv");
    ierr = IPOrthonormalizeBasis_Private(svd->ip,&svd->nini,&svd->IS,svd->V);CHKERRQ(ierr);
  }
  if (svd->ninil<0 && svd->U) { /* skip this if the solver is not using a left basis */
    svd->ninil = -svd->ninil;
    if (svd->ninil>svd->ncv) SETERRQ(PetscObjectComm((PetscObject)svd),1,"The number of left initial vectors is larger than ncv");
    ierr = IPOrthonormalizeBasis_Private(svd->ip,&svd->ninil,&svd->ISL,svd->U);CHKERRQ(ierr);
  }

  ierr = PetscLogEventEnd(SVD_SetUp,svd,0,0,0);CHKERRQ(ierr);
  svd->setupcalled = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SVDSetInitialSpace"
/*@
   SVDSetInitialSpace - Specify a basis of vectors that constitute the initial
   (right) space, that is, a rough approximation to the right singular subspace
   from which the solver starts to iterate.

   Collective on SVD and Vec

   Input Parameter:
+  svd   - the singular value solver context
.  n     - number of vectors
-  is    - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one SVDSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted singular space. Then, convergence may be faster.

   Level: intermediate
@*/
PetscErrorCode SVDSetInitialSpace(SVD svd,PetscInt n,Vec *is)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative"); 
  ierr = SlepcBasisReference_Private(n,is,&svd->nini,&svd->IS);CHKERRQ(ierr);
  if (n>0) svd->setupcalled = 0;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SVDSetInitialSpaceLeft"
/*@
   SVDSetInitialSpaceLeft - Specify a basis of vectors that constitute the initial
   left space, that is, a rough approximation to the left singular subspace
   from which the solver starts to iterate.

   Collective on SVD and Vec

   Input Parameter:
+  svd   - the singular value solver context
.  n     - number of vectors
-  is    - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one SVDSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted singular space. Then, convergence may be faster.

   Level: intermediate
@*/
PetscErrorCode SVDSetInitialSpaceLeft(SVD svd,PetscInt n,Vec *is)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative"); 
  ierr = SlepcBasisReference_Private(n,is,&svd->ninil,&svd->ISL);CHKERRQ(ierr);
  if (n>0) svd->setupcalled = 0;
  PetscFunctionReturn(0);
}

