/*
     This file contains some simple default routines for common operations.  
*/
#include "src/eps/epsimpl.h"   /*I "slepceps.h" I*/
#include "slepcblaslapack.h"

#undef __FUNCT__  
#define __FUNCT__ "EPSDestroy_Default"
PetscErrorCode EPSDestroy_Default(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (eps->data) {ierr = PetscFree(eps->data);CHKERRQ(ierr);}

  /* free work vectors */
  ierr = EPSDefaultFreeWork(eps);CHKERRQ(ierr);
  ierr = EPSFreeSolution(eps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "EPSBackTransform_Default"
PetscErrorCode EPSBackTransform_Default(EPS eps)
{
  PetscErrorCode ierr;
  int            i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  for (i=0;i<eps->nconv;i++) {
    ierr = STBackTransform(eps->OP,&eps->eigr[i],&eps->eigi[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "EPSComputeVectors_Default"
/*
  EPSComputeVectors_Default - Compute eigenvectors from the vectors
  provided by the eigensolver. This version just copies the vectors
  and is intended for solvers such as power that provide the eigenvector.
 */
PetscErrorCode EPSComputeVectors_Default(EPS eps)
{
  PetscErrorCode ierr;
  int            i;

  PetscFunctionBegin;
  for (i=0;i<eps->nconv;i++) {
    ierr = VecCopy(eps->V[i],eps->AV[i]);CHKERRQ(ierr);
    if (eps->solverclass == EPS_TWO_SIDE) {
      ierr = VecCopy(eps->W[i],eps->AW[i]);CHKERRQ(ierr);
    }
  }
  eps->evecsavailable = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "EPSComputeVectors_Schur"
/*
  EPSComputeVectors_Schur - Compute eigenvectors from the vectors
  provided by the eigensolver. This version is intended for solvers 
  that provide Schur vectors. Given the partial Schur decomposition
  OP*V=V*T, the following steps are performed:
      1) compute eigenvectors of T: T*Z=Z*D
      2) compute eigenvectors of OP: X=V*Z
  If left eigenvectors are required then also do Z'*Tl=D*Z', Y=W*Z
 */
PetscErrorCode EPSComputeVectors_Schur(EPS eps)
{
  PetscErrorCode ierr;
  int            i,mout,info,ncv=eps->ncv;
  PetscScalar    *Z,*work;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif
  
  PetscFunctionBegin;
  if (eps->ishermitian) {
    ierr = EPSComputeVectors_Default(eps);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

#if defined(PETSC_BLASLAPACK_ESSL_ONLY)
  SETERRQ(PETSC_ERR_SUP,"TREVC - Lapack routine is unavailable.");
#endif 

  ierr = PetscMalloc(ncv*ncv*sizeof(PetscScalar),&Z);CHKERRQ(ierr);
  ierr = PetscMalloc(3*ncv*sizeof(PetscScalar),&work);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscMalloc(ncv*sizeof(PetscReal),&rwork);CHKERRQ(ierr);
#endif

  /* right eigenvectors */
#if !defined(PETSC_USE_COMPLEX)
  LAtrevc_("R","A",PETSC_NULL,&ncv,eps->T,&ncv,PETSC_NULL,&ncv,Z,&ncv,&ncv,&mout,work,&info,1,1);
#else
  LAtrevc_("R","A",PETSC_NULL,&ncv,eps->T,&ncv,PETSC_NULL,&ncv,Z,&ncv,&ncv,&mout,work,rwork,&info,1,1);
#endif
  if (info) SETERRQ1(PETSC_ERR_LIB,"Error in Lapack xTREVC %i",info);

  for (i=0;i<eps->nconv;i++) {
    ierr = VecCopy(eps->V[i],eps->AV[i]);CHKERRQ(ierr);
  }
  ierr = EPSReverseProjection(eps,eps->AV,Z,0,ncv,eps->work);CHKERRQ(ierr);
   
  /* left eigenvectors */
  if (eps->solverclass == EPS_TWO_SIDE) {
#if !defined(PETSC_USE_COMPLEX)
    LAtrevc_("L","A",PETSC_NULL,&ncv,eps->Tl,&ncv,PETSC_NULL,&ncv,Z,&ncv,&ncv,&mout,work,&info,1,1);
#else
    LAtrevc_("L","A",PETSC_NULL,&ncv,eps->Tl,&ncv,PETSC_NULL,&ncv,Z,&ncv,&ncv,&mout,work,rwork,&info,1,1);
#endif
    if (info) SETERRQ1(PETSC_ERR_LIB,"Error in Lapack xTREVC %i",info);

    for (i=0;i<eps->nconv;i++) {
      ierr = VecCopy(eps->W[i],eps->AW[i]);CHKERRQ(ierr);
    }
    ierr = EPSReverseProjection(eps,eps->AW,Z,0,ncv,eps->work);CHKERRQ(ierr);
  }
   
  ierr = PetscFree(Z);CHKERRQ(ierr);
  ierr = PetscFree(work);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree(rwork);CHKERRQ(ierr);
#endif
  eps->evecsavailable = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "EPSDefaultGetWork"
/*
  EPSDefaultGetWork - Gets a number of work vectors.

  Input Parameters:
+ eps  - eigensolver context
- nw   - number of work vectors to allocate

  Notes:
  Call this only if no work vectors have been allocated.

 */
PetscErrorCode EPSDefaultGetWork(EPS eps, int nw)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (eps->nwork != nw) {
    if (eps->nwork > 0) {
      ierr = VecDestroyVecs(eps->work,eps->nwork); CHKERRQ(ierr);
    }
    eps->nwork = nw;
    ierr = VecDuplicateVecs(eps->vec_initial,nw,&eps->work); CHKERRQ(ierr);
    PetscLogObjectParents(eps,nw,eps->work);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "EPSDefaultFreeWork"
/*
  EPSDefaultFreeWork - Free work vectors.

  Input Parameters:
. eps  - eigensolver context

 */
PetscErrorCode EPSDefaultFreeWork(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (eps->work)  {
    ierr = VecDestroyVecs(eps->work,eps->nwork); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
