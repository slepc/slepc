/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The pdde_stability problem is a complex-symmetric QEP from the stability
   analysis of a discretized partial delay-differential equation. It requires
   complex scalars.
*/

static char help[] = "Stability analysis of a discretized partial delay-differential equation.\n\n"
  "The command line options are:\n"
  "  -m <m>, grid size, the matrices have dimension n=m*m.\n"
  "  -c <a0,b0,a1,b1,a2,b2,phi1>, comma-separated list of 7 real parameters.\n\n";

#include <slepcpep.h>

#define NMAT 3

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, the preferred value is the one with absolute value closest to 1.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscReal aa,ab;

  PetscFunctionBeginUser;
  aa = PetscAbsReal(SlepcAbsEigenvalue(ar,ai)-PetscRealConstant(1.0));
  ab = PetscAbsReal(SlepcAbsEigenvalue(br,bi)-PetscRealConstant(1.0));
  *r = aa > ab ? 1 : (aa < ab ? -1 : 0);
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       m=15,n,II,Istart,Iend,i,j,k;
  PetscReal      h,xi,xj,c[7] = { 2, .3, -2, .2, -2, -.3, -PETSC_PI/2 };
  PetscScalar    alpha,beta,gamma;
  PetscBool      flg,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
#if !PetscDefined(USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This example requires complex scalars");
#endif

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  n = m*m;
  h = PETSC_PI/(m+1);
  gamma = PetscExpScalar(PETSC_i*c[6]);
  gamma = gamma/PetscAbsScalar(gamma);
  k = 7;
  PetscCall(PetscOptionsGetRealArray(NULL,NULL,"-c",c,&k,&flg));
  PetscCheck(!flg || k==7,PETSC_COMM_WORLD,PETSC_ERR_USER,"The number of parameters -c should be 7, you provided %" PetscInt_FMT,k);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nPDDE stability, n=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n\n",n,m));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Compute the polynomial matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* initialize matrices */
  for (i=0;i<NMAT;i++) {
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A[i]));
    PetscCall(MatSetSizes(A[i],PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(A[i]));
  }
  PetscCall(MatGetOwnershipRange(A[0],&Istart,&Iend));

  /* A[1] has a pattern similar to the 2D Laplacian */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    xi = (i+1)*h; xj = (j+1)*h;
    alpha = c[0]+c[1]*PetscSinReal(xi)+gamma*(c[2]+c[3]*xi*(1.0-PetscExpReal(xi-PETSC_PI)));
    beta = c[0]+c[1]*PetscSinReal(xj)-gamma*(c[2]+c[3]*xj*(1.0-PetscExpReal(xj-PETSC_PI)));
    PetscCall(MatSetValue(A[1],II,II,alpha+beta-4.0/(h*h),INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A[1],II,II-1,1.0/(h*h),INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[1],II,II+1,1.0/(h*h),INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[1],II,II-m,1.0/(h*h),INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[1],II,II+m,1.0/(h*h),INSERT_VALUES));
  }

  /* A[0] and A[2] are diagonal */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    xi = (i+1)*h; xj = (j+1)*h;
    alpha = c[4]+c[5]*xi*(PETSC_PI-xi);
    beta = c[4]+c[5]*xj*(PETSC_PI-xj);
    PetscCall(MatSetValue(A[0],II,II,alpha,INSERT_VALUES));
    PetscCall(MatSetValue(A[2],II,II,beta,INSERT_VALUES));
  }

  /* assemble matrices */
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY));
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PEPSetOperators(pep,NMAT,A));
  PetscCall(PEPSetEigenvalueComparison(pep,MyEigenSort,NULL));
  PetscCall(PEPSetDimensions(pep,4,PETSC_DETERMINE,PETSC_DETERMINE));
  PetscCall(PEPSetFromOptions(pep));
  PetscCall(PEPSolve(pep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PEPDestroy(&pep));
  for (i=0;i<NMAT;i++) PetscCall(MatDestroy(&A[i]));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: complex

   test:
      suffix: 1
      args: -pep_type {{toar qarnoldi linear}} -pep_ncv 25 -terse
      requires: complex double

TEST*/
