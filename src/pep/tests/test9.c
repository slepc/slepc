/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates use of PEPSetEigenvalueComparison().\n\n"
  "Based on butterfly.c. The command line options are:\n"
  "  -m <m>, grid size, the dimension of the matrices is n=m*m.\n"
  "  -c <array>, problem parameters, must be 10 comma-separated real values.\n\n";

#include <slepcpep.h>

#define NMAT 5

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, eigenvalues are sorted by magnitude but those with
    positive real part are preferred.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscReal rea,reb;

  PetscFunctionBeginUser;
#if PetscDefined(USE_COMPLEX)
  rea = PetscRealPart(ar); reb = PetscRealPart(br);
#else
  rea = ar; reb = br;
#endif
  *r = rea<0.0? 1: (reb<0.0? -1: PetscSign(SlepcAbsEigenvalue(br,bi)-SlepcAbsEigenvalue(ar,ai)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n,m=8,k,II,Istart,Iend,i,j;
  PetscReal      c[10] = { 0.6, 1.3, 1.3, 0.1, 0.1, 1.2, 1.0, 1.0, 1.2, 1.0 };
  PetscBool      flg;
  PetscBool      terse;
  const char     *prefix;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  n = m*m;
  k = 10;
  PetscCall(PetscOptionsGetRealArray(NULL,NULL,"-c",c,&k,&flg));
  PetscCheck(!flg || k==10,PETSC_COMM_WORLD,PETSC_ERR_USER,"The number of parameters -c should be 10, you provided %" PetscInt_FMT,k);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nButterfly problem, n=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n\n",n,m));

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

  /* A0 */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    PetscCall(MatSetValue(A[0],II,II,4.0*c[0]/6.0+4.0*c[1]/6.0,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A[0],II,II-1,c[0]/6.0,INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[0],II,II+1,c[0]/6.0,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[0],II,II-m,c[1]/6.0,INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[0],II,II+m,c[1]/6.0,INSERT_VALUES));
  }

  /* A1 */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    if (j>0) PetscCall(MatSetValue(A[1],II,II-1,c[2],INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[1],II,II+1,-c[2],INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[1],II,II-m,c[3],INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[1],II,II+m,-c[3],INSERT_VALUES));
  }

  /* A2 */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    PetscCall(MatSetValue(A[2],II,II,-2.0*c[4]-2.0*c[5],INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A[2],II,II-1,c[4],INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[2],II,II+1,c[4],INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[2],II,II-m,c[5],INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[2],II,II+m,c[5],INSERT_VALUES));
  }

  /* A3 */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    if (j>0) PetscCall(MatSetValue(A[3],II,II-1,c[6],INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[3],II,II+1,-c[6],INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[3],II,II-m,c[7],INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[3],II,II+m,-c[7],INSERT_VALUES));
  }

  /* A4 */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    PetscCall(MatSetValue(A[4],II,II,2.0*c[8]+2.0*c[9],INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(A[4],II,II-1,-c[8],INSERT_VALUES));
    if (j<m-1) PetscCall(MatSetValue(A[4],II,II+1,-c[8],INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(A[4],II,II-m,-c[9],INSERT_VALUES));
    if (i<m-1) PetscCall(MatSetValue(A[4],II,II+m,-c[9],INSERT_VALUES));
  }

  /* assemble matrices */
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY));
  for (i=0;i<NMAT;i++) PetscCall(MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PEPSetOptionsPrefix(pep,"check_"));
  PetscCall(PEPAppendOptionsPrefix(pep,"myprefix_"));
  PetscCall(PEPGetOptionsPrefix(pep,&prefix));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"PEP prefix is currently: %s\n\n",prefix));

  PetscCall(PEPSetOperators(pep,NMAT,A));
  PetscCall(PEPSetEigenvalueComparison(pep,MyEigenSort,NULL));
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

   test:
      args: -check_myprefix_pep_nev 4 -terse
      requires: double

TEST*/
