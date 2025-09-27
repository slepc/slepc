/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves the same problem as in ex5, but with a user-defined sorting criterion.\n\n"
  "It is a standard nonsymmetric eigenproblem with real eigenvalues and the rightmost eigenvalue is known to be 1.\n"
  "This example illustrates how the user can set a custom spectrum selection.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/

PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx);
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);

int main(int argc,char **argv)
{
  Vec            v0;              /* initial vector */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation associated */
  PetscReal      tol=PETSC_SMALL;
  PetscScalar    target=0.5;
  PetscInt       N,m=15,nev;
  char           str[50];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL));
  N = m*(m+1)/2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%" PetscInt_FMT " (m=%" PetscInt_FMT ")\n",N,m));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-target",&target,NULL));
  PetscCall(SlepcSNPrintfScalar(str,sizeof(str),target,PETSC_FALSE));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Searching closest eigenvalues to the right of %s.\n\n",str));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatMarkovModel(m,A));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));

  /*
     Set the custom comparing routine in order to obtain the eigenvalues
     closest to the target on the right only
  */
  PetscCall(EPSSetEigenvalueComparison(eps,MyEigenSort,&target));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /*
     Set the preconditioner based on A - target * I
  */
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetShift(st,target));

  /*
     Set the initial vector. This is optional, if not done the initial
     vector is set to random values
  */
  PetscCall(MatCreateVecs(A,&v0,NULL));
  PetscCall(VecSet(v0,1.0));
  PetscCall(EPSSetInitialSpace(eps,1,&v0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&v0));
  PetscCall(SlepcFinalize());
  return 0;
}

PetscErrorCode MatMarkovModel(PetscInt m,Mat A)
{
  const PetscReal cst = 0.5/(PetscReal)(m-1);
  PetscReal       pd,pu;
  PetscInt        Istart,Iend,i,j,jmax,ix=0;

  PetscFunctionBeginUser;
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=1;i<=m;i++) {
    jmax = m-i+1;
    for (j=1;j<=jmax;j++) {
      ix = ix + 1;
      if (ix-1<Istart || ix>Iend) continue;  /* compute only owned rows */
      if (j!=jmax) {
        pd = cst*(PetscReal)(i+j-1);
        /* north */
        if (i==1) PetscCall(MatSetValue(A,ix-1,ix,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix,pd,INSERT_VALUES));
        /* east */
        if (j==1) PetscCall(MatSetValue(A,ix-1,ix+jmax-1,2*pd,INSERT_VALUES));
        else PetscCall(MatSetValue(A,ix-1,ix+jmax-1,pd,INSERT_VALUES));
      }
      /* south */
      pu = 0.5 - cst*(PetscReal)(i+j-3);
      if (j>1) PetscCall(MatSetValue(A,ix-1,ix-2,pu,INSERT_VALUES));
      /* west */
      if (i>1) PetscCall(MatSetValue(A,ix-1,ix-jmax-2,pu,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, the preferred value is the one closest to the target,
    but on the right side.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscScalar target = *(PetscScalar*)ctx;
  PetscReal   da,db;
  PetscBool   aisright,bisright;

  PetscFunctionBeginUser;
  if (PetscRealPart(target) < PetscRealPart(ar)) aisright = PETSC_TRUE;
  else aisright = PETSC_FALSE;
  if (PetscRealPart(target) < PetscRealPart(br)) bisright = PETSC_TRUE;
  else bisright = PETSC_FALSE;
  if (aisright == bisright) {
    /* both are on the same side of the target */
    da = SlepcAbsEigenvalue(ar-target,ai);
    db = SlepcAbsEigenvalue(br-target,bi);
    if (da < db) *r = -1;
    else if (da > db) *r = 1;
    else *r = 0;
  } else if (aisright && !bisright) *r = -1; /* 'a' is on the right */
  else *r = 1;  /* 'b' is on the right */
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test11_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur arnoldi lapack}} -st_type sinvert
      test:
         suffix: 1_ks_cayley
         args: -st_type cayley -st_cayley_antishift 1

   test:
      suffix: 2
      args: -target 0.77 -eps_type gd -eps_nev 4 -eps_tol 1e-7 -eps_gd_krylov_start -eps_gd_blocksize 3
      requires: double
      filter: sed -e "s/[+-]0\.0*i//g"

TEST*/
