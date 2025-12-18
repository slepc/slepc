/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test EPSSetArbitrarySelection.\n\n";

#include <slepceps.h>

PetscErrorCode MyArbitrarySelection(PetscScalar eigr,PetscScalar eigi,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx)
{
  Vec             xref = *(Vec*)ctx;

  PetscFunctionBeginUser;
  PetscCall(VecDot(xr,xref,rr));
  *rr = PetscAbsScalar(*rr);
  if (ri) *ri = 0.0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MyArbitrarySelectionDestroy(void **ctx)
{
  Vec             xref = *(Vec*)*ctx;

  PetscFunctionBeginUser;
  PetscCall(VecDestroy(&xref));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrices */
  EPS            eps;         /* eigenproblem solver context */
  PetscScalar    seigr,seigi;
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON;
  Vec            sxr,sxi;
  PetscInt       n=30,i,Istart,Iend,nconv;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nTridiagonal with zero diagonal, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create matrix tridiag([-1 0 -1])
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve eigenproblem and store some solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>0) {
    PetscCall(MatCreateVecs(A,&sxr,NULL));
    PetscCall(MatCreateVecs(A,&sxi,NULL));
    PetscCall(EPSGetEigenpair(eps,0,&seigr,&seigi,sxr,sxi));
    PetscCall(VecDestroy(&sxi));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve eigenproblem using an arbitrary selection
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(EPSSetArbitrarySelection(eps,MyArbitrarySelection,&sxr));
    PetscCall(EPSSetArbitrarySelectionContextDestroy(eps,MyArbitrarySelectionDestroy));
    PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE));
    PetscCall(EPSSolve(eps));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  } else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Problem: no eigenpairs converged.\n"));

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_max_it 5000 -st_pc_type jacobi
      output_file: output/test13_1.out
      filter: sed -e "s/-1.98975/-1.98974/"
      test:
         suffix: 1
         args: -eps_type {{krylovschur gd}}
      test:
         suffix: 1_jd
         args: -eps_type jd -st_ksp_type gmres
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 2
         args: -eps_non_hermitian -eps_type {{krylovschur gd jd}}
      test:
         suffix: 2_gd2
         args: -eps_non_hermitian -eps_type gd -eps_gd_double_expansion
         timeoutfactor: 2

TEST*/
