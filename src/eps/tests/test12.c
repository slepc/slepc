/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Diagonal eigenproblem. Illustrates use of shell preconditioner.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n"
  "  -seed <s>, where <s> = seed for random number generation.\n\n";

#include <slepceps.h>

typedef struct {
  PetscScalar target;
} PCCtx;

PetscErrorCode PCApply_User(PC pc,Vec x,Vec y)
{
  PetscInt    i,rstart,rend;
  PetscScalar *yarray;
  PCCtx       *ctx;

  PetscFunctionBeginUser;
  PetscCall(PCShellGetContext(pc,&ctx));
  PetscCall(VecCopy(x,y));
  PetscCall(VecGetOwnershipRange(y,&rstart,&rend));
  PetscCall(VecGetArray(y,&yarray));
  for (i=0;i<rend-rstart;i++) yarray[i] /= (i-rstart+1-ctx->target);
  PetscCall(VecRestoreArray(y,&yarray));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  Vec            v0;          /* initial vector */
  PetscRandom    rand;
  PetscReal      tol=PETSC_SMALL;
  PetscInt       n=30,i,Istart,Iend,seed=0x12345678;
  ST             st;
  KSP            ksp;
  PC             pc;
  PCCtx          ctx;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Eigenproblem, n=%" PetscInt_FMT "\n\n",n));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(A,i,i,i+1,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));
  PetscCall(EPSSetTolerances(eps,tol,PETSC_CURRENT));
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCSHELL));
  PetscCall(PCShellSetApply(pc,PCApply_User));
  PetscCall(EPSGetTarget(eps,&ctx.target));
  PetscCall(PCShellSetContext(pc,&ctx));

  /* set random initial vector */
  PetscCall(MatCreateVecs(A,&v0,NULL));
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&rand));
  PetscCall(PetscRandomSetFromOptions(rand));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-seed",&seed,NULL));
  PetscCall(PetscRandomSetSeed(rand,seed));
  PetscCall(PetscRandomSeed(rand));
  PetscCall(VecSetRandom(v0,rand));
  PetscCall(EPSSetInitialSpace(eps,1,&v0));
  /* call the solver */
  PetscCall(EPSSolve(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&v0));
  PetscCall(PetscRandomDestroy(&rand));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -eps_nev 4 -eps_target 31
      requires: !single
      output_file: output/test12_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur subspace arnoldi power}} -st_type sinvert
      test:
         suffix: 1_gd
         args: -eps_type {{gd jd}}
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion

TEST*/
