/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Simple 1-D nonlinear eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions.\n"
  "  -draw_sol, to draw the computed solution.\n\n";

/*
   Solve 1-D PDE
            -u'' = lambda*u
   on [0,1] subject to
            u(0)=0, u'(1)=u(1)*lambda*kappa/(kappa-lambda)
*/

#include <slepcnep.h>

/*
   User-defined routines
*/
PetscErrorCode FormInitialGuess(Vec);
PetscErrorCode FormFunction(NEP,PetscScalar,Mat,Mat,void*);
PetscErrorCode FormJacobian(NEP,PetscScalar,Mat,void*);
PetscErrorCode CheckSolution(PetscScalar,Vec,PetscReal*,void*);
PetscErrorCode FixSign(Vec);

/*
   User-defined application context
*/
typedef struct {
  PetscScalar kappa;   /* ratio between stiffness of spring and attached mass */
  PetscReal   h;       /* mesh spacing */
} ApplicationCtx;

int main(int argc,char **argv)
{
  NEP            nep;             /* nonlinear eigensolver context */
  Vec            x;               /* eigenvector */
  PetscScalar    lambda;          /* eigenvalue */
  Mat            F,J;             /* Function and Jacobian matrices */
  ApplicationCtx ctx;             /* user-defined context */
  NEPType        type;
  PetscInt       n=128,nev,i,its,maxit,nconv;
  PetscReal      re,im,tol,norm,error;
  PetscBool      draw_sol=PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Nonlinear Eigenproblem, n=%" PetscInt_FMT "\n\n",n));
  ctx.h = 1.0/(PetscReal)n;
  ctx.kappa = 1.0;
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-draw_sol",&draw_sol,NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear eigensolver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPCreate(PETSC_COMM_WORLD,&nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; set Function evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&F));
  PetscCall(MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(F));
  PetscCall(MatSeqAIJSetPreallocation(F,3,NULL));
  PetscCall(MatMPIAIJSetPreallocation(F,3,NULL,1,NULL));

  /*
     Set Function matrix data structure and default Function evaluation
     routine
  */
  PetscCall(NEPSetFunction(nep,F,F,FormFunction,&ctx));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD,&J));
  PetscCall(MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(J));
  PetscCall(MatSeqAIJSetPreallocation(J,3,NULL));
  PetscCall(MatMPIAIJSetPreallocation(J,3,NULL,1,NULL));

  /*
     Set Jacobian matrix data structure and default Jacobian evaluation
     routine
  */
  PetscCall(NEPSetJacobian(nep,J,FormJacobian,&ctx));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSetTolerances(nep,1e-9,PETSC_CURRENT));
  PetscCall(NEPSetDimensions(nep,1,PETSC_DETERMINE,PETSC_DETERMINE));

  /*
     Set solver parameters at runtime
  */
  PetscCall(NEPSetFromOptions(nep));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Initialize application
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Evaluate initial guess
  */
  PetscCall(MatCreateVecs(F,&x,NULL));
  PetscCall(FormInitialGuess(x));
  PetscCall(NEPSetInitialSpace(nep,1,&x));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(NEPSolve(nep));
  PetscCall(NEPGetIterationNumber(nep,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of NEP iterations = %" PetscInt_FMT "\n\n",its));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(NEPGetType(nep,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type));
  PetscCall(NEPGetDimensions(nep,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(NEPGetTolerances(nep,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get number of converged approximate eigenpairs
  */
  PetscCall(NEPGetConverged(nep,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %" PetscInt_FMT "\n\n",nconv));

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k              ||T(k)x||           error\n"
         "   ----------------- ------------------ ------------------\n"));
    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs (in this example they are always real)
      */
      PetscCall(NEPGetEigenpair(nep,i,&lambda,NULL,x,NULL));
      PetscCall(FixSign(x));
      /*
         Compute residual norm and error
      */
      PetscCall(NEPComputeError(nep,i,NEP_ERROR_RELATIVE,&norm));
      PetscCall(CheckSolution(lambda,x,&error,&ctx));

#if PetscDefined(USE_COMPLEX)
      re = PetscRealPart(lambda);
      im = PetscImaginaryPart(lambda);
#else
      re = lambda;
      im = 0.0;
#endif
      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g     %12g\n",(double)re,(double)im,(double)norm,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f         %12g     %12g\n",(double)re,(double)norm,(double)error));
      if (draw_sol) {
        PetscCall(PetscViewerDrawSetPause(PETSC_VIEWER_DRAW_WORLD,-1));
        PetscCall(VecView(x,PETSC_VIEWER_DRAW_WORLD));
      }
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

  PetscCall(NEPDestroy(&nep));
  PetscCall(MatDestroy(&F));
  PetscCall(MatDestroy(&J));
  PetscCall(VecDestroy(&x));
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormInitialGuess - Computes initial guess.

   Input/Output Parameter:
.  x - the solution vector
*/
PetscErrorCode FormInitialGuess(Vec x)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(x,1.0));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
   FormFunction - Computes Function matrix  T(lambda)

   Input Parameters:
.  nep    - the NEP context
.  lambda - the scalar argument
.  ctx    - optional user-defined context, as set by NEPSetFunction()

   Output Parameters:
.  fun - Function matrix
.  B   - optionally different preconditioning matrix
*/
PetscErrorCode FormFunction(NEP nep,PetscScalar lambda,Mat fun,Mat B,void *ctx)
{
  ApplicationCtx *user = (ApplicationCtx*)ctx;
  PetscScalar    A[3],c,d;
  PetscReal      h;
  PetscInt       i,n,j[3],Istart,Iend;
  PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;

  PetscFunctionBeginUser;
  /*
     Compute Function entries and insert into matrix
  */
  PetscCall(MatGetSize(fun,&n,NULL));
  PetscCall(MatGetOwnershipRange(fun,&Istart,&Iend));
  if (Istart==0) FirstBlock=PETSC_TRUE;
  if (Iend==n) LastBlock=PETSC_TRUE;
  h = user->h;
  c = user->kappa/(lambda-user->kappa);
  d = n;

  /*
     Interior grid points
  */
  for (i=(FirstBlock? Istart+1: Istart);i<(LastBlock? Iend-1: Iend);i++) {
    j[0] = i-1; j[1] = i; j[2] = i+1;
    A[0] = A[2] = -d-lambda*h/6.0; A[1] = 2.0*(d-lambda*h/3.0);
    PetscCall(MatSetValues(fun,1,&i,3,j,A,INSERT_VALUES));
  }

  /*
     Boundary points
  */
  if (FirstBlock) {
    i = 0;
    j[0] = 0; j[1] = 1;
    A[0] = 2.0*(d-lambda*h/3.0); A[1] = -d-lambda*h/6.0;
    PetscCall(MatSetValues(fun,1,&i,2,j,A,INSERT_VALUES));
  }

  if (LastBlock) {
    i = n-1;
    j[0] = n-2; j[1] = n-1;
    A[0] = -d-lambda*h/6.0; A[1] = d-lambda*h/3.0+lambda*c;
    PetscCall(MatSetValues(fun,1,&i,2,j,A,INSERT_VALUES));
  }

  /*
     Assemble matrix
  */
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  if (fun != B) {
    PetscCall(MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian - Computes Jacobian matrix  T'(lambda)

   Input Parameters:
.  nep    - the NEP context
.  lambda - the scalar argument
.  ctx    - optional user-defined context, as set by NEPSetJacobian()

   Output Parameters:
.  jac - Jacobian matrix
*/
PetscErrorCode FormJacobian(NEP nep,PetscScalar lambda,Mat jac,void *ctx)
{
  ApplicationCtx *user = (ApplicationCtx*)ctx;
  PetscScalar    A[3],c;
  PetscReal      h;
  PetscInt       i,n,j[3],Istart,Iend;
  PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;

  PetscFunctionBeginUser;
  /*
     Compute Jacobian entries and insert into matrix
  */
  PetscCall(MatGetSize(jac,&n,NULL));
  PetscCall(MatGetOwnershipRange(jac,&Istart,&Iend));
  if (Istart==0) FirstBlock=PETSC_TRUE;
  if (Iend==n) LastBlock=PETSC_TRUE;
  h = user->h;
  c = user->kappa/(lambda-user->kappa);

  /*
     Interior grid points
  */
  for (i=(FirstBlock? Istart+1: Istart);i<(LastBlock? Iend-1: Iend);i++) {
    j[0] = i-1; j[1] = i; j[2] = i+1;
    A[0] = A[2] = -h/6.0; A[1] = -2.0*h/3.0;
    PetscCall(MatSetValues(jac,1,&i,3,j,A,INSERT_VALUES));
  }

  /*
     Boundary points
  */
  if (FirstBlock) {
    i = 0;
    j[0] = 0; j[1] = 1;
    A[0] = -2.0*h/3.0; A[1] = -h/6.0;
    PetscCall(MatSetValues(jac,1,&i,2,j,A,INSERT_VALUES));
  }

  if (LastBlock) {
    i = n-1;
    j[0] = n-2; j[1] = n-1;
    A[0] = -h/6.0; A[1] = -h/3.0-c*c;
    PetscCall(MatSetValues(jac,1,&i,2,j,A,INSERT_VALUES));
  }

  /*
     Assemble matrix
  */
  PetscCall(MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
   CheckSolution - Given a computed solution (lambda,x) check if it
   satisfies the analytic solution.

   Input Parameters:
+  lambda - the computed eigenvalue
-  y      - the computed eigenvector

   Output Parameter:
.  error - norm of difference between the computed and exact eigenvector
*/
PetscErrorCode CheckSolution(PetscScalar lambda,Vec y,PetscReal *error,void *ctx)
{
  PetscScalar    *uu;
  PetscInt       i,n,Istart,Iend;
  PetscReal      nu,x;
  Vec            u;
  ApplicationCtx *user = (ApplicationCtx*)ctx;

  PetscFunctionBeginUser;
  nu = PetscSqrtReal(PetscRealPart(lambda));
  PetscCall(VecDuplicate(y,&u));
  PetscCall(VecGetSize(u,&n));
  PetscCall(VecGetOwnershipRange(y,&Istart,&Iend));
  PetscCall(VecGetArray(u,&uu));
  for (i=Istart;i<Iend;i++) {
    x = (i+1)*user->h;
    uu[i-Istart] = PetscSinReal(nu*x);
  }
  PetscCall(VecRestoreArray(u,&uu));
  PetscCall(VecNormalize(u,NULL));
  PetscCall(VecAXPY(u,-1.0,y));
  PetscCall(VecNorm(u,NORM_2,error));
  PetscCall(VecDestroy(&u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
   FixSign - Force the eigenfunction to be real and positive, since
   some eigensolvers may return the eigenvector multiplied by a
   complex number of modulus one.

   Input/Output Parameter:
.  x - the computed vector
*/
PetscErrorCode FixSign(Vec x)
{
  PetscMPIInt       rank;
  PetscScalar       sign=0.0;
  const PetscScalar *xx;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  if (!rank) {
    PetscCall(VecGetArrayRead(x,&xx));
    sign = *xx/PetscAbsScalar(*xx);
    PetscCall(VecRestoreArrayRead(x,&xx));
  }
  PetscCallMPI(MPI_Bcast(&sign,1,MPIU_SCALAR,0,PETSC_COMM_WORLD));
  PetscCall(VecScale(x,1.0/sign));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   build:
      requires: !single

   test:
      suffix: 1
      args: -nep_target 4
      filter: sed -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
      requires: !single

   testset:
      args: -nep_type nleigs -nep_target 10 -nep_nev 4
      test:
         suffix: 2
         filter: sed -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/0\.\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
         args: -rg_interval_endpoints 4,200
         requires: !single !complex
      test:
         suffix: 2_complex
         filter: sed -e "s/[+-]0.[0-9]*i//" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/0\.\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
         output_file: output/ex20_2.out
         args: -rg_interval_endpoints 4,200,-.1,.1 -nep_target_real
         requires: !single complex

   testset:
      args: -nep_type nleigs -nep_target 10 -nep_nev 4 -nep_ncv 18 -nep_two_sided {{0 1}} -nep_nleigs_full_basis
      test:
         suffix: 3
         filter: sed -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/0\.\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
         output_file: output/ex20_2.out
         args: -rg_interval_endpoints 4,200
         requires: !single !complex
      test:
         suffix: 3_complex
         filter: sed -e "s/[+-]0.[0-9]*i//" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/0\.\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
         output_file: output/ex20_2.out
         args: -rg_interval_endpoints 4,200,-.1,.1
         requires: !single complex

   test:
      suffix: 4
      args: -n 256 -nep_nev 2 -nep_target 10
      filter: sed -e "s/[+-]0\.0*i//g" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/ Number of NEP iterations = \([0-9]*\)/ Number of NEP iterations = /"
      requires: !single

TEST*/
