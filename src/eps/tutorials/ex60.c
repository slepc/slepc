/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "2-D discrete Gross-Pitaevskii model solved as a Nonlinear Eigenvalue Problem with eigenvector dependence (NEPv).\n\n"
  "Implemented by embedding an EPS inside a nonlinear SNES loop.\n\n"
  "The Gross-Pitaevskii (GP) equation models the ground state of a rotating Bose-Einstein Condensate (BEC).\n"
  "The goal is to self-consistently find the macroscopic wavefunction (eigenvector) and its corresponding\n"
  "chemical potential or ground-state energy (eigenvalue) under a magnetic trap.\n\n"
  "The nonlinear Hamiltonian matrix is defined as H(x) = A0 + beta * Diag(rho(x)), where:\n"
  "   A0   = base linear Hamiltonian containing 2-D kinetic energy, a harmonic trap potential, and angular momentum rotation.\n"
  "   rho  = condensate density vector computed from the ground-state wavefunction x as rho = |x|^2.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid points per dimension.\n"
  "  -beta <beta>, where <beta> = real scaling parameter controlling the nonlinearity and strength of the atomic interactions.\n"
  "  -symm <bool>: trap potential profile (true = symmetric, false = asymmetric).\n"
  "  -verbose, to print the converged Hamiltonian diagonal and its sum.\n\n";

#include <slepceps.h>

/*
  Context structure to store the necessary data regarding the discrete Gross-Pitaevskii (DGP) problem
  as well as the objects needed for the self-consistent field (SCF) iteration.
*/
typedef struct {
  // DGP context

  // Inputs
  PetscInt  N;        // Number of points per dimension
  PetscReal beta;     // Parameter controlling the nonlinearity
  PetscBool symm;     // True: Symmetric potential, False: Asymmetric potential

  // Constants
  PetscReal L;        // Length of the semi-domain [-L, L]
  PetscReal h;        // Step size
  PetscReal omega;    // Angular frequency

  // Working vectors and matrices
  Vec       rho;     // Vector to store the condensate density
  Vec       z;       // Intermediate vector to store the result z = beta * rho
  Mat       A0;      // Base linear Hamiltonian

  // SCF context
  EPS       eps;     // Solver for each H
  Mat       H;       // Hamiltonian
  Vec       x;       // Eigenvector
} DGPContext;

/*
  Initialize the context for the 2D discrete Gross-Pitaevskii problem.

  This function allocates the necessary memory and builds the base linear
  Hamiltonian A0 as

    A0 = h^2 * ( -0.5*L + V - omega*1i*L_z )

  where
    L is the discrete 2-D Laplacian,
    L_z is the discrete 2-D angular momentum operator (without scaling), and
    V is a diagonal matrix with the potential for each coordinate (x,y):
      (x.^2+y.^2)/2 for the symmetric case, (x.^2+100*y.^2)/2 for the non-symmetric case.

  Arguments:
    N       - Grid size per dimension
    beta    - Parameter controlling the repulsive interaction of the condensate
    symm    - Defines the shape of the magnetic trap (True = symmetric, False = asymmetric)
    ctx_out - Pointer to the memory address where the created context will be stored
*/
PetscErrorCode DGPCreate(MPI_Comm comm,PetscInt N,PetscReal beta,PetscBool symm,DGPContext **ctx_out)
{
  DGPContext  *ctx;

  PetscInt    Istart,Iend,i,II,j;
  PetscReal   h2;
  PetscReal   x,y,V;
  PetscScalar alpha;

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&ctx));

  // Inputs
  ctx->N=N;
  ctx->beta=beta;
  ctx->symm=symm;

  // Constants
  ctx->L=1.0;
  ctx->h=2.0*ctx->L/(ctx->N+1.0);
  ctx->omega=0.85;
  h2=ctx->h*ctx->h;

  // 1. Create A0
  PetscCall(MatCreate(comm,&ctx->A0));
  PetscCall(MatSetSizes(ctx->A0,PETSC_DECIDE,PETSC_DECIDE,N*N,N*N));
  PetscCall(MatSetFromOptions(ctx->A0));

  // 2. Fill A0
  PetscCall(MatGetOwnershipRange(ctx->A0,&Istart,&Iend));

  for (II=Istart;II<Iend;II++) {
    i=II/ctx->N; j=II-i*ctx->N;

    // Potential for mesh point with coordinates (x,y)
    x=-ctx->L+(j+1)*ctx->h;
    y=-ctx->L+(i+1)*ctx->h;
    if (ctx->symm) V=0.5*(x*x+y*y);
    else V=0.5*(x*x+100.0*y*y);

    alpha = -h2*ctx->omega/(2.0*ctx->h)*PETSC_i;
    PetscCall(MatSetValue(ctx->A0,II,II,0.5*4.0+h2*V,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(ctx->A0,II,II-ctx->N,-0.5+alpha*x,INSERT_VALUES));
    if (i<ctx->N-1) PetscCall(MatSetValue(ctx->A0,II,II+ctx->N,-0.5-alpha*x,INSERT_VALUES));
    if (j>0) PetscCall(MatSetValue(ctx->A0,II,II-1,-0.5-alpha*y,INSERT_VALUES));
    if (j<ctx->N-1) PetscCall(MatSetValue(ctx->A0,II,II+1,-0.5+alpha*y,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(ctx->A0,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(ctx->A0,MAT_FINAL_ASSEMBLY));

  // 3. Create vectors to store the condensate density
  PetscCall(MatCreateVecs(ctx->A0,&ctx->rho,&ctx->z));

  *ctx_out=ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Configure the SLEPc objects necessary for the SCF iteration.

  This function prepares the Hamiltonian matrix by cloning the structure of A0,
  initializes the vector to store the condensate state, and configures the eigenvalue
  solver (EPS) by defining the problem type and adjusting its tolerance.

  Arguments:
    ctx - Pointer to the previously initialized DGP context
    tol - Desired tolerance for the internal eigenvalue solver
*/
PetscErrorCode DGPSetupSCF(DGPContext *ctx,PetscReal tol)
{
  MPI_Comm comm;

  PetscFunctionBeginUser;
  comm = PetscObjectComm((PetscObject)ctx->A0);
  /* Prepare the H matrix */
  PetscCall(MatDuplicate(ctx->A0,MAT_DO_NOT_COPY_VALUES,&ctx->H));

  /* Configure the eigenvector x */
  PetscCall(MatCreateVecs(ctx->A0,&ctx->x,NULL));

  /* Configure SLEPc EPS */
  PetscCall(EPSCreate(comm,&ctx->eps));
  PetscCall(EPSSetProblemType(ctx->eps,EPS_HEP));
  PetscCall(EPSSetWhichEigenpairs(ctx->eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSetDimensions(ctx->eps,1,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(EPSSetTolerances(ctx->eps,tol,PETSC_DECIDE));
  PetscCall(EPSSetFromOptions(ctx->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Generate the initial guess X0 by solving the base linear Hamiltonian.

  This function solves the eigenvalue problem (A0 * x = E * x) and fills
  the vector x with the lowest energy eigenvector obtained. This
  initial guess improves the convergence of the SCF loop.

  Arguments:
    ctx - Pointer to the previously initialized DGP context
*/
PetscErrorCode DGPGenerateInitialGuess(DGPContext *ctx)
{
  PetscInt nconv;
  MPI_Comm comm;

  PetscFunctionBeginUser;
  comm = PetscObjectComm((PetscObject)ctx->A0);
  /* 1. Configure and solve the base linear problem A0 */
  PetscCall(EPSSetOperators(ctx->eps,ctx->A0,NULL));
  PetscCall(EPSSolve(ctx->eps));

  /* Safety check: guarantee that SLEPc found the solution */
  PetscCall(EPSGetConverged(ctx->eps,&nconv));
  PetscCheck(nconv>=1,comm,PETSC_ERR_NOT_CONVERGED,"SLEPc could not find the initial linear ground state.");

  /* 2. Extract the eigenvector and store it in x */
  PetscCall(EPSGetEigenvector(ctx->eps,0,ctx->x,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Compute the condensate density from the ground-state eigenvector.

  This function computes the density rho(x) as the squared magnitude
  of the components of the eigenvector X0. The formula used is: rho(x) = |x|.^2.

  Arguments:
    ctx     - Pointer to the initialized DGP context
    rho_out - Pre-created vector where the computed density will be stored
*/
PetscErrorCode DGPCalculateDensity(DGPContext *ctx,Vec rho_out)
{
  PetscInt          i,n_loc;
  const PetscScalar *x_local;
  PetscScalar       *rho_local;

  PetscFunctionBeginUser;
  PetscCall(VecSet(rho_out,0.0));
  PetscCall(VecGetLocalSize(rho_out,&n_loc));
  PetscCall(VecGetArray(rho_out,&rho_local));

  PetscCall(VecGetArrayRead(ctx->x,&x_local));

  for (i=0;i<n_loc;i++) {
    /* Compute the squared magnitude: rho = Real^2 + Imag^2 */
    PetscReal mod=PetscAbsScalar(x_local[i]);
    rho_local[i]=mod*mod;
  }

  PetscCall(VecRestoreArrayRead(ctx->x,&x_local));
  PetscCall(VecRestoreArray(rho_out,&rho_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Build the Hamiltonian from an input density.

  This function assembles the updated Hamiltonian matrix stored in the context
  using the formula: H = A0 + beta * Diag(rho_in).

  Arguments:
    ctx    - Pointer to the initialized DGP context (ctx->H is updated)
    rho_in - Vector with the input density
*/
PetscErrorCode DGPBuildHamiltonian(DGPContext *ctx,Vec rho_in)
{
  PetscFunctionBeginUser;
  // 1. Build H = A0 + beta * Diag(rho_in)
  PetscCall(MatCopy(ctx->A0,ctx->H,SAME_NONZERO_PATTERN));
  PetscCall(VecCopy(rho_in,ctx->z));

  PetscCall(VecScale(ctx->z,ctx->beta));
  PetscCall(MatDiagonalSet(ctx->H,ctx->z,ADD_VALUES));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Evaluation function (Callback) for the SNES nonlinear solver in DGP.

  In each SNES iteration, this function receives a proposed density (rho_in),
  builds the Hamiltonian, solves the eigenvalue equation, computes the
  resulting density, and returns the residual F = rho_out - rho_in.

  Arguments:
    snes     - The nonlinear solver context
    rho_in   - Condensate density proposed by SNES
    F        - Vector where the computed residual will be stored
    ctx_void - Pointer to the user's DGP context (DGPContext)
*/
PetscErrorCode DGPIterationSCF(SNES snes,Vec rho_in,Vec F,void *ctx_void)
{
  DGPContext *ctx=(DGPContext*)ctx_void;
  PetscInt   nconv;
  MPI_Comm   comm;

  PetscFunctionBeginUser;
  comm = PetscObjectComm((PetscObject)ctx->A0);
  /* 1. Build the physics (Hamiltonian) with the guess proposed by SNES */
  PetscCall(DGPBuildHamiltonian(ctx,rho_in));

  /* 2. Solve the eigenvalue problem with the current H */
  PetscCall(EPSSetOperators(ctx->eps,ctx->H,NULL));
  PetscCall(EPSSolve(ctx->eps));

  // Safety check: verify that SLEPc found the ground state
  PetscCall(EPSGetConverged(ctx->eps,&nconv));
  if (nconv<1) PetscCall(PetscPrintf(comm,"Warning: SLEPc did not find the ground state in this iteration.\n"));

  /* 3. Extract the new eigenvector and store it directly in ctx->x */
  PetscCall(EPSGetEigenvector(ctx->eps,0,ctx->x,NULL));

  /* 4. Compute the new density (rho_out) generated by this eigenvector */
  // Store the result in ctx->rho
  PetscCall(DGPCalculateDensity(ctx,ctx->rho));

  /* 5. Calculate the nonlinear residual: F = rho_out - rho_in */
  PetscCall(VecWAXPY(F,-1.0,ctx->rho,rho_in));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Free the memory associated with the DGP context.

  This function destroys all internal PETSc objects created during
  initialization and frees the memory of the main structure.

  Arguments:
    ctx - Pointer to the DGP context (set to NULL upon completion)
*/
PetscErrorCode DGPDestroy(DGPContext **ctx)
{
  PetscFunctionBeginUser;
  if (!*ctx) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(MatDestroy(&(*ctx)->A0));
  PetscCall(VecDestroy(&(*ctx)->rho));
  PetscCall(VecDestroy(&(*ctx)->z));

  PetscCall(MatDestroy(&(*ctx)->H));
  PetscCall(VecDestroy(&(*ctx)->x));
  PetscCall(EPSDestroy(&(*ctx)->eps));

  PetscCall(PetscFree(*ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  DGPContext          *ctx;
  Vec                 H_diag,rho_guess;
  PetscScalar         diagonal_sum,kr;
  PetscInt            n=3,maxit=100,its;
  PetscInt            nconv;
  PetscReal           beta=2.0,rtol=SLEPC_DEFAULT_TOL,stol=1e-12;
  SNES                snes;
  Vec                 F; // Vector to store the residual (F = rho_out - rho_in)
  PetscBool           symm=PETSC_TRUE,verbose=PETSC_FALSE;
  SNESConvergedReason reason;
  SNESLineSearch      linesearch;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCheck(PetscDefined(USE_COMPLEX),PETSC_COMM_WORLD,PETSC_ERR_SUP,"This example requires complex scalars");

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--- DGP SNES SCF ---\n"));

  /* 1. Read DGP-specific command line options */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-beta",&beta,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-symm",&symm,NULL));

  /* 2. Initialization of the DGP context and base linear Hamiltonian A0 */
  PetscCall(DGPCreate(PETSC_COMM_WORLD,n,beta,symm,&ctx));

  /* 3. Configure objects for the SCF (H matrix, EPS, Vec x) */
  PetscCall(DGPSetupSCF(ctx,rtol*0.1));

  /* 4. Generate the initial guess by solving the linear problem A0 */
  PetscCall(DGPGenerateInitialGuess(ctx));

  /* 5. Prepare the initial density (rho_guess) from the linear ground state */
  PetscCall(VecDuplicate(ctx->rho,&rho_guess));
  PetscCall(DGPCalculateDensity(ctx,rho_guess));

  /* 6. Prepare the residual vector F by cloning rho */
  PetscCall(VecDuplicate(ctx->rho,&F));

  /* 7. Configure the nonlinear solver engine (SNES) with the DGP callback */
  PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
  PetscCall(SNESSetFunction(snes,F,DGPIterationSCF,ctx));
  PetscCall(SNESSetTolerances(snes,PETSC_DETERMINE,rtol,stol,maxit,PETSC_DETERMINE));

  // Default -> NRICHARDSON with step lambda = 1.0 to simulate basic SCF
  PetscCall(SNESSetType(snes,SNESNRICHARDSON));
  PetscCall(SNESGetLineSearch(snes,&linesearch));
  PetscCall(SNESLineSearchSetType(linesearch,SNESLINESEARCHNONE));
  PetscCall(SNESLineSearchSetDamping(linesearch,1.0));

  PetscCall(SNESSetFromOptions(snes));

  PetscCall(SNESGetTolerances(snes,NULL,&rtol,NULL,&maxit,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Current parameters: n=%" PetscInt_FMT ", beta=%g, symm=%s, rtol=%g, maxit=%" PetscInt_FMT "\n",n,(double)beta,symm ? "True" : "False",(double)rtol,maxit));

  /* 8. SCF loop */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nStarting SCF iterations...\n"));
  PetscCall(SNESSolve(snes,NULL,rho_guess));

  PetscCall(SNESGetConvergedReason(snes,&reason));
  PetscCall(SNESGetIterationNumber(snes,&its));

  /* 9. Analysis of results */
  if (reason>0) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--> CONVERGENCE REACHED in %" PetscInt_FMT " iterations (Reason: %s).\n",its,SNESConvergedReasons[reason]));

    /* Extract the eigenvalue of the converged ground state */
    PetscCall(EPSGetConverged(ctx->eps,&nconv));
    if (nconv>0) {
      PetscCall(EPSGetEigenvalue(ctx->eps,0,&kr,NULL));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nEigenvalue found = %10.6f\n",(double)PetscRealPart(kr)));
    }

    PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));
    if (verbose) {
      PetscCall(MatCreateVecs(ctx->H,NULL,&H_diag));
      PetscCall(MatGetDiagonal(ctx->H,H_diag));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Diagonal of the final Hamiltonian:\n"));
      PetscCall(VecView(H_diag,PETSC_VIEWER_STDOUT_WORLD));

      PetscCall(VecSum(H_diag,&diagonal_sum));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Sum of the diagonal (Approx total energy): %g\n",(double)PetscRealPart(diagonal_sum)));
      PetscCall(VecDestroy(&H_diag));
    }

  } else {
    PetscCall(SNESGetConvergedReason(snes,&reason));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--> ERROR: The SCF did not converge (Reason: %s)\n",SNESConvergedReasons[reason]));
  }

  /* 10. Clean up memory */
  PetscCall(VecDestroy(&F));
  PetscCall(VecDestroy(&rho_guess));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DGPDestroy(&ctx));

  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: complex

   testset:
      filter: sed -e "s/rtol=1e-05/rtol=1e-08/" -e "s/rtol=1e-16/rtol=1e-08/" -e "s/[0-9]\{1,\} iterations/23 iterations/" -e "s/CONVERGED_SNORM_RELATIVE/CONVERGED_FNORM_RELATIVE/"
      output_file: output/ex60_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -snes_type anderson -snes_anderson_m 7
      test:
         suffix: 3
         args: -snes_type ngmres -npc_snes_type nrichardson -snes_npc_side right
      test:
         suffix: 4
         args: -snes_type composite -snes_composite_type additiveoptimal -snes_composite_sneses anderson,nrichardson -sub_0_snes_anderson_m 7 -sub_0_snes_anderson_beta 0.4

TEST*/
