/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "1-D discrete Kohn-Sham model solved as a Nonlinear Eigenvalue Problem with eigenvector dependence (NEPv).\n\n"
  "Implemented by embedding an EPS inside a nonlinear SNES loop.\n\n"
  "The Kohn-Sham (KS) equations simplify a complex, interacting many-electron system by mapping it\n"
  "onto an equivalent system of noninteracting particles. The goal is to find the single-particle\n"
  "wavefunctions (eigenvectors) and orbital energies (eigenvalues) that self-consistently describe the ground state.\n\n"
  "The nonlinear Hamiltonian matrix is defined as H(X) = L + alpha * Diag(L^-1 * rho(X)), where:\n"
  "   L   = 1-D discrete Laplacian operator matrix.\n"
  "   rho = electronic density vector computed from the current block of eigenvectors X as rho = diag(X * X').\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid points.\n"
  "  -k <k>, where <k> = number of eigenvectors to compute.\n"
  "  -alpha <alpha>, where <alpha> = real scaling parameter controlling the nonlinearity.\n"
  "  -verbose, to print the converged Hamiltonian diagonal and its sum.\n\n";

#include <slepceps.h>

/*
  Context structure to store the necessary data regarding the discrete Kohn-Sham (DKS) problem
  as well as the objects needed for the self-consistent field (SCF) iteration.
*/
typedef struct {
  // DKS context
  PetscInt  n;       // Spatial mesh size
  PetscInt  k;       // Number of eigenvectors to compute
  PetscReal alpha;   // Parameter controlling the nonlinearity
  Mat       L;       // 1D discrete Laplacian matrix
  KSP       ksp;     // Linear solver (to apply L^-1)

  Vec       rho;     // Vector to store the electronic density
  Vec       z;       // Intermediate vector to store the result z = L^-1 * rho

  // SCF context
  EPS       eps;     // Solver for each H
  Mat       H;       // Hamiltonian matrix
  BV        X;       // Block of eigenvectors
} DKSContext;

/*
  Initialize the context for the 1D discrete Kohn-Sham problem.

  This function allocates the necessary memory, builds the 1D Laplacian operator,
  and configures the linear solver (KSP) to apply L^-1.

  Arguments:
    comm    - MPI communicator
    n       - Spatial mesh size
    k       - Number of eigenvectors to compute
    alpha   - Parameter controlling the nonlinearity
    ctx_out - Output pointer where the created context will be stored
*/
PetscErrorCode DKSCreate(MPI_Comm comm,PetscInt n,PetscInt k,PetscReal alpha,DKSContext **ctx_out)
{
  DKSContext *ctx;
  PetscInt   Istart,Iend,i;
  PC         pc;

  PetscFunctionBeginUser;
  PetscCall(PetscNew(&ctx));
  ctx->n=n;
  ctx->k=k;
  ctx->alpha=alpha;

  // 1. Create and fill the Laplacian L
  PetscCall(MatCreate(comm,&ctx->L));
  PetscCall(MatSetSizes(ctx->L,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(ctx->L));

  PetscCall(MatGetOwnershipRange(ctx->L,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    PetscCall(MatSetValue(ctx->L,i,i,2.0,INSERT_VALUES));
    if (i>0) PetscCall(MatSetValue(ctx->L,i,i-1,-1.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(ctx->L,i,i+1,-1.0,INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(ctx->L,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(ctx->L,MAT_FINAL_ASSEMBLY));

  // 2. Configure the KSP
  PetscCall(KSPCreate(comm,&ctx->ksp));
  PetscCall(KSPSetOperators(ctx->ksp,ctx->L,ctx->L));
  PetscCall(KSPSetType(ctx->ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ctx->ksp,&pc));
  PetscCall(PCSetType(pc,PCCHOLESKY));

  if (PetscDefined(HAVE_MUMPS)) PetscCall(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));

  PetscCall(KSPSetFromOptions(ctx->ksp));

  // 3. Create internal work vectors
  PetscCall(MatCreateVecs(ctx->L,&ctx->rho,&ctx->z));

  *ctx_out=ctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Configure the SLEPc objects necessary for the SCF iteration.

  This function prepares the Hamiltonian matrix by cloning the structure of the Laplacian,
  initializes the basis vectors (BV) block for the eigenvectors, and configures the eigenvalue
  solver (EPS) by defining the problem type and adjusting its tolerance.

  Arguments:
    ctx - Pointer to the previously initialized DKS context
    tol - Desired tolerance for the internal eigenvalue solver
*/
PetscErrorCode DKSSetupSCF(DKSContext *ctx,PetscReal tol)
{
  MPI_Comm comm;

  PetscFunctionBeginUser;
  comm = PetscObjectComm((PetscObject)ctx->L);
  /* Prepare the H matrix */
  PetscCall(MatDuplicate(ctx->L,MAT_DO_NOT_COPY_VALUES,&ctx->H));

  /* Configure the eigenvector block X */
  PetscCall(BVCreate(comm,&ctx->X));
  PetscCall(BVSetSizesFromVec(ctx->X,ctx->rho,ctx->k));
  PetscCall(BVSetFromOptions(ctx->X));

  /* Configure SLEPc EPS */
  PetscCall(EPSCreate(comm,&ctx->eps));
  PetscCall(EPSSetProblemType(ctx->eps,EPS_HEP));
  PetscCall(EPSSetWhichEigenpairs(ctx->eps,EPS_SMALLEST_REAL));
  PetscCall(EPSSetDimensions(ctx->eps,ctx->k,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(EPSSetTolerances(ctx->eps,tol,PETSC_DECIDE));
  PetscCall(EPSSetFromOptions(ctx->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Generate the initial guess X0 using the exact eigenvectors of the 1D Laplacian.

  This function fills the vector block X with the initial guess, which greatly
  improves the convergence of the SCF loop.

  Formula used: vv = [1:n]'/(n+1)*pi; X0 = sin(vv * [1:k])*sqrt(2/(n+1));

  Arguments:
    ctx - Pointer to the previously initialized DKS context
*/
PetscErrorCode DKSGenerateInitialGuess(DKSContext *ctx)
{
  PetscInt    Istart,Iend,i,j;
  PetscReal   h_val,norm_factor,val;
  Vec         col;
  PetscScalar *x_local;

  PetscFunctionBeginUser;
  // Mathematical constants of the formula
  // X0 = sin( [1:n]' * [1:k] * pi/(n+1) ) * sqrt(2/(n+1));
  h_val=PETSC_PI/(ctx->n+1.0);
  norm_factor=PetscSqrtReal(2.0/(ctx->n+1.0));

  // Iterate over each column (eigenvector)
  for (j=0;j<ctx->k;j++) {
    PetscCall(BVGetColumn(ctx->X,j,&col));
    PetscCall(VecGetOwnershipRange(col,&Istart,&Iend));
    PetscCall(VecGetArray(col,&x_local));

    for (i=Istart;i<Iend;i++) {
      // i and j start at 0 in C, so we add 1 for the mathematical formula
      val=PetscSinReal((i+1.0)*(j+1.0)*h_val)*norm_factor;
      x_local[i-Istart]=val;
    }

    PetscCall(VecRestoreArray(col,&x_local));
    PetscCall(BVRestoreColumn(ctx->X,j,&col));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Calculate the electronic density from the eigenvectors.

  This function computes the density rho(X) as the sum of the squares of the
  components of each eigenvector. The formula used is: rho(X) = diag(X * X').

  Arguments:
    ctx     - Pointer to the initialized DKS context
    rho_out - Pre-created vector where the computed density will be stored
*/
PetscErrorCode DKSCalculateDensity(DKSContext *ctx,Vec rho_out)
{
  PetscInt          i,j,n_loc;
  Vec               col;
  const PetscScalar *x_local;
  PetscScalar       *rho_local;

  PetscFunctionBeginUser;
  PetscCall(VecSet(rho_out,0.0));
  PetscCall(VecGetLocalSize(rho_out,&n_loc));
  PetscCall(VecGetArray(rho_out,&rho_local));

  for (j=0;j<ctx->k;j++) {
    PetscCall(BVGetColumn(ctx->X,j,&col));
    PetscCall(VecGetArrayRead(col,&x_local));
    for (i=0;i<n_loc;i++) {
      rho_local[i]+=x_local[i]*PetscConj(x_local[i]);
    }

    PetscCall(VecRestoreArrayRead(col,&x_local));
    PetscCall(BVRestoreColumn(ctx->X,j,&col));
  }

  PetscCall(VecRestoreArray(rho_out,&rho_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Build the Hamiltonian from an input density.

  This function solves the linear system L * z = rho_in using the configured
  KSP solver. Then, it assembles the updated Hamiltonian matrix stored in the context
  using the formula: H = L + alpha * Diag(z), where z = L^-1 * rho_in.

  Arguments:
    ctx    - Pointer to the initialized DKS context (ctx->H and ctx->z are updated)
    rho_in - Vector with the proposed input electronic density
*/
PetscErrorCode DKSBuildHamiltonian(DKSContext *ctx,Vec rho_in)
{
  PetscFunctionBeginUser;
  // 1. Solve L * z = rho_in
  PetscCall(KSPSolve(ctx->ksp,rho_in,ctx->z));

  // 2. Build H = L + alpha * Diag(z)
  PetscCall(MatCopy(ctx->L,ctx->H,SAME_NONZERO_PATTERN));
  PetscCall(VecScale(ctx->z,ctx->alpha));
  PetscCall(MatDiagonalSet(ctx->H,ctx->z,ADD_VALUES));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Evaluation function (Callback) for the SNES nonlinear solver.

  In each SNES iteration, this function receives a proposed density (rho_in),
  builds the Hamiltonian, solves the eigenvalue equation, computes the
  resulting density, and returns the residual F = rho_out - rho_in.

  Arguments:
    snes     - The nonlinear solver context
    rho_in   - Input electronic density proposed by SNES
    F        - Vector where the computed residual will be stored
    ctx_void - Pointer to the user's DKS context (DKSContext)
*/
PetscErrorCode DKSIterationSCF(SNES snes,Vec rho_in,Vec F,void *ctx_void)
{
  DKSContext *ctx=(DKSContext*)ctx_void;
  PetscInt   j,nconv;
  Vec        xr,col;
  MPI_Comm   comm;

  PetscFunctionBeginUser;
  comm = PetscObjectComm((PetscObject)ctx->L);
  /* 1. Build the physics (Hamiltonian) with the guess proposed by SNES */
  PetscCall(DKSBuildHamiltonian(ctx,rho_in));

  /* 2. Solve with the current H */
  // Pass our newly built H matrix to SLEPc
  PetscCall(EPSSetOperators(ctx->eps,ctx->H,NULL));
  PetscCall(EPSSolve(ctx->eps));

  // (Safety check: verify that SLEPc has not failed internally)
  PetscCall(EPSGetConverged(ctx->eps,&nconv));
  if (nconv<ctx->k) PetscCall(PetscPrintf(comm,"Warning: SLEPc only converged %" PetscInt_FMT " out of %" PetscInt_FMT " eigenvalues.\n",nconv,ctx->k));

  /* 3. Extract the new eigenvectors and store them in ctx->X */
  PetscCall(MatCreateVecs(ctx->H,&xr,NULL));

  for (j=0;j<ctx->k;j++) {
    PetscCall(EPSGetEigenvector(ctx->eps,j,xr,NULL)); // Extracts the j-th vector
    PetscCall(BVGetColumn(ctx->X,j,&col));          // Gets column j of X
    PetscCall(VecCopy(xr,col));                     // Copies the data
    PetscCall(BVRestoreColumn(ctx->X,j,&col));      // Restores the column
  }

  /* 4. Calculate the new density generated by these electrons */
  // We use ctx->rho as a temporary work vector to store rho_out
  PetscCall(DKSCalculateDensity(ctx,ctx->rho));

  /* 5. Calculate the nonlinear residual: F = rho_out - rho_in */
  PetscCall(VecWAXPY(F,-1.0,ctx->rho,rho_in));

  /* Cleanup of the temporary memory of this iteration */
  PetscCall(VecDestroy(&xr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Free the memory associated with the DKS context.

  This function destroys all internal PETSc objects created during
  initialization and frees the memory of the main structure.

  Arguments:
    ctx - Pointer to the DKS context (set to NULL upon completion)
*/
PetscErrorCode DKSDestroy(DKSContext **ctx)
{
  PetscFunctionBeginUser;
  if (!*ctx) PetscFunctionReturn(PETSC_SUCCESS);

  PetscCall(MatDestroy(&(*ctx)->L));
  PetscCall(KSPDestroy(&(*ctx)->ksp));
  PetscCall(VecDestroy(&(*ctx)->rho));
  PetscCall(VecDestroy(&(*ctx)->z));

  PetscCall(MatDestroy(&(*ctx)->H));
  PetscCall(BVDestroy(&(*ctx)->X));
  PetscCall(EPSDestroy(&(*ctx)->eps));

  PetscCall(PetscFree(*ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc,char **argv)
{
  DKSContext          *ctx;
  Vec                 H_diag,rho_guess;
  PetscScalar         diagonal_sum,kr;
  PetscInt            n=5,k=3,maxit=100,its;
  PetscInt            nconv;
  PetscReal           alpha=0.5,rtol=SLEPC_DEFAULT_TOL,stol=1e-12;
  SNES                snes;
  Vec                 F; // Vector to store the residual (F = rho_out - rho_in)
  PetscBool           verbose=PETSC_FALSE;
  SNESConvergedReason reason;
  SNESLineSearch      linesearch;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--- DKS SNES SCF ---\n"));

  /* 1. Read options from the command line (if provided by the user) */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-alpha",&alpha,NULL));

  /* 2. Initialization of the DKS context */
  PetscCall(DKSCreate(PETSC_COMM_WORLD,n,k,alpha,&ctx));

  /* 3. Configure general objects for the SCF (Matrices, EPS, eigenvectors) */
  PetscCall(DKSSetupSCF(ctx,rtol*0.1));

  /* 4. Generate the initial guess for the eigenvectors and their density */
  PetscCall(DKSGenerateInitialGuess(ctx));

  /* 5. Calculate the initial density (rho_guess) from X0 */
  // Since SNES works with densities, the density is extracted from our X0
  PetscCall(VecDuplicate(ctx->rho,&rho_guess));
  PetscCall(DKSCalculateDensity(ctx,rho_guess));

  /* 6. Prepare the residual vector F by cloning the structure of rho */
  PetscCall(VecDuplicate(ctx->rho,&F));

  /* 7. Create and set up the nonlinear solver engine (SNES) */
  PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
  PetscCall(SNESSetFunction(snes,F,DKSIterationSCF,ctx));
  PetscCall(SNESSetTolerances(snes,PETSC_DETERMINE,rtol,stol,maxit,PETSC_DETERMINE));

  // Default -> NRICHARDSON with step lambda = 1.0 to simulate basic SCF
  PetscCall(SNESSetType(snes,SNESNRICHARDSON));
  PetscCall(SNESGetLineSearch(snes,&linesearch));
  PetscCall(SNESLineSearchSetType(linesearch,SNESLINESEARCHNONE));
  PetscCall(SNESLineSearchSetDamping(linesearch,1.0));

  PetscCall(SNESSetFromOptions(snes));

  PetscCall(SNESGetTolerances(snes,NULL,&rtol,NULL,&maxit,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Current parameters: n=%" PetscInt_FMT ", k=%" PetscInt_FMT ", alpha=%g, rtol=%g, maxit=%" PetscInt_FMT "\n",n,k,(double)alpha,(double)rtol,maxit));

  /* 8. SCF Loop */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nStarting SCF iterations...\n"));
  PetscCall(SNESSolve(snes,NULL,rho_guess));

  PetscCall(SNESGetConvergedReason(snes,&reason));
  PetscCall(SNESGetIterationNumber(snes,&its));

  /* 9. Result analysis */
  if (reason>0) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--> CONVERGENCE REACHED in %" PetscInt_FMT " iterations (Reason: %s).\n",its,SNESConvergedReasons[reason]));

    /* Show the final result */

    PetscCall(EPSGetConverged(ctx->eps,&nconv));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n--- Eigenvalues (%" PetscInt_FMT " found) ---\n",nconv));

    for (PetscInt i=0;i<nconv;i++) {
      PetscCall(EPSGetEigenvalue(ctx->eps,i,&kr,NULL));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Eigenvalue[%" PetscInt_FMT "] = %10.6f\n",i,(double)PetscRealPart(kr)));
    }

    PetscCall(PetscOptionsHasName(NULL,NULL,"-verbose",&verbose));

    if (verbose) {
      PetscCall(MatCreateVecs(ctx->H,NULL,&H_diag));
      PetscCall(MatGetDiagonal(ctx->H,H_diag));

      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Diagonal of the converged Hamiltonian:\n"));
      PetscCall(VecView(H_diag,PETSC_VIEWER_STDOUT_WORLD));

      PetscCall(VecSum(H_diag,&diagonal_sum));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Sum of the diagonal: %g\n",(double)PetscRealPart(diagonal_sum)));
      PetscCall(VecDestroy(&H_diag));
    }

  } else {
    PetscCall(SNESGetConvergedReason(snes,&reason));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"--> ERROR: SCF did not converge (Reason: %s)\n",SNESConvergedReasons[reason]));
  }

  /* 10. Clean up memory */
  PetscCall(VecDestroy(&F));
  PetscCall(VecDestroy(&rho_guess));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DKSDestroy(&ctx));

  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      filter: sed -e "s/1.364212/1.364211/" -e "s/4.864809/4.864808/" -e "s/1.921852/1.921853/" -e "s/1.921854/1.921853/" -e "s/2.931104/2.931103/" -e "s/3.957516/3.957515/" -e "s/rtol=1e-05/rtol=1e-08/" -e "s/rtol=1e-16/rtol=1e-08/" -e "s/[0-9]\{1,\} iterations/8 iterations/" -e "s/CONVERGED_SNORM_RELATIVE/CONVERGED_FNORM_RELATIVE/"
      output_file: output/ex59_1.out
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
