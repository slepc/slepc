/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This is a nonlinear eigenvalue problem. When p=2, it is reduced to a linear Laplace eigenvalue
   problem.

   -\nabla\cdot(|\nabla u|^{p-2} \nabla u) = k |u|^{p-2} u in (0,1)x(0,1),

   u = 0 on the entire boundary.

   The code is implemented based on DMPlex using Q1 FEM on a quadrilateral mesh. In this code, we consider p=3.

   Contributed  by Fande Kong fdkong.jd@gmail.com
*/

static char help[] = "Nonlinear inverse iteration for A(x)*x=lambda*B(x)*x.\n\n";

#include <slepceps.h>
#include <petscdmplex.h>
#include <petscds.h>

PetscErrorCode CreateSquareMesh(MPI_Comm,DM*);
PetscErrorCode SetupDiscretization(DM);
PetscErrorCode FormJacobianA(SNES,Vec,Mat,Mat,void*);
PetscErrorCode FormFunctionA(SNES,Vec,Vec,void*);
PetscErrorCode MatMult_A(Mat A,Vec x,Vec y);
PetscErrorCode FormJacobianB(SNES,Vec,Mat,Mat,void*);
PetscErrorCode FormFunctionB(SNES,Vec,Vec,void*);
PetscErrorCode MatMult_B(Mat A,Vec x,Vec y);
PetscErrorCode FormFunctionAB(SNES,Vec,Vec,Vec,void*);
PetscErrorCode BoundaryGlobalIndex(DM,const char*,IS*);
PetscErrorCode FormNorm(SNES,Vec,PetscReal*,void*);

typedef struct {
  IS    bdis; /* global indices for boundary DoFs */
  SNES  snes;
  EPS   eps;
} AppCtx;

int main(int argc,char **argv)
{
  DM             dm;
  MPI_Comm       comm;
  AppCtx         user;
  EPS            eps;  /* eigenproblem solver context */
  ST             st;
  EPSType        type;
  Mat            A,B,P;
  Vec            v0;
  PetscContainer container;
  PetscInt       nev,nconv,m,n,M,N;
  PetscBool      nonlin,flg=PETSC_FALSE,update;
  SNES           snes;
  PetscReal      tol,relerr;
  PetscBool      use_shell_matrix=PETSC_FALSE,test_init_sol=PETSC_FALSE,use_custom_norm=PETSC_FALSE,sign_normalization=PETSC_TRUE;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  comm = PETSC_COMM_WORLD;
  /* Create a quadrilateral mesh on domain (0,1)x(0,1) */
  PetscCall(CreateSquareMesh(comm,&dm));
  /* Setup basis function */
  PetscCall(SetupDiscretization(dm));
  PetscCall(BoundaryGlobalIndex(dm,"marker",&user.bdis));
  /* Check if we are going to use shell matrices */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-use_shell_matrix",&use_shell_matrix,NULL));
  if (use_shell_matrix) {
    PetscCall(DMCreateMatrix(dm,&P));
    PetscCall(MatGetLocalSize(P,&m,&n));
    PetscCall(MatGetSize(P,&M,&N));
    PetscCall(MatCreateShell(comm,m,n,M,N,&user,&A));
    PetscCall(MatShellSetOperation(A,MATOP_MULT,(PetscErrorCodeFn*)MatMult_A));
    PetscCall(MatCreateShell(comm,m,n,M,N,&user,&B));
    PetscCall(MatShellSetOperation(B,MATOP_MULT,(PetscErrorCodeFn*)MatMult_B));
  } else {
    PetscCall(DMCreateMatrix(dm,&A));
    PetscCall(MatDuplicate(A,MAT_COPY_VALUES,&B));
  }
  /* Check whether we should use a custom normalization */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-use_custom_norm",&use_custom_norm,NULL));
  /* Check whether we should normalize Bx by the sign of its first nonzero element */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-sign_normalization",&sign_normalization,NULL));

  /*
     Compose callback functions and context that will be needed by the solver
  */
  PetscCall(PetscObjectComposeFunction((PetscObject)A,"formFunction",FormFunctionA));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-form_function_ab",&flg,NULL));
  if (flg) PetscCall(PetscObjectComposeFunction((PetscObject)A,"formFunctionAB",FormFunctionAB));
  PetscCall(PetscObjectComposeFunction((PetscObject)A,"formJacobian",FormJacobianA));
  PetscCall(PetscObjectComposeFunction((PetscObject)B,"formFunction",FormFunctionB));
  if (use_custom_norm) PetscCall(PetscObjectComposeFunction((PetscObject)B,"formNorm",FormNorm));
  PetscCall(PetscContainerCreate(comm,&container));
  PetscCall(PetscContainerSetPointer(container,&user));
  PetscCall(PetscObjectCompose((PetscObject)A,"formFunctionCtx",(PetscObject)container));
  PetscCall(PetscObjectCompose((PetscObject)A,"formJacobianCtx",(PetscObject)container));
  PetscCall(PetscObjectCompose((PetscObject)B,"formFunctionCtx",(PetscObject)container));
  if (use_custom_norm) PetscCall(PetscObjectCompose((PetscObject)B,"formNormCtx",(PetscObject)container));
  PetscCall(PetscContainerDestroy(&container));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(comm,&eps));
  PetscCall(EPSSetOperators(eps,A,B));
  PetscCall(EPSSetProblemType(eps,EPS_GNHEP));
  user.eps = eps;
  /*
     Use nonlinear inverse iteration
  */
  PetscCall(EPSSetType(eps,EPSPOWER));
  PetscCall(EPSPowerSetNonlinear(eps,PETSC_TRUE));
  /* Set the Bx sign normalization (or not) */
  PetscCall(EPSPowerSetSignNormalization(eps,sign_normalization));
  /*
    Attach DM to SNES
  */
  PetscCall(EPSPowerGetSNES(eps,&snes));
  user.snes = snes;
  PetscCall(SNESSetDM(snes,dm));
  PetscCall(EPSSetFromOptions(eps));

  /* Set a preconditioning matrix to ST */
  if (use_shell_matrix) {
    PetscCall(EPSGetST(eps,&st));
    PetscCall(STSetPreconditionerMat(st,P));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-test_init_sol",&test_init_sol,NULL));
  if (nconv && test_init_sol) {
    PetscScalar   k;
    PetscReal     norm0;
    PetscInt      nits;

    PetscCall(MatCreateVecs(A,&v0,NULL));
    PetscCall(EPSGetEigenpair(eps,0,&k,NULL,v0,NULL));
    PetscCall(EPSSetInitialSpace(eps,1,&v0));
    PetscCall(VecDestroy(&v0));
    /* Norm of the previous residual */
    PetscCall(SNESGetFunctionNorm(snes,&norm0));
    /* Make the tolerance smaller than the last residual
       SNES will converge right away if the initial is setup correctly */
    PetscCall(SNESSetTolerances(snes,norm0*1.2,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT,PETSC_CURRENT));
    PetscCall(EPSSolve(eps));
    /* Number of Newton iterations supposes to be zero */
    PetscCall(SNESGetIterationNumber(snes,&nits));
    if (nits) PetscCall(PetscPrintf(comm," Number of Newton iterations %" PetscInt_FMT " should be zero \n",nits));
  }

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetType(eps,&type));
  PetscCall(EPSGetTolerances(eps,&tol,NULL));
  PetscCall(EPSPowerGetNonlinear(eps,&nonlin));
  PetscCall(EPSPowerGetUpdate(eps,&update));
  PetscCall(PetscPrintf(comm," Solution method: %s%s\n\n",type,nonlin?(update?" (nonlinear with monolithic update)":" (nonlinear)"):""));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(comm," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* print eigenvalue and error */
  PetscCall(EPSGetConverged(eps,&nconv));
  if (nconv>0) {
    PetscScalar   k;
    PetscReal     na,nb;
    Vec           a,b,eigen;
    PetscCall(DMCreateGlobalVector(dm,&a));
    PetscCall(VecDuplicate(a,&b));
    PetscCall(VecDuplicate(a,&eigen));
    PetscCall(EPSGetEigenpair(eps,0,&k,NULL,eigen,NULL));
    PetscCall(FormFunctionA(snes,eigen,a,&user));
    PetscCall(FormFunctionB(snes,eigen,b,&user));
    PetscCall(VecAXPY(a,-k,b));
    PetscCall(VecNorm(a,NORM_2,&na));
    if (use_custom_norm) PetscCall(FormNorm(snes,b,&nb,&user));
    else PetscCall(VecNorm(b,NORM_2,&nb));
    relerr = na/(nb*PetscAbsScalar(k));
    if (relerr<10*tol) PetscCall(PetscPrintf(comm,"k: %g, relative error below tol\n",(double)PetscRealPart(k)));
    else PetscCall(PetscPrintf(comm,"k: %g, relative error: %g\n",(double)PetscRealPart(k),(double)relerr));
    PetscCall(VecDestroy(&a));
    PetscCall(VecDestroy(&b));
    PetscCall(VecDestroy(&eigen));
  } else PetscCall(PetscPrintf(comm,"Solver did not converge\n"));

  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  if (use_shell_matrix) PetscCall(MatDestroy(&P));
  PetscCall(DMDestroy(&dm));
  PetscCall(EPSDestroy(&eps));
  PetscCall(ISDestroy(&user.bdis));
  PetscCall(SlepcFinalize());
  return 0;
}

/* <|u|u, v> */
static void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscScalar cof = PetscAbsScalar(u[0]);

  f0[0] = cof*u[0];
}

/* <|\nabla u| \nabla u, \nabla v> */
static void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt    d;
  PetscScalar cof = 0;
  for (d = 0; d < dim; ++d)  cof += u_x[d]*u_x[d];

  cof = PetscSqrtScalar(cof);

  for (d = 0; d < dim; ++d) f1[d] = u_x[d]*cof;
}

/* approximate  Jacobian for   <|u|u, v> */
static void g0_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = 1.0*PetscAbsScalar(u[0]);
}

/* approximate  Jacobian for   <|\nabla u| \nabla u, \nabla v> */
static void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;

  for (d = 0; d < dim; ++d) g3[d*dim+d] = 1.0;
}

PetscErrorCode SetupDiscretization(DM dm)
{
  PetscFE        fe;
  MPI_Comm       comm;

  PetscFunctionBeginUser;
  /* Create finite element */
  PetscCall(PetscObjectGetComm((PetscObject)dm,&comm));
  PetscCall(PetscFECreateDefault(comm,2,1,PETSC_FALSE,NULL,-1,&fe));
  PetscCall(PetscObjectSetName((PetscObject)fe,"u"));
  PetscCall(DMSetField(dm,0,NULL,(PetscObject)fe));
  PetscCall(DMCreateDS(dm));
  PetscCall(PetscFEDestroy(&fe));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CreateSquareMesh(MPI_Comm comm,DM *dm)
{
  PetscInt       cells[] = {8,8};
  PetscInt       dim = 2;
  DM             pdm;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscCall(DMPlexCreateBoxMesh(comm,dim,PETSC_FALSE,cells,NULL,NULL,NULL,PETSC_TRUE,0,PETSC_TRUE,dm));
  PetscCall(DMSetFromOptions(*dm));
  PetscCall(DMSetUp(*dm));
  PetscCallMPI(MPI_Comm_size(comm,&size));
  if (size > 1) {
    PetscCall(DMPlexDistribute(*dm,0,NULL,&pdm));
    PetscCall(DMDestroy(dm));
    *dm = pdm;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode BoundaryGlobalIndex(DM dm,const char labelname[],IS *bdis)
{
  IS             bdpoints;
  PetscInt       nindices,*indices,numDof,offset,npoints,i,j;
  const PetscInt *bdpoints_indices;
  DMLabel        bdmarker;
  PetscSection   gsection;

  PetscFunctionBegin;
  PetscCall(DMGetGlobalSection(dm,&gsection));
  PetscCall(DMGetLabel(dm,labelname,&bdmarker));
  PetscCall(DMLabelGetStratumIS(bdmarker,1,&bdpoints));
  PetscCall(ISGetLocalSize(bdpoints,&npoints));
  PetscCall(ISGetIndices(bdpoints,&bdpoints_indices));
  nindices = 0;
  for (i=0;i<npoints;i++) {
    PetscCall(PetscSectionGetDof(gsection,bdpoints_indices[i],&numDof));
    if (numDof<=0) continue;
    nindices += numDof;
  }
  PetscCall(PetscCalloc1(nindices,&indices));
  nindices = 0;
  for (i=0;i<npoints;i++) {
    PetscCall(PetscSectionGetDof(gsection,bdpoints_indices[i],&numDof));
    if (numDof<=0) continue;
    PetscCall(PetscSectionGetOffset(gsection,bdpoints_indices[i],&offset));
    for (j=0;j<numDof;j++) indices[nindices++] = offset+j;
  }
  PetscCall(ISRestoreIndices(bdpoints,&bdpoints_indices));
  PetscCall(ISDestroy(&bdpoints));
  PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)dm),nindices,indices,PETSC_OWN_POINTER,bdis));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FormJacobian(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  DM             dm;
  Vec            Xloc;

  PetscFunctionBegin;
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetLocalVector(dm,&Xloc));
  PetscCall(VecZeroEntries(Xloc));
  PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,Xloc));
  PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,Xloc));
  CHKMEMQ;
  PetscCall(DMPlexSNESComputeJacobianFEM(dm,Xloc,A,B,ctx));
  if (A!=B) {
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  }
  CHKMEMQ;
  PetscCall(DMRestoreLocalVector(dm,&Xloc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormJacobianA(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  DM             dm;
  PetscDS        prob;
  PetscWeakForm  wf;
  AppCtx         *userctx = (AppCtx *)ctx;

  PetscFunctionBegin;
  PetscCall(MatSetOption(B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE));
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetDS(dm,&prob));
  PetscCall(PetscDSGetWeakForm(prob, &wf));
  PetscCall(PetscWeakFormClearIndex(wf, NULL, 0, 0, 0, PETSC_WF_G3, 0));
  PetscCall(PetscWeakFormSetIndexJacobian(wf, NULL, 0, 0, 0, 0, 0, NULL, 0, NULL, 0, NULL, 0, g3_uu));
  PetscCall(FormJacobian(snes,X,A,B,ctx));
  PetscCall(MatZeroRowsIS(B,userctx->bdis,1.0,NULL,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormJacobianB(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  DM             dm;
  PetscDS        prob;
  PetscWeakForm  wf;
  AppCtx         *userctx = (AppCtx *)ctx;

  PetscFunctionBegin;
  PetscCall(MatSetOption(B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE));
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetDS(dm,&prob));
  PetscCall(PetscDSGetWeakForm(prob, &wf));
  PetscCall(PetscWeakFormClearIndex(wf, NULL, 0, 0, 0, PETSC_WF_G3, 0));
  PetscCall(PetscWeakFormSetIndexJacobian(wf, NULL, 0, 0, 0, 0, 0, g0_uu, 0, NULL, 0, NULL, 0, NULL));
  PetscCall(FormJacobian(snes,X,A,B,ctx));
  PetscCall(MatZeroRowsIS(B,userctx->bdis,0.0,NULL,NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormFunctionAB(SNES snes,Vec x,Vec Ax,Vec Bx,void *ctx)
{
  PetscFunctionBegin;
  /*
   * In real applications, users should have a generic formFunctionAB which
   * forms Ax and Bx simultaneously for a more efficient calculation.
   * In this example, we just call FormFunctionA+FormFunctionB to mimic how
   * to use FormFunctionAB
   */
  PetscCall(FormFunctionA(snes,x,Ax,ctx));
  PetscCall(FormFunctionB(snes,x,Bx,ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FormFunction(SNES snes,Vec X,Vec F,void *ctx)
{
  DM             dm;
  Vec            Xloc,Floc;

  PetscFunctionBegin;
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetLocalVector(dm,&Xloc));
  PetscCall(DMGetLocalVector(dm,&Floc));
  PetscCall(VecZeroEntries(Xloc));
  PetscCall(VecZeroEntries(Floc));
  PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,Xloc));
  PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,Xloc));
  CHKMEMQ;
  PetscCall(DMPlexSNESComputeResidualFEM(dm,Xloc,Floc,ctx));
  CHKMEMQ;
  PetscCall(VecZeroEntries(F));
  PetscCall(DMLocalToGlobalBegin(dm,Floc,ADD_VALUES,F));
  PetscCall(DMLocalToGlobalEnd(dm,Floc,ADD_VALUES,F));
  PetscCall(DMRestoreLocalVector(dm,&Xloc));
  PetscCall(DMRestoreLocalVector(dm,&Floc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormFunctionA(SNES snes,Vec X,Vec F,void *ctx)
{
  DM             dm;
  PetscDS        prob;
  PetscWeakForm  wf;
  PetscInt       nindices,iStart,iEnd,i;
  AppCtx         *userctx = (AppCtx *)ctx;
  PetscScalar    *array,value;
  const PetscInt *indices;
  PetscInt       vecstate;

  PetscFunctionBegin;
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetDS(dm,&prob));
  /* hook functions */
  PetscCall(PetscDSGetWeakForm(prob, &wf));
  PetscCall(PetscWeakFormClearIndex(wf, NULL, 0, 0, 0, PETSC_WF_F0, 0));
  PetscCall(PetscWeakFormSetIndexResidual(wf, NULL, 0, 0, 0, 0, NULL, 0, f1_u));
  PetscCall(FormFunction(snes,X,F,ctx));
  /* Boundary condition */
  PetscCall(VecLockGet(X,&vecstate));
  if (vecstate>0) PetscCall(VecLockReadPop(X));
  PetscCall(VecGetOwnershipRange(X,&iStart,&iEnd));
  PetscCall(VecGetArray(X,&array));
  PetscCall(ISGetLocalSize(userctx->bdis,&nindices));
  PetscCall(ISGetIndices(userctx->bdis,&indices));
  for (i=0;i<nindices;i++) {
    value = array[indices[i]-iStart] - 0.0;
    PetscCall(VecSetValue(F,indices[i],value,INSERT_VALUES));
  }
  PetscCall(ISRestoreIndices(userctx->bdis,&indices));
  PetscCall(VecRestoreArray(X,&array));
  if (vecstate>0) PetscCall(VecLockReadPush(X));
  PetscCall(VecAssemblyBegin(F));
  PetscCall(VecAssemblyEnd(F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormNorm(SNES snes,Vec Bx,PetscReal *norm,void *ctx)
{
  PetscFunctionBegin;
  PetscCall(VecNorm(Bx,NORM_2,norm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_A(Mat A,Vec x,Vec y)
{
  AppCtx         *userctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(A,&userctx));
  PetscCall(FormFunctionA(userctx->snes,x,y,userctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormFunctionB(SNES snes,Vec X,Vec F,void *ctx)
{
  DM             dm;
  PetscDS        prob;
  PetscWeakForm  wf;
  PetscInt       nindices,iStart,iEnd,i;
  AppCtx         *userctx = (AppCtx *)ctx;
  PetscScalar    value;
  const PetscInt *indices;

  PetscFunctionBegin;
  PetscCall(SNESGetDM(snes,&dm));
  PetscCall(DMGetDS(dm,&prob));
  /* hook functions */
  PetscCall(PetscDSGetWeakForm(prob, &wf));
  PetscCall(PetscWeakFormClearIndex(wf, NULL, 0, 0, 0, PETSC_WF_F1, 0));
  PetscCall(PetscWeakFormSetIndexResidual(wf, NULL, 0, 0, 0, 0, f0_u, 0, NULL));
  PetscCall(FormFunction(snes,X,F,ctx));
  /* Boundary condition */
  PetscCall(VecGetOwnershipRange(F,&iStart,&iEnd));
  PetscCall(ISGetLocalSize(userctx->bdis,&nindices));
  PetscCall(ISGetIndices(userctx->bdis,&indices));
  for (i=0;i<nindices;i++) {
    value = 0.0;
    PetscCall(VecSetValue(F,indices[i],value,INSERT_VALUES));
  }
  PetscCall(ISRestoreIndices(userctx->bdis,&indices));
  PetscCall(VecAssemblyBegin(F));
  PetscCall(VecAssemblyEnd(F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMult_B(Mat B,Vec x,Vec y)
{
  AppCtx         *userctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&userctx));
  PetscCall(FormFunctionB(userctx->snes,x,y,userctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      requires: double
      args: -petscspace_degree 1 -petscspace_poly_tensor -checkfunctionlist 0
      output_file: output/ex34_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -eps_power_update -form_function_ab {{0 1}}
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 3
         args: -use_shell_matrix -eps_power_snes_mf_operator 1
      test:
         suffix: 4
         args: -use_shell_matrix -eps_power_update -init_eps_power_snes_mf_operator 1 -eps_power_snes_mf_operator 1 -form_function_ab {{0 1}}
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 5
         args: -use_shell_matrix -eps_power_update -init_eps_power_snes_mf_operator 1 -eps_power_snes_mf_operator 1 -form_function_ab {{0 1}} -test_init_sol 1
         filter: sed -e "s/ with monolithic update//"

      test:
         suffix: 6
         args: -use_shell_matrix -eps_power_update -init_eps_power_snes_mf_operator 1 -eps_power_snes_mf_operator 1 -form_function_ab {{0 1}} -eps_monitor_all
         output_file: output/ex34_6.out
         filter: sed -e "s/\([+-].*i\)//g" -e "1,3s/[0-9]//g" -e "/[45] EPS/d"
      test:
         suffix: 7
         args: -use_custom_norm -sign_normalization 1 -eps_power_snes_mf_operator 1
      test:
         suffix: 8
         args: -use_custom_norm -sign_normalization 1 -eps_power_update -form_function_ab {{0 1}} -eps_power_snes_mf_operator 1 -init_eps_power_snes_mf_operator 1
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 9
         requires: !complex
         args: -use_custom_norm {{0 1}} -sign_normalization 0 -eps_power_snes_mf_operator 1
      test:
         suffix: 10
         requires: !complex
         args: -use_custom_norm {{0 1}} -sign_normalization 0 -eps_power_update -form_function_ab {{0 1}} -eps_power_snes_mf_operator 1 -init_eps_power_snes_mf_operator 1
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 11
         requires: complex
         args: -use_custom_norm {{0 1}} -sign_normalization 0 -eps_power_snes_type nrichardson -eps_power_snes_atol 1e-12 -eps_power_snes_linesearch_atol 1e-22
      test:
         suffix: 12
         requires: complex
         args: -use_custom_norm {{0 1}} -sign_normalization 0 -eps_power_update -init_eps_power_snes_type nrichardson -init_eps_max_it 2 -eps_power_snes_mf_operator 1
         filter: sed -e "s/ with monolithic update//"
TEST*/
