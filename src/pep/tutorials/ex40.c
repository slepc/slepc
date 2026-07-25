/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Checking the definite property in quadratic symmetric eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices.\n"
  "  -transform... whether to transform to a hyperbolic problem or not.\n"
  "  -nonhyperbolic... to test with a modified (definite) problem that is not hyperbolic.\n\n";

#include <slepcpep.h>

/*
  This example is based on spring.c, for fixed values mu=1,tau=10,kappa=5

  The transformations are based on the method proposed in [Niendorf and Voss, LAA 2010].
*/

PetscErrorCode QEPDefiniteTransformGetMatrices(PEP,PetscBool,PetscReal,PetscReal,Mat[3]);
PetscErrorCode QEPDefiniteTransformMap(PetscBool,PetscReal,PetscReal,PetscInt,PetscScalar*,PetscBool);
PetscErrorCode QEPDefiniteCheckError(Mat*,PEP,PetscBool,PetscReal,PetscReal);
PetscErrorCode TransformMatricesMoebius(Mat[3],MatStructure,PetscReal,PetscReal,PetscReal,PetscReal,Mat[3]);

int main(int argc,char **argv)
{
  Mat            M,C,K,*Op,A[3],At[3],B[3]; /* problem matrices */
  PEP            pep;        /* polynomial eigenproblem solver context */
  ST             st;         /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PEPProblemType type;
  PetscBool      terse,transform=PETSC_FALSE,nohyp=PETSC_FALSE;
  PetscInt       n=100,Istart,Iend,i,def=0,hyp;
  PetscReal      muu=1,tau=10,kappa=5,inta,intb;
  PetscReal      alpha,beta,xi,mu,at[2]={0.0,0.0},c=.857,s;
  PetscScalar    target,targett,ats[2];

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nPEP example that checks definite property, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&K));
  PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(K));

  PetscCall(MatGetOwnershipRange(K,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(K,i,i-1,-kappa,INSERT_VALUES));
    PetscCall(MatSetValue(K,i,i,kappa*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(K,i,i+1,-kappa,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

  /* C is a tridiagonal */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&C));
  PetscCall(MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(C));

  PetscCall(MatGetOwnershipRange(C,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i>0) PetscCall(MatSetValue(C,i,i-1,-tau,INSERT_VALUES));
    PetscCall(MatSetValue(C,i,i,tau*3.0,INSERT_VALUES));
    if (i<n-1) PetscCall(MatSetValue(C,i,i+1,-tau,INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY));

  /* M is a diagonal matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD,&M));
  PetscCall(MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(M));
  PetscCall(MatGetOwnershipRange(M,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) PetscCall(MatSetValue(M,i,i,muu,INSERT_VALUES));
  PetscCall(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY));

  PetscCall(PetscOptionsGetBool(NULL,NULL,"-nonhyperbolic",&nohyp,NULL));
  A[0] = K; A[1] = C; A[2] = M;
  if (nohyp) {
    s = c*.6;
    PetscCall(TransformMatricesMoebius(A,UNKNOWN_NONZERO_PATTERN,c,s,-s,c,At));
    for (i=0;i<3;i++) PetscCall(MatDestroy(&A[i]));
    Op = At;
  } else Op = A;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  PetscCall(PEPCreate(PETSC_COMM_WORLD,&pep));
  PetscCall(PEPSetProblemType(pep,PEP_HERMITIAN));
  PetscCall(PEPSetType(pep,PEPSTOAR));
  /*
     Set operators and set problem type
  */
  PetscCall(PEPSetOperators(pep,3,Op));

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */
  PetscCall(PEPGetST(pep,&st));
  PetscCall(STGetKSP(st,&ksp));
  PetscCall(KSPSetType(ksp,KSPPREONLY));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCCHOLESKY));

  /*
     Use MUMPS if available.
     Note that in complex scalars we cannot use MUMPS for spectrum slicing,
     because MatGetInertia() is not available in that case.
  */
#if PetscDefined(HAVE_MUMPS) && !PetscDefined(USE_COMPLEX)
  PetscCall(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));
  /*
     Add several MUMPS options (see ex43.c for a better way of setting them in program):
     '-st_mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
  */
  PetscCall(PetscOptionsInsertString(NULL,"-st_mat_mumps_icntl_13 1 -st_mat_mumps_icntl_24 1 -st_mat_mumps_cntl_3 1e-12"));
#endif

  /*
     Set solver parameters at runtime
  */
  PetscCall(PEPSetFromOptions(pep));

  PetscCall(PetscOptionsGetBool(NULL,NULL,"-transform",&transform,NULL));
  if (transform) {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Check if the problem is definite
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(PEPCheckDefiniteQEP(pep,&xi,&mu,&def,&hyp));
    switch (def) {
      case 1:
        if (hyp==1) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Hyperbolic Problem xi=%g\n",(double)xi));
        else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Definite Problem xi=%g mu=%g\n",(double)xi,(double)mu));
        break;
      case -1:
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Not Definite Problem\n"));
        break;
      default:
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Cannot determine definiteness\n"));
        break;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Transform the QEP to have a definite inner product in the linearization
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if (def==1) {
      PetscCall(QEPDefiniteTransformGetMatrices(pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,B));
      PetscCall(PEPSetOperators(pep,3,B));
      PetscCall(PEPGetTarget(pep,&target));
      targett = target;
      PetscCall(QEPDefiniteTransformMap(hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,1,&targett,PETSC_FALSE));
      PetscCall(PEPSetTarget(pep,targett));
      PetscCall(PEPGetProblemType(pep,&type));
      PetscCall(PEPSetProblemType(pep,PEP_HYPERBOLIC));
      PetscCall(PEPSTOARGetLinearization(pep,&alpha,&beta));
      PetscCall(PEPSTOARSetLinearization(pep,1.0,0.0));
      PetscCall(PEPGetInterval(pep,&inta,&intb));
      if (inta!=intb) {
        ats[0] = inta; ats[1] = intb;
        PetscCall(QEPDefiniteTransformMap(hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,2,ats,PETSC_FALSE));
        at[0] = PetscRealPart(ats[0]); at[1] = PetscRealPart(ats[1]);
        if (at[0]<at[1]) PetscCall(PEPSetInterval(pep,at[0],at[1]));
        else PetscCall(PEPSetInterval(pep,PETSC_MIN_REAL,at[1]));
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPSolve(pep));

  /* show detailed info unless -terse option is given by user */
  if (def!=1) {
    PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
    if (terse) PetscCall(PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL));
    else {
      PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
      PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
      PetscCall(PEPErrorView(pep,PEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
      PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
    }
  } else {
    /* Map the solution */
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(QEPDefiniteCheckError(Op,pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu));
    for (i=0;i<3;i++) PetscCall(MatDestroy(B+i));
  }
  if (at[0]>at[1]) {
    PetscCall(PEPSetInterval(pep,at[0],PETSC_MAX_REAL));
    PetscCall(PEPSolve(pep));
    PetscCall(PEPConvergedReasonView(pep,PETSC_VIEWER_STDOUT_WORLD));
    /* Map the solution */
    PetscCall(QEPDefiniteCheckError(Op,pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu));
  }
  if (def==1) {
    PetscCall(PEPSetTarget(pep,target));
    PetscCall(PEPSetProblemType(pep,type));
    PetscCall(PEPSTOARSetLinearization(pep,alpha,beta));
    if (inta!=intb) PetscCall(PEPSetInterval(pep,inta,intb));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PEPDestroy(&pep));
  for (i=0;i<3;i++) PetscCall(MatDestroy(Op+i));
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformMap_Initial - map a scalar value with a certain Moebius transform

                   a theta + b
         lambda = --------------
                   c theta + d

  Input:
    xi,mu: real values such that Q(xi)<0 and Q(mu)>0
    hyperbolic: if true the problem is assumed hyperbolic (mu is not used)
  Input/Output:
    val (array of length n)
    if backtransform=true returns lambda from theta, else returns theta from lambda
*/
static PetscErrorCode QEPDefiniteTransformMap_Initial(PetscBool hyperbolic,PetscReal xi,PetscReal mu,PetscInt n,PetscScalar *val,PetscBool backtransform)
{
  PetscInt  i;
  PetscReal a,b,c,d,s;

  PetscFunctionBegin;
  if (hyperbolic) { a = 1.0; b = xi; c =0.0; d = 1.0; }
  else { a = mu; b = mu*xi-1; c = 1.0; d = xi+mu; }
  if (!backtransform) { s = a; a = -d; d = -s; }
  for (i=0;i<n;i++) {
    if (PetscRealPart(val[i]) >= PETSC_MAX_REAL || PetscRealPart(val[i]) <= PETSC_MIN_REAL) val[i] = a/c;
    else if (val[i] == -d/c) val[i] = PETSC_MAX_REAL;
    else val[i] = (a*val[i]+b)/(c*val[i]+d);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformMap - perform the mapping if the problem is hyperbolic, otherwise
  modify the value of xi in advance
*/
PetscErrorCode QEPDefiniteTransformMap(PetscBool hyperbolic,PetscReal xi,PetscReal mu,PetscInt n,PetscScalar *val,PetscBool backtransform)
{
  PetscReal      xit;
  PetscScalar    alpha;

  PetscFunctionBegin;
  xit = xi;
  if (!hyperbolic) {
    alpha = xi;
    PetscCall(QEPDefiniteTransformMap_Initial(PETSC_FALSE,0.0,mu,1,&alpha,PETSC_FALSE));
    xit = PetscRealPart(alpha);
  }
  PetscCall(QEPDefiniteTransformMap_Initial(hyperbolic,xit,mu,n,val,backtransform));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
  TransformMatricesMoebius - transform the coefficient matrices of a QEP

  Input:
    A: coefficient matrices of the original QEP
    a,b,c,d: parameters of the Moebius transform
    str: structure flag for MatAXPY operations
  Output:
    B: transformed matrices
*/
PetscErrorCode TransformMatricesMoebius(Mat A[3],MatStructure str,PetscReal a,PetscReal b,PetscReal c,PetscReal d,Mat B[3])
{
  PetscInt       i,k;
  PetscReal      cf[9];

  PetscFunctionBegin;
  for (i=0;i<3;i++) PetscCall(MatDuplicate(A[2],MAT_COPY_VALUES,&B[i]));
  /* Ct = b*b*A+b*d*B+d*d*C */
  cf[0] = d*d; cf[1] = b*d; cf[2] = b*b;
  /* Bt = 2*a*b*A+(b*c+a*d)*B+2*c*d*C*/
  cf[3] = 2*c*d; cf[4] = b*c+a*d; cf[5] = 2*a*b;
  /* At = a*a*A+a*c*B+c*c*C */
  cf[6] = c*c; cf[7] = a*c; cf[8] = a*a;
  for (k=0;k<3;k++) {
    PetscCall(MatScale(B[k],cf[k*3+2]));
    for (i=0;i<2;i++) PetscCall(MatAXPY(B[k],cf[3*k+i],A[i],str));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformGetMatrices - given a PEP of degree 2, transform the three
  matrices with TransformMatricesMoebius

  Input:
    pep: polynomial eigenproblem to be transformed, with Q(.) being the quadratic polynomial
    xi,mu: real values such that Q(xi)<0 and Q(mu)>0
    hyperbolic: if true the problem is assumed hyperbolic (mu is not used)
  Output:
    T: coefficient matrices of the transformed polynomial
*/
PetscErrorCode QEPDefiniteTransformGetMatrices(PEP pep,PetscBool hyperbolic,PetscReal xi,PetscReal mu,Mat T[3])
{
  MatStructure   str;
  ST             st;
  PetscInt       i;
  PetscReal      a,b,c,d;
  PetscScalar    xit;
  Mat            A[3];

  PetscFunctionBegin;
  for (i=2;i>=0;i--) PetscCall(PEPGetOperators(pep,i,&A[i]));
  if (hyperbolic) { a = 1.0; b = xi; c =0.0; d = 1.0; }
  else {
    xit = xi;
    PetscCall(QEPDefiniteTransformMap_Initial(PETSC_FALSE,0.0,mu,1,&xit,PETSC_FALSE));
    a = mu; b = mu*PetscRealPart(xit)-1.0; c = 1.0; d = PetscRealPart(xit)+mu;
  }
  PetscCall(PEPGetST(pep,&st));
  PetscCall(STGetMatStructure(st,&str));
  PetscCall(TransformMatricesMoebius(A,str,a,b,c,d,T));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
  Auxiliary function to compute the residual norm of an eigenpair of a QEP defined
  by coefficient matrices A
*/
static PetscErrorCode PEPResidualNorm(Mat *A,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  PetscInt       i,nmat=3;
  PetscScalar    vals[3];
  Vec            u,w;
#if !PetscDefined(USE_COMPLEX)
  Vec            ui,wi;
  PetscReal      ni;
  PetscBool      imag;
  PetscScalar    ivals[3];
#endif

  PetscFunctionBegin;
  u = z[0]; w = z[1];
  PetscCall(VecSet(u,0.0));
#if !PetscDefined(USE_COMPLEX)
  ui = z[2]; wi = z[3];
#endif
  vals[0] = 1.0;
  vals[1] = kr;
  vals[2] = kr*kr-ki*ki;
#if !PetscDefined(USE_COMPLEX)
  ivals[0] = 0.0;
  ivals[1] = ki;
  ivals[2] = 2.0*kr*ki;
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON))
    imag = PETSC_FALSE;
  else {
    imag = PETSC_TRUE;
    PetscCall(VecSet(ui,0.0));
  }
#endif
  for (i=0;i<nmat;i++) {
    if (vals[i]!=0.0) {
      PetscCall(MatMult(A[i],xr,w));
      PetscCall(VecAXPY(u,vals[i],w));
    }
#if !PetscDefined(USE_COMPLEX)
    if (imag) {
      if (ivals[i]!=0 || vals[i]!=0) {
        PetscCall(MatMult(A[i],xi,wi));
        if (vals[i]==0) PetscCall(MatMult(A[i],xr,w));
      }
      if (ivals[i]!=0) {
        PetscCall(VecAXPY(u,-ivals[i],wi));
        PetscCall(VecAXPY(ui,ivals[i],w));
      }
      if (vals[i]!=0) PetscCall(VecAXPY(ui,vals[i],wi));
    }
#endif
  }
  PetscCall(VecNorm(u,NORM_2,norm));
#if !PetscDefined(USE_COMPLEX)
  if (imag) {
    PetscCall(VecNorm(ui,NORM_2,&ni));
    *norm = SlepcAbsEigenvalue(*norm,ni);
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteCheckError - check and print the residual norm of a transformed PEP

  Input:
    A: coefficient matrices of the original problem
    pep: solver containing the computed solution of the transformed problem
    xi,mu,hyperbolic: parameters used in transformation
*/
PetscErrorCode QEPDefiniteCheckError(Mat *A,PEP pep,PetscBool hyperbolic,PetscReal xi,PetscReal mu)
{
  PetscScalar    er,ei;
  PetscReal      re,im,error;
  Vec            vr,vi,w[4];
  PetscInt       i,nconv;
  BV             bv;
  char           ex[30],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  PetscCall(PetscSNPrintf(ex,sizeof(ex),"||P(k)x||/||kx||"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%s            k             %s\n%s",sep,ex,sep));
  PetscCall(PEPGetConverged(pep,&nconv));
  PetscCall(PEPGetBV(pep,&bv));
  PetscCall(BVCreateVec(bv,w));
  PetscCall(VecDuplicate(w[0],&vr));
  PetscCall(VecDuplicate(w[0],&vi));
  for (i=1;i<4;i++) PetscCall(VecDuplicate(w[0],w+i));
  for (i=0;i<nconv;i++) {
    PetscCall(PEPGetEigenpair(pep,i,&er,&ei,vr,vi));
    PetscCall(QEPDefiniteTransformMap(hyperbolic,xi,mu,1,&er,PETSC_TRUE));
    PetscCall(PEPResidualNorm(A,er,0.0,vr,vi,w,&error));
    error /= SlepcAbsEigenvalue(er,0.0);
#if PetscDefined(USE_COMPLEX)
    re = PetscRealPart(er);
    im = PetscImaginaryPart(ei);
#else
    re = er;
    im = ei;
#endif
    if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  % 9f%+9fi      %12g\n",(double)re,(double)im,(double)error));
    else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"    % 12f           %12g\n",(double)re,(double)error));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%s",sep));
  for (i=0;i<4;i++) PetscCall(VecDestroy(w+i));
  PetscCall(VecDestroy(&vi));
  PetscCall(VecDestroy(&vr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      requires: !single
      args: -pep_nev 3 -nonhyperbolic -pep_target 2
      output_file: output/ex40_1.out
      filter: grep -v "Definite" | sed -e "s/iterations [0-9]\([0-9]*\)/iterations xx/g" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      test:
         suffix: 1
         requires: !complex
      test:
         suffix: 1_complex
         requires: complex !mumps
      test:
         suffix: 1_transform
         requires: !complex
         args: -transform
      test:
         suffix: 1_transform_complex
         requires: complex !mumps
         args: -transform

TEST*/
