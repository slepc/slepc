/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Use the matrix exponential to compute rightmost eigenvalues.\n\n"
  "Same problem as ex9.c but with explicitly created matrix. The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>
#include <slepcmfn.h>

/*
   This example computes the eigenvalues with largest real part of the
   following matrix

        A = [ tau1*T+(beta-1)*I     alpha^2*I
                  -beta*I        tau2*T-alpha^2*I ],

   where

        T = tridiag{1,-2,1}
        h = 1/(n+1)
        tau1 = delta1/(h*L)^2
        tau2 = delta2/(h*L)^2

   but it builds A explicitly, as opposed to ex9.c
*/

/* Routines for shell spectral transformation */
PetscErrorCode STApply_Exp(ST,Vec,Vec);
PetscErrorCode STBackTransform_Exp(ST,PetscInt,PetscScalar*,PetscScalar*);
PetscErrorCode STDestroy_Exp(void**);

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation context */
  MFN            mfn;             /* matrix function solver object to compute exp(A)*v */
  FN             f;
  EPSType        type;
  PetscScalar    alpha,beta,tau1,tau2,delta1,delta2,L,h;
  PetscInt       n=30,i,Istart,Iend,nev;
  PetscBool      isShell,terse;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
#if defined(PETSC_HAVE_COMPLEX)
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model with matrix exponential, n=%" PetscInt_FMT "\n\n",n));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  alpha  = 2.0;
  beta   = 5.45;
  delta1 = 0.008;
  delta2 = 0.004;
  L      = 0.51302;

  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-alpha",&alpha,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-beta",&beta,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL));
  PetscCall(PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL));

  h = 1.0 / (PetscReal)(n+1);
  tau1 = delta1 / ((h*L)*(h*L));
  tau2 = delta2 / ((h*L)*(h*L));

  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,2*n,2*n));
  PetscCall(MatSetFromOptions(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {
    if (i<n) {  /* upper blocks */
      if (i>0) PetscCall(MatSetValue(A,i,i-1,tau1,INSERT_VALUES));
      if (i<n-1) PetscCall(MatSetValue(A,i,i+1,tau1,INSERT_VALUES));
      PetscCall(MatSetValue(A,i,i,-2.0*tau1+beta-1.0,INSERT_VALUES));
      PetscCall(MatSetValue(A,i,i+n,alpha*alpha,INSERT_VALUES));
    } else {  /* lower blocks */
      if (i>n) PetscCall(MatSetValue(A,i,i-1,tau2,INSERT_VALUES));
      if (i<2*n-1) PetscCall(MatSetValue(A,i,i+1,tau2,INSERT_VALUES));
      PetscCall(MatSetValue(A,i,i,-2.0*tau2-alpha*alpha,INSERT_VALUES));
      PetscCall(MatSetValue(A,i,i-n,-beta,INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  PetscCall(EPSGetST(eps,&st));
  PetscCall(STSetType(st,STSHELL));
  PetscCall(EPSSetFromOptions(eps));

  /*
     Initialize shell spectral transformation
  */
  PetscCall(PetscObjectTypeCompare((PetscObject)st,STSHELL,&isShell));
  if (isShell) {

    /* Create the MFN object to be used by the spectral transform */
    PetscCall(MFNCreate(PETSC_COMM_WORLD,&mfn));
    PetscCall(MFNSetOperator(mfn,A));
    PetscCall(MFNGetFN(mfn,&f));
    PetscCall(FNSetType(f,FNEXP));
    PetscCall(FNSetScale(f,0.03,1.0));  /* this can be set with -fn_scale */
    PetscCall(MFNSetFromOptions(mfn));

    /* Set callback functions */
    PetscCall(STShellSetApply(st,STApply_Exp));
    PetscCall(STShellSetBackTransform(st,STBackTransform_Exp));
    PetscCall(STShellSetContext(st,mfn));
    PetscCall(STShellSetContextDestroy(st,STDestroy_Exp));
    PetscCall(PetscObjectSetName((PetscObject)st,"STEXP"));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL));
  else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
#else
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This example requires C99 complex numbers");
#endif
  PetscCall(SlepcFinalize());
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   STBackTransform_Exp - Undoes the exp(A) transformation by taking logarithms.

   Input Parameters:
+  st - spectral transformation context
-  n  - number of eigenvalues to transform

   Input/Output Parameters:
+  eigr - pointer to real part of eigenvalues
-  eigi - pointer to imaginary part of eigenvalues
*/
PetscErrorCode STBackTransform_Exp(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
#if defined(PETSC_HAVE_COMPLEX)
  PetscInt       j;
  MFN            mfn;
  FN             fn;
  PetscScalar    tau,eta;
#if !defined(PETSC_USE_COMPLEX)
  PetscComplex   theta,lambda;
#endif

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&mfn));
  PetscCall(MFNGetFN(mfn,&fn));
  PetscCall(FNGetScale(fn,&tau,&eta));
  for (j=0;j<n;j++) {
#if defined(PETSC_USE_COMPLEX)
    eigr[j] = PetscLogComplex(eigr[j]/eta)/tau;
#else
    theta   = PetscCMPLX(eigr[j],eigi[j])/eta;
    lambda  = PetscLogComplex(theta)/tau;
    eigr[j] = PetscRealPartComplex(lambda);
    eigi[j] = PetscImaginaryPartComplex(lambda);
#endif
  }
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}

/*
   STApply_Exp - Applies the operator exp(tau*A) to a given vector using an MFN object.

   Input Parameters:
+  st - spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApply_Exp(ST st,Vec x,Vec y)
{
  MFN            mfn;

  PetscFunctionBeginUser;
  PetscCall(STShellGetContext(st,&mfn));
  PetscCall(MFNSolve(mfn,x,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   STDestroy_Exp - This routine destroys the shell ST context.

   Input Parameter:
.  ctx - user-defined spectral transformation context
*/
PetscErrorCode STDestroy_Exp(void **ctx)
{
  MFN       mfn = (MFN)*ctx;

  PetscFunctionBeginUser;
  PetscCall(MFNDestroy(&mfn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -eps_nev 4 -mfn_ncv 16 -terse
      requires: c99_complex !single
      filter: sed -e "s/-2/+2/g"
      output_file: output/ex36_1.out
      test:
         suffix: 1
         requires: !__float128
      test:
         suffix: 1_quad
         args: -eps_tol 1e-14
         requires: __float128

   test:
      suffix: 2
      args: -n 56 -eps_nev 2 -st_type sinvert -eps_target -390 -eps_target_magnitude -eps_type power
      args: -eps_power_shift_type {{constant rayleigh}} -eps_two_sided {{0 1}} -eps_tol 1e-14 -terse
      requires: c99_complex !single
      filter: sed -e "s/[+-]0\.0*i//g"

   test:
      suffix: 3
      args: -n 100 -st_type sinvert -eps_type ciss -rg_type ellipse -rg_ellipse_center 0 -rg_ellipse_radius 6 -eps_all -eps_tol 1e-6 -terse
      requires: c99_complex !single
      filter: sed -e "s/-3.37036-3.55528i, -3.37036+3.55528i/-3.37036+3.55528i, -3.37036-3.55528i/" -e "s/-1.79853-3.03216i, -1.79853+3.03216i/-1.79853+3.03216i, -1.79853-3.03216i/" -e "s/-0.67471-2.52856i, -0.67471+2.52856i/-0.67471+2.52856i, -0.67471-2.52856i/" -e "s/0.00002-2.13950i, 0.00002+2.13950i/0.00002+2.13950i, 0.00002-2.13950i/"

TEST*/
