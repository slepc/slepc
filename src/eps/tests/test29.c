/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the computation of left eigenvectors for generalized eigenproblems.\n\n"
  "The command line options are:\n"
  "  -f1 <filename> -f2 <filename>, PETSc binary files containing A and B\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode ComputeResidualNorm(Mat,Mat,PetscBool,PetscScalar,PetscScalar,Vec,Vec,Vec*,PetscReal*);

int main(int argc,char **argv)
{
  Mat            A,B;
  EPS            eps;
  EPSType        type;
  PetscInt       i,nconv;
  PetscBool      twosided,flg;
  PetscReal      nrmr,nrml=0.0,re,im,lev;
  PetscScalar    *kr,*ki;
  Vec            t,*xr,*xi,*yr,*yi,*z;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem stored in file.\n\n"));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-f1",filename,sizeof(filename),&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Must indicate a file name for matrix A with the -f1 option");

#if defined(PETSC_USE_COMPLEX)
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...\n"));
#else
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...\n"));
#endif
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatLoad(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(PetscOptionsGetString(NULL,NULL,"-f2",filename,sizeof(filename),&flg));
  if (flg) {
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&B));
    PetscCall(MatSetFromOptions(B));
    PetscCall(MatLoad(B,viewer));
    PetscCall(PetscViewerDestroy(&viewer));
  } else {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Matrix B was not provided, setting B=I\n\n"));
    B = NULL;
  }
  PetscCall(MatCreateVecs(A,NULL,&t));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,B));

  /* use a two-sided algorithm to compute left eigenvectors as well */
  PetscCall(EPSSetTwoSided(eps,PETSC_TRUE));

  /* allow user to change settings at run time */
  PetscCall(EPSSetFromOptions(eps));
  PetscCall(EPSGetTwoSided(eps,&twosided));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));

  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get number of converged approximate eigenpairs
  */
  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));
  PetscCall(PetscMalloc2(nconv,&kr,nconv,&ki));
  PetscCall(VecDuplicateVecs(t,3,&z));
  PetscCall(VecDuplicateVecs(t,nconv,&xr));
  PetscCall(VecDuplicateVecs(t,nconv,&xi));
  if (twosided) {
    PetscCall(VecDuplicateVecs(t,nconv,&yr));
    PetscCall(VecDuplicateVecs(t,nconv,&yi));
  }

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k            ||Ax-kBx||         ||y'A-y'Bk||\n"
         "   ---------------- ------------------ ------------------\n"));

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      PetscCall(EPSGetEigenpair(eps,i,&kr[i],&ki[i],xr[i],xi[i]));
      if (twosided) PetscCall(EPSGetLeftEigenvector(eps,i,yr[i],yi[i]));
      /*
         Compute the residual norms associated to each eigenpair
      */
      PetscCall(ComputeResidualNorm(A,B,PETSC_FALSE,kr[i],ki[i],xr[i],xi[i],z,&nrmr));
      if (twosided) PetscCall(ComputeResidualNorm(A,B,PETSC_TRUE,kr[i],ki[i],yr[i],yi[i],z,&nrml));

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr[i]);
      im = PetscImaginaryPart(kr[i]);
#else
      re = kr[i];
      im = ki[i];
#endif
      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %8f%+8fi %12g       %12g\n",(double)re,(double)im,(double)nrmr,(double)nrml));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12g\n",(double)re,(double)nrmr,(double)nrml));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
    /*
       Check bi-orthogonality of eigenvectors
    */
    if (twosided) {
      PetscCall(VecCheckOrthogonality(xr,nconv,yr,nconv,B,NULL,&lev));
      if (lev<100*PETSC_MACHINE_EPSILON) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors < 100*eps\n\n"));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors: %g\n\n",(double)lev));
    }
  }

  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(VecDestroy(&t));
  PetscCall(PetscFree2(kr,ki));
  PetscCall(VecDestroyVecs(3,&z));
  PetscCall(VecDestroyVecs(nconv,&xr));
  PetscCall(VecDestroyVecs(nconv,&xi));
  if (twosided) {
    PetscCall(VecDestroyVecs(nconv,&yr));
    PetscCall(VecDestroyVecs(nconv,&yi));
  }
  PetscCall(SlepcFinalize());
  return 0;
}

/*
   ComputeResidualNorm - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     trans - whether A' must be used instead of A
     kr,ki - eigenvalue
     xr,xi - eigenvector
     z     - three work vectors (the second one not referenced in complex scalars)
*/
PetscErrorCode ComputeResidualNorm(Mat A,Mat B,PetscBool trans,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  Vec            u,w=NULL;
  PetscScalar    alpha;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v;
  PetscReal      ni,nr;
#endif
  PetscErrorCode (*matmult)(Mat,Vec,Vec) = trans? MatMultHermitianTranspose: MatMult;

  PetscFunctionBegin;
  u = z[0];
  if (B) w = z[2];

#if !defined(PETSC_USE_COMPLEX)
  v = z[1];
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON)) {
#endif
    PetscCall((*matmult)(A,xr,u));                          /* u=A*x */
    if (PetscAbsScalar(kr) > PETSC_MACHINE_EPSILON) {
      if (B) PetscCall((*matmult)(B,xr,w));             /* w=B*x */
      else w = xr;
      alpha = trans? -PetscConj(kr): -kr;
      PetscCall(VecAXPY(u,alpha,w));                        /* u=A*x-k*B*x */
    }
    PetscCall(VecNorm(u,NORM_2,norm));
#if !defined(PETSC_USE_COMPLEX)
  } else {
    PetscCall((*matmult)(A,xr,u));                          /* u=A*xr */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      if (B) PetscCall((*matmult)(B,xr,v));             /* v=B*xr */
      else PetscCall(VecCopy(xr,v));
      PetscCall(VecAXPY(u,-kr,v));                          /* u=A*xr-kr*B*xr */
      if (B) PetscCall((*matmult)(B,xi,w));             /* w=B*xi */
      else w = xi;
      PetscCall(VecAXPY(u,trans?-ki:ki,w));                 /* u=A*xr-kr*B*xr+ki*B*xi */
    }
    PetscCall(VecNorm(u,NORM_2,&nr));
    PetscCall((*matmult)(A,xi,u));                          /* u=A*xi */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      PetscCall(VecAXPY(u,-kr,w));                          /* u=A*xi-kr*B*xi */
      PetscCall(VecAXPY(u,trans?ki:-ki,v));                 /* u=A*xi-kr*B*xi-ki*B*xr */
    }
    PetscCall(VecNorm(u,NORM_2,&ni));
    *norm = SlepcAbsEigenvalue(nr,ni);
  }
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

   testset:
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_nev 4 -st_type sinvert -eps_target -190000
      filter: grep -v "method" | sed -e "s/[+-]0\.0*i//g" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 1
      test:
         suffix: 1_rqi
         args: -eps_type power -eps_power_shift_type rayleigh -eps_nev 2 -eps_target -2000
      test:
         suffix: 1_rqi_singular
         args: -eps_type power -eps_power_shift_type rayleigh -eps_nev 1 -eps_target -195500

   test:
      suffix: 2
      args: -f1 ${DATAFILESPATH}/matrices/complex/mhd1280a.petsc -f2 ${DATAFILESPATH}/matrices/complex/mhd1280b.petsc -eps_nev 6 -eps_tol 1e-11
      filter: sed -e "s/-892/+892/" | sed -e "s/-759/+759/" | sed -e "s/-674/+674/" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)
      timeoutfactor: 2

   test:
      suffix: 3
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_nev 4 -st_type sinvert -eps_target -250000
      filter: grep -v "method" | sed -e "s/[+-]0\.0*i//g" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g" | grep -v bi-orthog
      requires: double !complex !defined(PETSC_USE_64BIT_INDICES)

TEST*/
