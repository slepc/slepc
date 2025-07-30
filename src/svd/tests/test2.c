/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test SVD with different builds with a matrix loaded from a file"
  " (matrices available in PETSc's distribution).\n\n";

#include <slepcsvd.h>

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  SVD            svd;             /* singular value problem solver context */
  char           filename[PETSC_MAX_PATH_LEN],path[PETSC_MAX_PATH_LEN];
  const char     *prefix,*scalar,*ints,*floats;
  PetscReal      tol=PETSC_SMALL;
  PetscViewer    viewer;
  PetscBool      flg;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrix for which the SVD must be computed
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined(PETSC_USE_COMPLEX)
  prefix = "nh";
  scalar = "complex";
#else
  prefix = "ns";
  scalar = "real";
#endif
#if defined(PETSC_USE_64BIT_INDICES)
  ints   = "int64";
#else
  ints   = "int32";
#endif
#if defined(PETSC_USE_REAL_DOUBLE)
  floats = "float64";
#elif defined(PETSC_USE_REAL_SINGLE)
  floats = "float32";
#endif

  PetscCall(PetscOptionsGetString(NULL,NULL,"-path",path,sizeof(path),&flg));
  if (flg) PetscCall(PetscSNPrintf(filename,sizeof(filename),"%s/%s-%s-%s-%s",path,prefix,scalar,ints,floats));
  else PetscCall(PetscSNPrintf(filename,sizeof(filename),"%s/share/petsc/datafiles/matrices/%s-%s-%s-%s",PETSC_DIR,prefix,scalar,ints,floats));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nReading matrix from binary file...\n\n"));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatLoad(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Create the SVD solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(SVDCreate(PETSC_COMM_WORLD,&svd));
  PetscCall(SVDSetOperators(svd,A,NULL));
  PetscCall(SVDSetTolerances(svd,tol,PETSC_CURRENT));
  PetscCall(SVDSetFromOptions(svd));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Compute the singular triplets and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(SVDSolve(svd));
  PetscCall(SVDErrorView(svd,SVD_ERROR_RELATIVE,NULL));
  PetscCall(SVDDestroy(&svd));
  PetscCall(MatDestroy(&A));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: !__float128

   test:
      args: -svd_nsv 7 -svd_type {{lanczos trlanczos cross cyclic lapack randomized}} -path ${DATAFILESPATH}/matrices/petsc-small
      requires: !single datafilespath

   testset:
      args: -svd_nsv 7 -svd_mpd 11 -svd_type primme -path ${DATAFILESPATH}/matrices/petsc-small
      requires: primme !single datafilespath
      output_file: output/test2_1.out
      test:
         suffix: 1_primme
      test:
         suffix: 1_primme_args
         args: -svd_primme_blocksize 2 -svd_primme_method hybrid

TEST*/
