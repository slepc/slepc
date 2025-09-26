/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a GNHEP problem with contour integral.\n\n"
  "Based on ex7.\n"
  "The command line options are:\n"
  "  -f1 <filename> -f2 <filename>, PETSc binary files containing A and B.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  EPS               eps;
  RG                rg;
  Mat               A,B;
  PetscBool         flg;
  EPSCISSExtraction extr;
  EPSCISSQuadRule   quad;
  char              filename[PETSC_MAX_PATH_LEN];
  PetscViewer       viewer;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nGNHEP problem with contour integral\n\n"));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
  PetscCall(EPSSetOperators(eps,A,B));
  PetscCall(EPSSetProblemType(eps,EPS_GNHEP));
  PetscCall(EPSSetTolerances(eps,1e-9,PETSC_CURRENT));

  /* set CISS solver with various options */
  PetscCall(EPSSetType(eps,EPSCISS));
  PetscCall(EPSCISSSetExtraction(eps,EPS_CISS_EXTRACTION_HANKEL));
  PetscCall(EPSCISSSetQuadRule(eps,EPS_CISS_QUADRULE_CHEBYSHEV));
  PetscCall(EPSCISSSetUseST(eps,PETSC_TRUE));
  PetscCall(EPSGetRG(eps,&rg));
  PetscCall(RGSetType(rg,RGINTERVAL));
  PetscCall(RGIntervalSetEndpoints(rg,-3000.0,0.0,0.0,0.0));

  PetscCall(EPSSetFromOptions(eps));

  PetscCall(EPSSolve(eps));
  PetscCall(PetscObjectTypeCompare((PetscObject)eps,EPSCISS,&flg));
  if (flg) {
    PetscCall(EPSCISSGetExtraction(eps,&extr));
    PetscCall(EPSCISSGetQuadRule(eps,&quad));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solved with CISS using %s extraction with %s quadrature rule\n\n",EPSCISSExtractions[extr],EPSCISSQuadRules[quad]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSErrorView(eps,EPS_ERROR_BACKWARD,NULL));
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&B));
  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   testset:
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc
      output_file: output/test25_1.out
      test:
         suffix: 1
         requires: double !complex !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 1_cuda
         args: -mat_type aijcusparse -st_pc_factor_mat_solver_type cusparse
         requires: cuda double !complex !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 1_hip
         args: -mat_type aijhipsparse -st_pc_factor_mat_solver_type hipsparse
         requires: hip double !complex !defined(PETSC_USE_64BIT_INDICES)

   testset:
      args: -f1 ${DATAFILESPATH}/matrices/complex/qc324.petsc -rg_type ellipse -rg_ellipse_center 1-0.09i -rg_ellipse_radius 0.15 -rg_ellipse_vscale 0.1
      output_file: output/test25_2.out
      test:
         suffix: 2
         requires: double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 2_cuda
         args: -mat_type aijcusparse -st_pc_factor_mat_solver_type cusparse
         requires: cuda double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 2_hip
         args: -mat_type aijhipsparse -st_pc_factor_mat_solver_type hipsparse
         requires: hip double complex datafilespath !defined(PETSC_USE_64BIT_INDICES)

TEST*/
