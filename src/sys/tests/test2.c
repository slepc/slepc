/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests functions intended to be used from a debugger.\n\n";

#include <slepcsys.h>

int main(int argc,char **argv)
{
#if defined(PETSC_USE_DEBUG)
#if defined(PETSC_USE_COMPLEX)
  PetscScalar Xr[]={1.0,-0.5,0.625,1.25,-0.125,-5.5};
#else
  PetscScalar Xr[]={1.0,-0.5,0.625,1.25,-0.125,-5.5},Xi[]={0.0,0.0,0.0,0.0,0.0,0.0};
#endif
#endif

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,NULL,help));

#if defined(PETSC_USE_DEBUG)
#if defined(PETSC_USE_COMPLEX)
  PetscCall(SlepcDebugViewMatrix(2,3,Xr,NULL,2,"M",NULL));
#else
  PetscCall(SlepcDebugViewMatrix(2,3,Xr,Xi,2,"M",NULL));
#endif
#endif

  PetscCall(SlepcFinalize());
  return 0;
}

/*TEST

   build:
      requires: debug

   test:
      args: -help
      filter: sed -e "s/\(Development Git.*\)/version/" | sed -e "s/\(Release Version.*\)/version/" | sed -e "s/\(linked from.*\)/linked from PATH/"

TEST*/
