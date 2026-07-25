/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepcst.h>
#include <slepcds.h>
#include <slepcfn.h>
#include <slepcbv.h>
#include <slepcrg.h>

#if PetscDefined(HAVE_DYNAMIC_LIBRARIES)

#if PetscDefined(USE_SINGLE_LIBRARY)
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepceps(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcnep(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcpep(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsvd(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcmfn(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepclme(void);
#endif

/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the basic objects ST, FN, DS, BV, RG.
 */
#if PetscDefined(USE_SINGLE_LIBRARY)
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepc(void)
#else
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsys(void)
#endif
{
  PetscFunctionBegin;
  PetscCall(STInitializePackage());
  PetscCall(DSInitializePackage());
  PetscCall(FNInitializePackage());
  PetscCall(BVInitializePackage());
  PetscCall(RGInitializePackage());

#if PetscDefined(USE_SINGLE_LIBRARY)
  PetscCall(PetscDLLibraryRegister_slepceps());
  PetscCall(PetscDLLibraryRegister_slepcnep());
  PetscCall(PetscDLLibraryRegister_slepcpep());
  PetscCall(PetscDLLibraryRegister_slepcsvd());
  PetscCall(PetscDLLibraryRegister_slepcmfn());
  PetscCall(PetscDLLibraryRegister_slepclme());
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */
