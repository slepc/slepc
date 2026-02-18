#if !defined(SLEPC4PY_CUSTOM_H)
#define SLEPC4PY_CUSTOM_H

#if !defined(PETSC_ERR_PYTHON)
#define PETSC_ERR_PYTHON ((PetscErrorCode)(-1))
#endif

static PetscErrorCode SlepcInitializePackageAll(void)
{
  PetscFunctionBegin;
  PetscCall(EPSInitializePackage());
  PetscCall(SVDInitializePackage());
  PetscCall(PEPInitializePackage());
  PetscCall(NEPInitializePackage());
  PetscCall(LMEInitializePackage());
  PetscCall(MFNInitializePackage());
  PetscCall(LMEInitializePackage());
  PetscCall(STInitializePackage());
  PetscCall(BVInitializePackage());
  PetscCall(DSInitializePackage());
  PetscCall(FNInitializePackage());
  PetscCall(RGInitializePackage());
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif/*SLEPC4PY_CUSTOM_H*/
