!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcnepdef
        use slepcsysdef
        use slepcpepdef
#include <../ftn/nep/slepcall.h>
        end module slepcnepdef

        module slepcnep
        use slepcnepdef
        use slepcpep
#include <../ftn/nep/slepcall.h90>
#include <../src/nep/ftn-mod/slepcnep.h90>

        contains

#include <../ftn/nep/slepcall.hf90>
        end module slepcnep

! The following module imports all the functionality of SLEPc and PETSc
        module slepc
        use slepcnep
        use slepcmfn
        use petsc
        end module
