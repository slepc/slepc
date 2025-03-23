!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcepsdef
        use slepcsysdef
#include <../ftn/eps/slepcall.h>
        end module slepcepsdef

        module slepceps
        use petscsnes
        use slepcepsdef
        use slepcst
        use slepcds
        use slepclme
#include <../ftn/eps/slepcall.h90>
#include <../src/eps/ftn-mod/slepceps.h90>

        contains

#include <../ftn/eps/slepcall.hf90>

        end module slepceps
