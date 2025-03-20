!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcrgdef
        use slepcsysdef
#include <../ftn/sys/classes/rg/slepcall.h>
        end module slepcrgdef

        module slepcrg
        use slepcrgdef
        use slepcsys
#include <../ftn/sys/classes/rg/slepcall.h90>
#include <../src/sys/classes/rg/ftn-mod/slepcrg.h90>

        contains

#include <../ftn/sys/classes/rg/slepcall.hf90>
        end module slepcrg
