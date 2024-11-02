!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcstdef
        use slepcsysdef
#include <../ftn/sys/classes/st/slepcall.h>
        end module slepcstdef

        module slepcst
        use petscksp
        use slepcstdef
        use slepcbv
#include <../ftn/sys/classes/st/slepcall.h90>

        contains

#include <../ftn/sys/classes/st/slepcall.hf90>
        end module slepcst
