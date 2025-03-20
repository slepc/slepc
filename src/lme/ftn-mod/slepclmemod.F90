!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepclmedef
        use slepcsysdef
#include <../ftn/lme/slepcall.h>
        end module slepclmedef

        module slepclme
        use slepclmedef
        use slepcbv
#include <../ftn/lme/slepcall.h90>

        contains

#include <../ftn/lme/slepcall.hf90>
        end module slepclme
