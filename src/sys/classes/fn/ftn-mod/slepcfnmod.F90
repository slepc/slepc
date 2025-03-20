!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcfndef
        use slepcsysdef
#include <../ftn/sys/classes/fn/slepcall.h>
        end module slepcfndef

        module slepcfn
        use slepcfndef
        use slepcsys
#include <../ftn/sys/classes/fn/slepcall.h90>
#include <../src/sys/classes/fn/ftn-mod/slepcfn.h90>

        contains

#include <../ftn/sys/classes/fn/slepcall.hf90>
        end module slepcfn
