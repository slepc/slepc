!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcbvdef
        use slepcsysdef
#include <../ftn/sys/classes/bv/slepcall.h>
        end module slepcbvdef

        module slepcbv
        use slepcbvdef
        use slepcsys
#include <../ftn/sys/classes/bv/slepcall.h90>

        contains

#include <../ftn/sys/classes/bv/slepcall.hf90>
        end module slepcbv
