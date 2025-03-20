!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcmfndef
        use slepcsysdef
#include <../ftn/mfn/slepcall.h>
        end module slepcmfndef

        module slepcmfn
        use slepcmfndef
        use slepcbv
        use slepcfn
#include <../ftn/mfn/slepcall.h90>

        contains

#include <../ftn/mfn/slepcall.hf90>
        end module slepcmfn
