!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcsvddef
        use slepcsysdef
#include <../ftn/svd/slepcall.h>
        end module slepcsvddef

        module slepcsvd
        use slepcsvddef
        use slepceps
#include <../ftn/svd/slepcall.h90>
#include <../src/svd/ftn-mod/slepcsvd.h90>

        contains

#include <../ftn/svd/slepcall.hf90>
        end module slepcsvd
