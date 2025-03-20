!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcdsdef
        use slepcsysdef
#include <../ftn/sys/classes/ds/slepcall.h>
        end module slepcdsdef

        module slepcds
        use slepcdsdef
        use slepcfn
        use slepcrg
#include <../ftn/sys/classes/ds/slepcall.h90>
#include <../src/sys/classes/ds/ftn-mod/slepcds.h90>

        contains

#include <../ftn/sys/classes/ds/slepcall.hf90>
        end module slepcds
