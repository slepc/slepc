!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcpepdef
        use slepcsysdef
#include <../ftn/pep/slepcall.h>
        end module slepcpepdef

        module slepcpep
        use slepcpepdef
        use slepceps
#include <../ftn/pep/slepcall.h90>
#include <../src/pep/ftn-mod/slepcpep.h90>

        contains

#include <../ftn/pep/slepcall.hf90>
        end module slepcpep
