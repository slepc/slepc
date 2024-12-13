!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the EPS object in SLEPc
!
#if !defined(SLEPCEPSDEF_H)
#define SLEPCEPSDEF_H

#include "slepc/finclude/slepcsys.h"
#include "slepc/finclude/slepcst.h"
#include "slepc/finclude/slepcbv.h"
#include "slepc/finclude/slepcds.h"
#include "slepc/finclude/slepcrg.h"
#include "slepc/finclude/slepclme.h"
#include "petsc/finclude/petscsnes.h"

#define EPS type(tEPS)

#define EPSType                character*(80)
#define EPSConvergedReason     PetscEnum
#define EPSErrorType           PetscEnum
#define EPSProblemType         PetscEnum
#define EPSWhich               PetscEnum
#define EPSExtraction          PetscEnum
#define EPSBalance             PetscEnum
#define EPSConv                PetscEnum
#define EPSStop                PetscEnum
#define EPSPowerShiftType      PetscEnum
#define EPSKrylovSchurBSEType  PetscEnum
#define EPSLanczosReorthogType PetscEnum
#define EPSPRIMMEMethod        PetscEnum
#define EPSCISSQuadRule        PetscEnum
#define EPSCISSExtraction      PetscEnum
#define EPSEVSLDOSMethod       PetscEnum
#define EPSEVSLDamping         PetscEnum

#define EPSPOWER       'power'
#define EPSSUBSPACE    'subspace'
#define EPSARNOLDI     'arnoldi'
#define EPSLANCZOS     'lanczos'
#define EPSKRYLOVSCHUR 'krylovschur'
#define EPSGD          'gd'
#define EPSJD          'jd'
#define EPSRQCG        'rqcg'
#define EPSLOBPCG      'lobpcg'
#define EPSCISS        'ciss'
#define EPSLYAPII      'lyapii'
#define EPSLAPACK      'lapack'
#define EPSARPACK      'arpack'
#define EPSTRLAN       'trlan'
#define EPSBLOPEX      'blopex'
#define EPSPRIMME      'primme'
#define EPSFEAST       'feast'
#define EPSSCALAPACK   'scalapack'
#define EPSELEMENTAL   'elemental'
#define EPSEVSL        'evsl'
#define EPSCHASE       'chase'

#endif
