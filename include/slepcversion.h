/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#ifndef SLEPCVERSION_H
#define SLEPCVERSION_H

#define SLEPC_VERSION_RELEASE    0
#define SLEPC_VERSION_MAJOR      3
#define SLEPC_VERSION_MINOR      22
#define SLEPC_VERSION_SUBMINOR   2
#define SLEPC_RELEASE_DATE       "September 29, 2024"
#define SLEPC_VERSION_DATE       "unknown"

#if !defined (SLEPC_VERSION_GIT)
#define SLEPC_VERSION_GIT        "unknown"
#endif

#if !defined(SLEPC_VERSION_DATE_GIT)
#define SLEPC_VERSION_DATE_GIT   "unknown"
#endif

#define SLEPC_VERSION_EQ(MAJOR,MINOR,SUBMINOR) \
 ((SLEPC_VERSION_MAJOR == (MAJOR)) &&       \
  (SLEPC_VERSION_MINOR == (MINOR)) &&       \
  (SLEPC_VERSION_SUBMINOR == (SUBMINOR)) && \
  (SLEPC_VERSION_RELEASE  == 1))

#define SLEPC_VERSION_ SLEPC_VERSION_EQ

#define SLEPC_VERSION_LT(MAJOR,MINOR,SUBMINOR)          \
  (SLEPC_VERSION_RELEASE == 1 &&                        \
   (SLEPC_VERSION_MAJOR < (MAJOR) ||                    \
    (SLEPC_VERSION_MAJOR == (MAJOR) &&                  \
     (SLEPC_VERSION_MINOR < (MINOR) ||                  \
      (SLEPC_VERSION_MINOR == (MINOR) &&                \
       (SLEPC_VERSION_SUBMINOR < (SUBMINOR)))))))

#define SLEPC_VERSION_LE(MAJOR,MINOR,SUBMINOR) \
  (SLEPC_VERSION_LT(MAJOR,MINOR,SUBMINOR) || \
   SLEPC_VERSION_EQ(MAJOR,MINOR,SUBMINOR))

#define SLEPC_VERSION_GT(MAJOR,MINOR,SUBMINOR) \
  (0 == SLEPC_VERSION_LE(MAJOR,MINOR,SUBMINOR))

#define SLEPC_VERSION_GE(MAJOR,MINOR,SUBMINOR) \
  (0 == SLEPC_VERSION_LT(MAJOR,MINOR,SUBMINOR))

#endif
