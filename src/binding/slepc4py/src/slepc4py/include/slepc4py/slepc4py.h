#ifndef SLEPC4PY_H
#define SLEPC4PY_H

#include <Python.h>
#include <slepc.h>

#include "../../SLEPc_api.h"

static int import_slepc4py(void) {
  if (import_slepc4py__SLEPc() < 0) goto bad;
  return 0;
 bad:
  return -1;
}

#endif /* !SLEPC4PY_H */
