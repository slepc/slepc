#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

include ${PETSC_DIR}/lib/petsc/conf/rules_util.mk

vermin_slepc:
	@vermin --violations -t=3.4- ${VERMIN_OPTIONS} ${SLEPC_DIR}/config
	@vermin --violations -t=3.6- ${VERMIN_OPTIONS} ${SLEPC_DIR}/src/binding/slepc4py

lint_slepc:
	${PYTHON3} ${SLEPC_DIR}/lib/slepc/bin/maint/slepcClangLinter.py $(LINTER_OPTIONS)

help-lint_slepc:
	@${PYTHON3} ${SLEPC_DIR}/lib/slepc/bin/maint/slepcClangLinter.py --help
	-@echo "Basic usage:"
	-@echo "   make lint_slepc <options>"
	-@echo
	-@echo "Options:"
	-@echo "  LINTER_OPTIONS=\"--linter_options ...\"  See above for available options"
	-@echo

