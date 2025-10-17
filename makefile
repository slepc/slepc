#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# This is the top level makefile for compiling SLEPc.
#   * make help - useful messages on functionality
#   * make all  - compile the SLEPc libraries and utilities
#   * make check - runs a quick test that the libraries are built correctly and SLEPc applications can run
#
#   * make install - for use with ./configure is run with the --prefix=directory option
#   * make test - runs a comprehensive test suite (requires gnumake)
#   * make alldoc - build the entire SLEPc documentation (locally)
#   * a variety of rules that print library properties useful for building applications (use make help)
#   * a variety of rules for SLEPc developers
#
# gmakefile - manages the compiling SLEPc in parallel
# gmakefile.test - manages running the comprehensive test suite
#
# This makefile does not require GNUmake
ALL: all

# Include the rest of makefiles
include ./${PETSC_ARCH}/lib/slepc/conf/slepcvariables
include ${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/slepcvariables  # required in prefix builds
include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules_util.mk

# This makefile doesn't really do any work. Sub-makes still benefit from parallelism.
.NOTPARALLEL:

OMAKE_SELF = $(OMAKE) -f makefile
OMAKE_SELF_PRINTDIR = $(OMAKE_PRINTDIR) -f makefile

# ******** Rules for make all **************************************************************************

.PHONY: all
all:
	+@${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} chk_slepcdir | tee ${PETSC_ARCH}/lib/slepc/conf/make.log
	@ln -sf ${PETSC_ARCH}/lib/slepc/conf/make.log make.log
	+@(${OMAKE_SELF_PRINTDIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} all-local; echo "$$?" > ${PETSC_ARCH}/lib/slepc/conf/error.log) 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;
	+@if [ "`cat ${PETSC_ARCH}/lib/slepc/conf/error.log 2> /dev/null`" != "0" ]; then \
	    grep -E '(out of memory allocating.*after a total of|gfortran: fatal error: Killed signal terminated program f951|f95: fatal error: Killed signal terminated program f951)' ${PETSC_ARCH}/lib/slepc/conf/make.log | tee ${PETSC_ARCH}/lib/slepc/conf/memoryerror.log > /dev/null; \
	    if test -s ${PETSC_ARCH}/lib/slepc/conf/memoryerror.log; then \
              printf ${PETSC_TEXT_HILIGHT}"**************************ERROR*************************************\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
              echo "  Error during compile, you need to increase the memory allocated to the VM and rerun " 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
              printf "********************************************************************"${PETSC_TEXT_NORMAL}"\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;\
            else \
              printf ${PETSC_TEXT_HILIGHT}"*******************************ERROR************************************\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
              echo "  Error during compile, check ${PETSC_ARCH}/lib/slepc/conf/make.log" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
              echo "  Send all contents of ./${PETSC_ARCH}/lib/slepc/conf to slepc-maint@upv.es" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;\
              printf "************************************************************************"${PETSC_TEXT_NORMAL}"\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
            fi; \
          elif [ "${SLEPC_INSTALLDIR}" = "${SLEPC_DIR:/=}/${PETSC_ARCH}" ]; then \
           echo "=========================================";\
           echo "Now to check if the libraries are working do:";\
           echo "${MAKE_USER} SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} check";\
           echo "=========================================";\
         else \
           echo "=========================================";\
           echo "Now to install the library do:";\
           echo "${MAKE_USER} SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} install";\
           echo "=========================================";\
         fi
	@echo "Finishing make run at `date +'%a, %d %b %Y %H:%M:%S %z'`" >> ${PETSC_ARCH}/lib/slepc/conf/make.log
	@if [ "`cat ${PETSC_ARCH}/lib/slepc/conf/error.log 2> /dev/null`" != "0" ]; then exit 1; fi

.PHONY: all-local
all-local: info slepc_libs ${SLEPC_POST_BUILDS}

${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/files:
	@touch -t 197102020000 ${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/files

${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles:
	@${MKDIR} -p ${SLEPC_DIR}/${PETSC_ARCH}/tests && touch -t 197102020000 ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles

slepc_libs: ${SLEPC_DIR}/${PETSC_ARCH}/lib/slepc/conf/files ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	+@r=`echo "${MAKEFLAGS}" | grep ' -j'`; \
        if [ "$$?" = 0 ]; then make_j=""; else make_j="-j${MAKE_NP}"; fi; \
        r=`echo "${MAKEFLAGS}" | grep ' -l'`; \
        if [ "$$?" = 0 ]; then make_l=""; else make_l="-l${MAKE_LOAD}"; fi; \
        cmd="${OMAKE_PRINTDIR} -f gmakefile $${make_j} $${make_l} ${MAKE_PAR_OUT_FLG} V=${V} slepc_libs"; \
        cd ${SLEPC_DIR} && echo $${cmd} && exec $${cmd}

.PHONY: chk_slepcdir
chk_slepcdir:
	@mypwd=`pwd`; cd ${SLEPC_DIR} 2>&1 > /dev/null; true_SLEPC_DIR=`pwd`; cd $${mypwd} 2>&1 >/dev/null; \
        newpwd=`echo $${mypwd} | sed "s+$${true_SLEPC_DIR}+DUMMY+g"`;\
        hasslepc=`echo $${mypwd} | sed "s+slepc-+DUMMY+g"`;\
        if [ $${mypwd} = $${newpwd} -a $${hasslepc} != $${mypwd} ]; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************Warning*************************\n" ; \
          echo "Your SLEPC_DIR may not match the directory you are in";\
          echo "SLEPC_DIR " $${true_SLEPC_DIR} "Current directory" $${mypwd};\
          printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ; \
        fi

.PHONY: fortranbindings
fortranbindings: deletefortranbindings
	@${PYTHON} config/utils/generatefortranbindings.py --slepc-dir=${SLEPC_DIR} --petsc-dir=${PETSC_DIR} --petsc-arch=${PETSC_ARCH}

.PHONY: deletefortranbindings
deletefortranbindings:
	-@find src -type d -name ftn-auto* | xargs rm -rf
	-@if [ -n "${PETSC_ARCH}" ] && [ -d ${PETSC_ARCH} ] && [ -d ${PETSC_ARCH}/src ]; then \
          find ${PETSC_ARCH}/src -type d -name ftn-auto* | xargs rm -rf ;\
        fi

.PHONY: reconfigure
reconfigure: allclean
	@unset MAKEFLAGS && ${PYTHON} ${PETSC_ARCH}/lib/slepc/conf/reconfigure-${PETSC_ARCH}.py

# ******** Rules for make check ************************************************************************

RUN_TEST = ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR}

.PHONY: check
check: check_body ${SLEPC_POST_CHECKS}

.PHONY: check_install
check_install: check

.PHONY: check_body
check_body:
	-@echo "Running SLEPc check examples to verify correct installation"
	-@echo "Using SLEPC_DIR=${SLEPC_DIR}, PETSC_DIR=${PETSC_DIR}, and PETSC_ARCH=${PETSC_ARCH}"
	@if [ "${PETSC_WITH_BATCH}" != "" ]; then \
           echo "Running with batch filesystem, cannot run make check"; \
        elif [ "${MPIEXEC}" = "/bin/false" ]; then \
           echo "*mpiexec not found*. cannot run make check"; \
        else \
          ${RM} -f check_error; \
          ${RUN_TEST} OMP_NUM_THREADS=1 PETSC_OPTIONS="${EXTRA_OPTIONS} ${PETSC_TEST_OPTIONS}" PATH="${PETSC_DIR}/${PETSC_ARCH}/lib:${SLEPC_DIR}/${PETSC_ARCH}/lib:${PATH}" check_build 2>&1 | tee ./${PETSC_ARCH}/lib/slepc/conf/check.log; \
          if [ -f check_error ]; then \
            echo "Error while running make check"; \
            ${RM} -f check_error; \
            exit 1; \
          fi; \
          ${RM} -f check_error; \
        fi

.PHONY: check_build
check_build:
	+@cd src/eps/tests >/dev/null; ${RUN_TEST} clean-legacy
	+@cd src/eps/tests >/dev/null; ${RUN_TEST} testtest10
	+@if [ ! "${MPI_IS_MPIUNI}" ]; then cd src/eps/tests >/dev/null; ${RUN_TEST} testtest10_mpi; fi
	+@if [ "`grep -E '^#define PETSC_USE_FORTRAN_BINDINGS 1' ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null`" = "#define PETSC_USE_FORTRAN_BINDINGS 1" ] || [ "`grep -E '^#define PETSC_USE_FORTRAN_BINDINGS 1' ${PETSC_DIR}/include/petscconf.h 2>/dev/null`" = "#define PETSC_USE_FORTRAN_BINDINGS 1" ]; then \
           cd src/eps/tests >/dev/null; ${RUN_TEST} testtest7f; \
         fi
	+@if [ "`grep -E '^#define PETSC_HAVE_CUDA 1' ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null`" = "#define PETSC_HAVE_CUDA 1" ] || [ "`grep -E '^#define PETSC_HAVE_CUDA 1' ${PETSC_DIR}/include/petscconf.h 2>/dev/null`" = "#define PETSC_HAVE_CUDA 1" ]; then \
           cd src/eps/tests >/dev/null; ${RUN_TEST} testtest10_cuda; \
         fi
	+@if [ "`grep -E '^#define PETSC_HAVE_HIP 1' ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h 2>/dev/null`" = "#define PETSC_HAVE_HIP 1" ] || [ "`grep -E '^#define PETSC_HAVE_HIP 1' ${PETSC_DIR}/include/petscconf.h 2>/dev/null`" = "#define PETSC_HAVE_HIP 1" ]; then \
           cd src/eps/tests >/dev/null; ${RUN_TEST} testtest10_hip; \
         fi
	+@if [ "`grep -E '^#define SLEPC_HAVE_BLOPEX 1' ${SLEPC_DIR}/${PETSC_ARCH}/include/slepcconf.h`" = "#define SLEPC_HAVE_BLOPEX 1" ]; then \
           cd src/eps/tests >/dev/null; ${RUN_TEST} testtest5_blopex; \
         fi
	+@cd src/eps/tests >/dev/null; ${RUN_TEST} clean-legacy
	-@echo "Completed SLEPc check examples"

# ******** Rules for make install **********************************************************************

.PHONY: install
install:
	@${PYTHON} ./config/install.py ${SLEPC_DIR} ${PETSC_DIR} ${SLEPC_INSTALLDIR} -destDir=${DESTDIR} ${PETSC_ARCH} ${AR_LIB_SUFFIX} ${RANLIB}
	+${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} SLEPC_INSTALL=$@ install-builtafterslepc

# A smaller install with fewer extras
.PHONY: install-lib
install-lib:
	@${PYTHON} ./config/install.py ${SLEPC_DIR} ${PETSC_DIR} ${SLEPC_INSTALLDIR} -destDir=${DESTDIR} -no-examples ${PETSC_ARCH} ${AR_LIB_SUFFIX} ${RANLIB}
	+${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} SLEPC_INSTALL=$@ install-builtafterslepc

.PHONY: install-builtafterslepc
install-builtafterslepc:
	@if [ "${SLEPC_POST_INSTALLS}" != "" ]; then ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} SLEPC_INSTALL=${PETSC_INSTALL} ${SLEPC_POST_INSTALLS}; fi
	@echo "*** Install of SLEPc (and any other packages) complete ***"

# ******** Rules for running the full test suite *******************************************************

.PHONY: chk_in_slepcdir
chk_in_slepcdir:
	@if [ ! -f include/slepcversion.h ]; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR **********************************************\n" ; \
          echo " This target should be invoked in top level SLEPc source dir!"; \
          printf "****************************************************************************"${PETSC_TEXT_NORMAL}"\n" ;  false; fi

TESTMODE = testexamples
ALLTESTS_CHECK_FAILURES = no
ALLTESTS_MAKEFILE = ${SLEPC_DIR}/gmakefile.test
VALGRIND=0
.PHONY: alltests
alltests: chk_in_slepcdir ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	-@${RM} -rf ${PETSC_ARCH}/lib/slepc/conf/alltests.log alltests.log
	+@if [ -f ${SLEPC_DIR}/share/slepc/examples/gmakefile.test ] ; then \
            ALLTESTS_MAKEFILE=${SLEPC_DIR}/share/slepc/examples/gmakefile.test ; \
            ALLTESTSLOG=alltests.log ;\
          else \
            ALLTESTS_MAKEFILE=${SLEPC_DIR}/gmakefile.test; \
            ALLTESTSLOG=${PETSC_ARCH}/lib/slepc/conf/alltests.log ;\
            ln -s $${ALLTESTSLOG} alltests.log ;\
          fi; \
          ${OMAKE} allgtest ALLTESTS_MAKEFILE=$${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} MPIEXEC="${MPIEXEC}" DATAFILESPATH=${DATAFILESPATH} VALGRIND=${VALGRIND} 2>&1 | tee $${ALLTESTSLOG};\
          if [ x${ALLTESTS_CHECK_FAILURES} = xyes ]; then \
            cat $${ALLTESTSLOG} | grep -E '(^not ok|not remade because of errors|^# No tests run)' | wc -l | grep '^[ ]*0$$' > /dev/null; \
          fi;

.PHONY: allgtests-tap
allgtests-tap: allgtest-tap
	+@${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} check-test-errors

.PHONY: allgtest-tap
allgtest-tap: ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	+@MAKEFLAGS="-j$(MAKE_TEST_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)" ${OMAKE} ${MAKE_SHUFFLE_FLG} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test OUTPUT=1

.PHONY: allgtest
allgtest: ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	+@MAKEFLAGS="-j$(MAKE_TEST_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)" ${OMAKE} -k -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test V=0 2>&1 | grep -E -v '^(ok [^#]*(# SKIP|# TODO|$$)|[A-Za-z][A-Za-z0-9_]*\.(c|F|cxx|F90).$$)'

.PHONY: test
test:
	+${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test
.PHONY: cleantest
cleantest:
	+${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} cleantest

# ******** Rules for cleaning **************************************************************************

.PHONY: deletelibs
deletelibs:
	-${RM} -r ${SLEPC_LIB_DIR}/libslepc*.*

.PHONY: deleteshared
deleteshared:
	@for LIBNAME in ${SHLIBS}; \
        do \
           if [ -d ${SLEPC_INSTALLDIR}/lib/$${LIBNAME}$${LIB_NAME_SUFFIX}.dylib.dSYM ]; then \
             echo ${RM} -rf ${SLEPC_INSTALLDIR}/lib/$${LIBNAME}$${LIB_NAME_SUFFIX}.dylib.dSYM; \
             ${RM} -rf ${SLEPC_INSTALLDIR}/lib/$${LIBNAME}$${LIB_NAME_SUFFIX}.dylib.dSYM; \
           fi; \
           echo ${RM} ${SLEPC_INSTALLDIR}/lib/$${LIBNAME}$${LIB_NAME_SUFFIX}.${SL_LINKER_SUFFIX}; \
           ${RM} ${SLEPC_INSTALLDIR}/lib/$${LIBNAME}$${LIB_NAME_SUFFIX}.${SL_LINKER_SUFFIX}; \
        done
	@if [ -f ${SLEPC_INSTALLDIR}/lib/so_locations ]; then \
          echo ${RM} ${SLEPC_INSTALLDIR}/lib/so_locations; \
          ${RM} ${SLEPC_INSTALLDIR}/lib/so_locations; \
        fi

.PHONY: deletemods
deletemods:
	-${RM} -f ${SLEPC_DIR}/${PETSC_ARCH}/include/slepc*.mod

.PHONY: allclean
allclean:
	-@${OMAKE} -f gmakefile clean

.PHONY: clean
clean:: allclean

#********* Rules for printing library properties useful for building applications **********************

.PHONY: getversion_slepc
getversion_slepc:
	-@${SLEPC_DIR}/lib/slepc/bin/slepcversion

.PHONY: getlinklibs_slepc
getlinklibs_slepc:
	-@${OMAKE} -f gmakefile gmakegetlinklibs_slepc

.PHONY: getincludedirs_slepc
getincludedirs_slepc:
	-@${OMAKE} -f gmakefile gmakegetincludedirs_slepc

.PHONY: info
info:
	+@${OMAKE} -f gmakefile gmakeinfo

.PHONY: check_usermakefile
check_usermakefile:
	-@echo "Testing compile with user makefile"
	-@echo "Using SLEPC_DIR=${SLEPC_DIR}, PETSC_DIR=${PETSC_DIR}, and PETSC_ARCH=${PETSC_ARCH}"
	@cd src/eps/tutorials; ${RUN_TEST} clean-legacy
	@cd src/eps/tutorials; ${OMAKE} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} -f ${SLEPC_DIR}/share/slepc/Makefile.user ex10
	@if [ -f ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h ]; then \
           grep -E "^#define PETSC_USE_FORTRAN_BINDINGS 1" ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h | tee .ftn.log > /dev/null; \
         elif [ -f ${PETSC_DIR}/include/petscconf.h ]; then \
           grep -E "^#define PETSC_USE_FORTRAN_BINDINGS 1" ${PETSC_DIR}/include/petscconf.h | tee .ftn.log > /dev/null; \
         fi; \
         if test -s .ftn.log; then \
          cd src/eps/tutorials; ${OMAKE} SLEPC_DIR=${SLEPC_DIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} -f ${SLEPC_DIR}/share/slepc/Makefile.user ex10f; \
         fi; ${RM} .ftn.log;
	@cd src/eps/tutorials; ${RUN_TEST} clean-legacy
	-@echo "Completed compile with user makefile"

# ******** Rules for generating tag files **************************************************************

.PHONY: alletags
alletags:
	-@${PYTHON} lib/slepc/bin/maint/generateetags.py
	-@find config -type f -name "*.py" |grep -v SCCS | xargs etags -o TAGS_PYTHON

# ******** Rules for building documentation ************************************************************

.PHONY: alldoc
alldoc: doc_html

.PHONY: chk_c2html
chk_c2html:
	@if [ ${C2HTML}foo = foo ] ; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
          echo "Require c2html for html docs. Please reconfigure PETSc with --download-c2html=1"; \
          printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; fi

.PHONY: chk_doctext
chk_doctext:
	@if [ ${DOCTEXT}foo = foo ] ; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
          echo "Require sowing for html docs. Please reconfigure PETSc with --download-sowing=1"; \
          printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; fi

# Build just PDF manual + prerequisites
.PHONY: doc_pdf
doc_pdf:
	${OMAKE_SELF} -C doc latexpdf PETSC_DIR=${PETSC_DIR}

# Builds .html versions of the source
.PHONY: doc_html
doc_html: chk_c2html chk_doctext
	${OMAKE_SELF} -C doc website PETSC_DIR=${PETSC_DIR}

# Builds only .html version of the source
.PHONY: doc_html_only
doc_html_only: chk_c2html chk_doctext
	${OMAKE_SELF} -C doc html_only PETSC_DIR=${PETSC_DIR}

# Deletes documentation
.PHONY: alldocclean
alldocclean:
	-@${OMAKE_SELF} -C doc clean PETSC_DIR=${PETSC_DIR}

# ******** Rules for checking coding standards *********************************************************

# Run fortitude Fortran linter
.PHONY: fortitude
fortitude:
	-@fortitude check --line-length 1000 --ignore C003,C121,S241 --verbose --fix --preview

# Compare ABI/API of two versions of PETSc library with the old one defined by PETSC_{DIR,ARCH}_ABI_OLD
.PHONY: abitest
abitest:
	@if [ "x${SLEPC_DIR_ABI_OLD}" = "x" ] || [ "x${PETSC_ARCH_ABI_OLD}" = "x" ] || [ "x${PETSC_DIR_ABI_OLD}" = "x" ]; \
         then printf "You must set environment variables SLEPC_DIR_ABI_OLD, PETSC_ARCH_ABI_OLD, and PETSC_DIR_ABI_OLD to run abitest\n"; \
           exit 1; \
        fi;
	-@echo "Comparing ABI/API of the following two SLEPc versions (you must have already configured and built them using GCC and with -g):"
	-@echo "========================================================================================="
	-@echo "    Old: SLEPC_DIR_ABI_OLD  = ${SLEPC_DIR_ABI_OLD}"
	-@echo "         PETSC_ARCH_ABI_OLD = ${PETSC_ARCH_ABI_OLD}"
	-@echo "         PETSC_DIR_ABI_OLD  = ${PETSC_DIR_ABI_OLD}"
	-@cd ${SLEPC_DIR_ABI_OLD}; echo "         Branch             = "`git rev-parse --abbrev-ref HEAD`
	-@echo "    New: SLEPC_DIR          = ${SLEPC_DIR}"
	-@echo "         PETSC_ARCH         = ${PETSC_ARCH}"
	-@echo "         PETSC_DIR          = ${PETSC_DIR}"
	-@echo "         Branch             = "`git rev-parse --abbrev-ref HEAD`
	-@echo "========================================================================================="
	-@$(PYTHON) ${SLEPC_DIR}/lib/slepc/bin/maint/abicheck.py -old_dir ${SLEPC_DIR_ABI_OLD} -old_arch ${PETSC_ARCH_ABI_OLD} -old_petsc_dir ${PETSC_DIR_ABI_OLD} -new_dir ${SLEPC_DIR} -new_arch ${PETSC_ARCH} -new_petsc_dir ${PETSC_DIR} -report_format html

