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
#   * make docs - build the entire SLEPc documentation (locally)
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
include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules_doc.mk
include ${SLEPC_DIR}/lib/slepc/conf/slepc_rules_util.mk

# This makefile doesn't really do any work. Sub-makes still benefit from parallelism.
.NOTPARALLEL:

OMAKE_SELF = $(OMAKE) -f makefile
OMAKE_SELF_PRINTDIR = $(OMAKE_PRINTDIR) -f makefile

# ******** Rules for make all **************************************************************************

all:
	+@${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} chk_slepcdir | tee ${PETSC_ARCH}/lib/slepc/conf/make.log
	@ln -sf ${PETSC_ARCH}/lib/slepc/conf/make.log make.log
	+@${OMAKE_SELF_PRINTDIR} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} all-local 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;
	@grep -E '(out of memory allocating.*after a total of|gfortran: fatal error: Killed signal terminated program f951|f95: fatal error: Killed signal terminated program f951)' ${PETSC_ARCH}/lib/slepc/conf/make.log | tee ${PETSC_ARCH}/lib/slepc/conf/memoryerror.log > /dev/null
	@grep -v "using error" ${PETSC_ARCH}/lib/slepc/conf/make.log > /dev/null | grep -E -i "( error | error: |no such file or directory)" | tee ./${PETSC_ARCH}/lib/slepc/conf/error.log
	+@if test -s ${PETSC_ARCH}/lib/slepc/conf/memoryerror.log; then \
           printf ${PETSC_TEXT_HILIGHT}"**************************ERROR*************************************\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           echo "  Error during compile, you need to increase the memory allocated to the VM and rerun " 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           printf "********************************************************************"${PETSC_TEXT_NORMAL}"\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;\
         elif test -s ${PETSC_ARCH}/lib/slepc/conf/error.log; then \
           printf ${PETSC_TEXT_HILIGHT}"*******************************ERROR************************************\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           echo "  Error during compile, check ${PETSC_ARCH}/lib/slepc/conf/make.log" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
           echo "  Send all contents of ./${PETSC_ARCH}/lib/slepc/conf to slepc-maint@upv.es" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log;\
           printf "************************************************************************"${PETSC_TEXT_NORMAL}"\n" 2>&1 | tee -a ${PETSC_ARCH}/lib/slepc/conf/make.log; \
         elif [ "${SLEPC_INSTALLDIR}" = "${SLEPC_DIR:/=}/${PETSC_ARCH}" ]; then \
           echo "Now to check if the library is working do:";\
           echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} check";\
           echo "=========================================";\
         else \
           echo "Now to install the library do:";\
           echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} install";\
           echo "=========================================";\
         fi
	@echo "Finishing make run at `date +'%a, %d %b %Y %H:%M:%S %z'`" >> ${PETSC_ARCH}/lib/slepc/conf/make.log
	@if test -s ./${PETSC_ARCH}/lib/slepc/conf/error.log; then exit 1; fi

all-local: info slepc_libs slepc4py-build

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

allfortranstubs: deletefortranstubs
	@${PYTHON} lib/slepc/bin/maint/generatefortranstubs.py --slepc-dir=${SLEPC_DIR} --petsc-dir=${PETSC_DIR} --petsc-arch=${PETSC_ARCH} --bfort=${BFORT} --mode=generate --verbose=${V}
	-@${PYTHON} lib/slepc/bin/maint/generatefortranstubs.py --slepc-dir=${SLEPC_DIR} --petsc-dir=${PETSC_DIR} --petsc-arch=${PETSC_ARCH} --mode=merge --verbose=${V}

#copy of allfortranstubs with PETSC_ARCH=''
allfortranstubsinplace: deletefortranstubs
	@${PYTHON} lib/slepc/bin/maint/generatefortranstubs.py --slepc-dir=${SLEPC_DIR} --petsc-dir=${PETSC_DIR} --petsc-arch='' --bfort=${BFORT} --mode=generate --verbose=${V}
	-@${PYTHON} lib/slepc/bin/maint/generatefortranstubs.py --slepc-dir=${SLEPC_DIR} --petsc-dir=${PETSC_DIR} --petsc-arch='' --mode=merge --verbose=${V}

deletefortranstubs:
	-@find src -type d -name ftn-auto* | xargs rm -rf
	-@if [ -n "${PETSC_ARCH}" ] && [ -d ${PETSC_ARCH} ] && [ -d ${PETSC_ARCH}/src ]; then \
          find ${PETSC_ARCH}/src -type d -name ftn-auto* | xargs rm -rf ;\
        fi

reconfigure: allclean
	@unset MAKEFLAGS && ${PYTHON} ${PETSC_ARCH}/lib/slepc/conf/reconfigure-${PETSC_ARCH}.py

# ******** Rules for make check ************************************************************************

RUN_TEST = ${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR}

check_install: check
check:
	-@echo "Running SLEPc check examples to verify correct installation"
	-@echo "Using SLEPC_DIR=${SLEPC_DIR}, PETSC_DIR=${PETSC_DIR}, and PETSC_ARCH=${PETSC_ARCH}"
	@if [ "${PETSC_WITH_BATCH}" != "" ]; then \
           echo "Running with batch filesystem, cannot run make check"; \
        elif [ "${MPIEXEC}" = "/bin/false" ]; then \
           echo "*mpiexec not found*. cannot run make check"; \
        else \
          ${RM} -f check_error; \
          ${RUN_TEST} PETSC_OPTIONS="${PETSC_OPTIONS} ${PETSC_TEST_OPTIONS}" PATH="${PETSC_DIR}/${PETSC_ARCH}/lib:${SLEPC_DIR}/${PETSC_ARCH}/lib:${PATH}" check_build 2>&1 | tee ./${PETSC_ARCH}/lib/slepc/conf/check.log; \
          if [ -f check_error ]; then \
            echo "Error while running make check"; \
            ${RM} -f check_error; \
            exit 1; \
          fi; \
          ${RM} -f check_error; \
        fi

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

install:
	@${PYTHON} ./config/install.py ${SLEPC_DIR} ${PETSC_DIR} ${SLEPC_INSTALLDIR} -destDir=${DESTDIR} ${PETSC_ARCH} ${AR_LIB_SUFFIX} ${RANLIB}
	+${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} install-builtafterslepc

# A smaller install with fewer extras
install-lib:
	@${PYTHON} ./config/install.py ${SLEPC_DIR} ${PETSC_DIR} ${SLEPC_INSTALLDIR} -destDir=${DESTDIR} -no-examples ${PETSC_ARCH} ${AR_LIB_SUFFIX} ${RANLIB}
	+${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} install-builtafterslepc

install-builtafterslepc:
	+${OMAKE_SELF} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} slepc4py-install

# ******** Rules for running the full test suite *******************************************************

chk_in_slepcdir:
	@if [ ! -f include/slepcversion.h ]; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR **********************************************\n" ; \
          echo " This target should be invoked in top level SLEPc source dir!"; \
          printf "****************************************************************************"${PETSC_TEXT_NORMAL}"\n" ;  false; fi

TESTMODE = testexamples
ALLTESTS_CHECK_FAILURES = no
ALLTESTS_MAKEFILE = ${SLEPC_DIR}/gmakefile.test
VALGRIND=0
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

allgtests-tap: allgtest-tap
	+@${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} check-test-errors

allgtest-tap: ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	+@MAKEFLAGS="-j$(MAKE_TEST_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)" ${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test OUTPUT=1

allgtest: ${SLEPC_DIR}/${PETSC_ARCH}/tests/testfiles
	+@MAKEFLAGS="-j$(MAKE_TEST_NP) -l$(MAKE_LOAD) $(MAKEFLAGS)" ${OMAKE} -k -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test V=0 2>&1 | grep -E -v '^(ok [^#]*(# SKIP|# TODO|$$)|[A-Za-z][A-Za-z0-9_]*\.(c|F|cxx|F90).$$)'

test:
	+${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} test
cleantest:
	+${OMAKE} -f ${ALLTESTS_MAKEFILE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEPC_DIR} cleantest

# ******** Rules for cleaning **************************************************************************

deletelibs:
	-${RM} -r ${SLEPC_LIB_DIR}/libslepc*.*
deletemods:
	-${RM} -f ${SLEPC_DIR}/${PETSC_ARCH}/include/slepc*.mod

allclean:
	-@${OMAKE} -f gmakefile clean

clean:: allclean

#********* Rules for printing library properties useful for building applications **********************

getversion_slepc:
	-@${SLEPC_DIR}/lib/slepc/bin/slepcversion

getlinklibs_slepc:
	-@${OMAKE} -f gmakefile gmakegetlinklibs_slepc

getincludedirs_slepc:
	-@${OMAKE} -f gmakefile gmakegetincludedirs_slepc

info:
	+@${OMAKE} -f gmakefile gmakeinfo

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

alletags:
	-@${PYTHON} lib/slepc/bin/maint/generateetags.py
	-@find config -type f -name "*.py" |grep -v SCCS | xargs etags -o TAGS_PYTHON

# ******** Rules for building documentation ************************************************************

alldoc: allcite allpdf alldoc_pre alldoc_post docsetdate

chk_loc:
	@if [ ${LOC}foo = foo ] ; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR **********************************************\n" ; \
          echo " Please specify LOC variable for eg: make allmanpages LOC=/sandbox/slepc "; \
          printf "****************************************************************************"${PETSC_TEXT_NORMAL}"\n" ;  false; fi
	@${MKDIR} ${LOC}/manualpages

chk_c2html:
	@if [ ${C2HTML}foo = foo ] ; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
          echo "Require c2html for html docs. Please reconfigure PETSc with --download-c2html=1"; \
          printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; fi

# Build just citations
allcite: chk_loc deletemanualpages petsc_manualpages_buildcite
	-${OMAKE_SELF} ACTION=slepc_manualpages_buildcite slepc_tree_src LOC=${LOC}
	-@cat ${LOC}/docs/manualpages/manualpages.cit > ${LOC}/docs/manualpages/htmlmap
	-@cat ${LOC}/docs/manualpages/petscmanualpages.cit >> ${LOC}/docs/manualpages/htmlmap
	-@cat ${PETSC_DIR}/doc/manualpages/mpi.www.index >> ${LOC}/docs/manualpages/htmlmap

# Build just PDF manual + prerequisites
allpdf:
	-cd docs/manual; ${OMAKE_SELF} slepc.pdf clean; mv slepc.pdf ../../docs

# Build just manual pages + prerequisites
allmanpages: chk_loc allcite
	-${RM} ${SLEPC_DIR}/${PETSC_ARCH}/manualpages.err
	-${OMAKE_SELF} ACTION=slepc_manualpages slepc_tree_src LOC=${LOC}
	cat ${SLEPC_DIR}/${PETSC_ARCH}/manualpages.err
	@a=`cat ${SLEPC_DIR}/${PETSC_ARCH}/manualpages.err | wc -l`; test ! $$a -gt 0

# Build just manual examples + prerequisites
allmanexamples: chk_loc allmanpages
	-${OMAKE_SELF} ACTION=slepc_manexamples tree LOC=${LOC}

# Build everything that goes into 'doc' dir except html sources
alldoc_pre: chk_loc allcite allmanpages allmanexamples
	-${PYTHON} ${SLEPC_DIR}/lib/slepc/bin/maint/wwwindex.py ${SLEPC_DIR} ${LOC} "src/docs/manualpages-sec"
	-@echo "<html>" > singleindex.html
	-@echo "<head>" >> singleindex.html
	-@echo "  <title>Subroutine Index</title>" >> singleindex.html
	-@echo "  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" >> singleindex.html
	-@echo "  <link rel=\"stylesheet\" href=\"/slepc.css\" type=\"text/css\">" >> singleindex.html
	-@echo "</head>" >> singleindex.html
	-@echo "<body>" >> singleindex.html
	-@cat ${LOC}/docs/manualpages/singleindex.html >> singleindex.html
	-@sed -e 's/CC3333/883300/' singleindex.html > ${LOC}/docs/manualpages/singleindex.html
	-@${RM} singleindex.html

# Builds .html versions of the source
alldoc_post: chk_loc chk_c2html allcite
	-${OMAKE_SELF} ACTION=slepc_html PETSC_DIR=${PETSC_DIR} tree LOC=${LOC}
	cp ${LOC}/docs/manual.html ${LOC}/docs/index.html

# modify all generated html files and add in version number, date, canonical URL info.
docsetdate:
	@echo "Updating generated html files with slepc version, date, canonical URL info";\
        version_release=`grep '^#define SLEPC_VERSION_RELEASE ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_major=`grep '^#define SLEPC_VERSION_MAJOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_minor=`grep '^#define SLEPC_VERSION_MINOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        version_subminor=`grep '^#define SLEPC_VERSION_SUBMINOR ' include/slepcversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
        if  [ $${version_release} = 0 ]; then \
          slepcversion=slepc-main; \
          export slepcversion; \
        elif [ $${version_release} = 1 ]; then \
          slepcversion=slepc-$${version_major}.$${version_minor}.$${version_subminor}; \
          export slepcversion; \
        else \
          echo "Unknown SLEPC_VERSION_RELEASE: $${version_release}"; \
          exit; \
        fi; \
        datestr=`git log -1 --pretty=format:%ci | cut -d ' ' -f 1`; \
        export datestr; \
        gitver=`git describe --match "v*"`; \
        export gitver; \
        find include src docs/manualpages -type f -name \*.html \
          -exec perl -pi -e 's^(<body.*>)^$$1\n   <div id=\"version\" align=right><b>$$ENV{slepcversion} $$ENV{datestr}</b></div>\n   <div id="bugreport" align=right><a href="mailto:slepc-maint\@upv.es?subject=Typo or Error in Documentation &body=Please describe the typo or error in the documentation: $$ENV{slepcversion} $$ENV{gitver} {} "><small>Report Typos and Errors</small></a></div>^i' {} \; \
          -exec perl -pi -e 's^(<head>)^$$1 <link rel="canonical" href="https://slepc.upv.es/documentation/current/{}" />^i' {} \; ; \
        echo "Done fixing version number, date, canonical URL info"

# Deletes documentation
alldocclean: deletemanualpages allcleanhtml
deletemanualpages: chk_loc
	-@if [ -d ${LOC} -a -d ${LOC}/docs/manualpages ]; then \
          ${RM} -rf ${LOC}/docs/manualpages ;\
          ${RM} -f ${LOC}/docs/slepc.pdf ;\
        fi
allcleanhtml:
	-${OMAKE_SELF} ACTION=cleanhtml PETSC_DIR=${PETSC_DIR} tree

# ******** Rules for checking coding standards *********************************************************

countfortranfunctions:
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto` \
        `find ${SLEPC_DIR}/src -name ftn-custom`; do cd $$D; \
        grep -E '^void' *.c | \
        cut -d'(' -f1 | tr -s  ' ' | cut -d' ' -f3 | uniq | grep -E -v "(^$$|Petsc)" | \
        sed "s/_$$//"; done | sort > /tmp/countfortranfunctions

countcfunctions:
	-@ ls ${SLEPC_DIR}/include/*.h | grep -v slepcblaslapack.h | \
        xargs grep extern | grep "(" | tr -s ' ' | \
        cut -d'(' -f1 | cut -d' ' -f3 | grep -v "\*" | tr -s '\012' |  \
        tr 'A-Z' 'a-z' |  sort > /tmp/countcfunctions

difffortranfunctions: countfortranfunctions countcfunctions
	-@echo -------------- Functions missing in the Fortran interface ---------------------
	-@${DIFF} /tmp/countcfunctions /tmp/countfortranfunctions | grep "^<" | cut -d' ' -f2
	-@echo ----------------- Functions missing in the C interface ------------------------
	-@${DIFF} /tmp/countcfunctions /tmp/countfortranfunctions | grep "^>" | cut -d' ' -f2
	-@${RM}  /tmp/countcfunctions /tmp/countfortranfunctions

checkbadfortranstubs:
	-@echo "========================================="
	-@echo "Functions with MPI_Comm as an Argument"
	-@echo "========================================="
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
        grep '^void' *.c | grep 'MPI_Comm' | \
        tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3; done
	-@echo "========================================="
	-@echo "Functions with a String as an Argument"
	-@echo "========================================="
	-@for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
        grep '^void' *.c | grep 'char \*' | \
        tr -s ' ' | tr -s ':' ' ' |cut -d'(' -f1 | cut -d' ' -f1,3; done
	-@echo "========================================="
	-@echo "Functions with Pointers to PETSc Objects as Argument"
	-@echo "========================================="
	-@_p_OBJ=`grep _p_ ${PETSC_DIR}/include/*.h | tr -s ' ' | \
        cut -d' ' -f 3 | tr -s '\012' | grep -v '{' | cut -d'*' -f1 | \
        sed "s/_p_//g" | tr -s '\012 ' ' *|' ` ; \
        _p_OBJS=`grep _p_ ${SLEPC_DIR}/include/*.h | tr -s ' ' | \
        cut -d' ' -f 3 | tr -s '\012' | grep -v '{' | cut -d'*' -f1 | \
        sed "s/_p_//g" | tr -s '\012 ' ' *|' ` ; \
        for D in `find ${SLEPC_DIR}/src -name ftn-auto`; do cd $$D; \
        for OBJ in $$_p_OBJ $$_p_OBJS; do \
        grep "$$OBJ \*" *.c | tr -s ' ' | tr -s ':' ' ' | \
        cut -d'(' -f1 | cut -d' ' -f1,4; \
        done; done

# Compare ABI/API of two versions of PETSc library with the old one defined by PETSC_{DIR,ARCH}_ABI_OLD
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

.PHONY: info all deletelibs allclean alletags alldoc allcleanhtml countfortranfunctions install

