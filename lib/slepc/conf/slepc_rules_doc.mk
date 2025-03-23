#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

# The following additional variables are used by PETSc documentation targets
#
# SOURCEALL  - sources and includes
# SOURCED    - sources/includes [but not Fortran - for doc parsing]
#
# Note that EXAMPLESALL is only used in the tutorial directories and SOURCED only in the non-tutorials and tests directories
#
SOURCEALL   = `ls *.c *.cxx *.F *.F90 *.cu *.cpp *.h *.hpp 2> /dev/null`
SOURCED     = `ls *.c *.cxx           *.cu *.cpp *.h *.hpp 2> /dev/null`
EXAMPLESALL = `(ls *.c *.cxx *.F *.F90 *.cu *.cpp | sort -V) 2> /dev/null`

# Performs the specified action on all source/include directories except output; used by c2html and cleanhtml
tree: ${ACTION}
	-@for dir in `ls -d */ 2> /dev/null` foo ;  do \
            if [[ $${dir} != "doc/" && $${dir} != "output/" ]]; then \
              if [[ -f $${dir}makefile ]]; then \
	        (cd $$dir ; ${OMAKE} ACTION=${ACTION} PETSC_ARCH=${PETSC_ARCH}  LOC=${LOC} tree) ; \
              fi; \
           fi; \
	 done
#
# Performs the specified action on all directories
slepc_tree: ${ACTION}
	-@for dir in `ls -d */ 2> /dev/null` foo ; do \
            if [[ -f $${dir}makefile ]]; then \
              (cd $$dir ; ${OMAKE} ACTION=${ACTION} PETSC_ARCH=${PETSC_ARCH} LOC=${LOC} slepc_tree); \
            fi; \
          done

# Performs the specified action on all source/include directories except tutorial and test
slepc_tree_src: ${ACTION}
	-@for dir in `ls -d */ 2> /dev/null` foo ;  do \
            if [[ $${dir} != "tests/" && $${dir} != "tutorials/" && $${dir} != "doc/" && $${dir} != "output/" ]]; then \
              if [[ -f $${dir}makefile ]]; then \
                (cd $$dir ; ${OMAKE} ACTION=${ACTION} PETSC_ARCH=${PETSC_ARCH}  LOC=${LOC} slepc_tree_src) ; \
              fi; \
           fi; \
         done

slepc_manualpages:
	-@slepc_dir=$$(realpath ${SLEPC_DIR}); LOCDIR=$$(pwd | sed s"?$${slepc_dir}/??g")/; \
        if [ "${MANSEC}" = "" ] ; then \
          for f in ${SOURCED}; do \
            LMANSEC=`grep SUBMANSEC $${f} | sed s'?[ ]*/\*[ ]*SUBMANSEC[ ]*=[ ]*\([a-zA-Z]*\)[ ]*\*/?\1?'g`; \
            if [ "$${LMANSEC}" = "" ] ; then LMANSEC="MissingSUBMANSEC"; fi; \
            DOCTEXT_PATH=${PETSC_DIR}/doc/manualpages/doctext; export DOCTEXT_PATH; \
            ${DOCTEXT} -html \
                   -mpath ${LOC}/docs/manualpages/$${LMANSEC} -heading SLEPc \
                   -defn ${SLEPC_DIR}/src/docs/doctext/html.def \
                   -locdir $${LOCDIR} -mapref ${LOC}/docs/manualpages/htmlmap -Wargdesc $${f} 2>&1 | tee -a ${SLEPC_DIR}/${PETSC_ARCH}/manualpages.err; \
            if [ -f "${LOC}/docs/manualpages/$${LMANSEC}" ]; then chmod g+w "${LOC}"/docs/manualpages/$${LMANSEC}/*; fi; \
          done; \
        else \
          if [ "${SUBMANSEC}" = "" ] ; then LMANSEC=${MANSEC}; else LMANSEC=${SUBMANSEC}; fi; \
          DOCTEXT_PATH=${PETSC_DIR}/doc/manualpages/doctext; export DOCTEXT_PATH; \
          ${DOCTEXT} -html \
                 -mpath ${LOC}/docs/manualpages/$${LMANSEC} -heading SLEPc \
                 -defn ${SLEPC_DIR}/src/docs/doctext/html.def \
                 -locdir $${LOCDIR} -mapref ${LOC}/docs/manualpages/htmlmap -Wargdesc ${SOURCED} 2>&1 | tee -a ${SLEPC_DIR}/${PETSC_ARCH}/manualpages.err; \
          if [ -f "${LOC}/docs/manualpages/$${LMANSEC}" ]; then chmod g+w "${LOC}"/docs/manualpages/$${LMANSEC}/*; fi; \
        fi

slepc_manexamples:
	-@base=$$(basename $$(pwd)); \
        if [ "$${base}" = "tutorials" -o "$${base}" = "nlevp" -o "$${base}" = "cnetwork" ] ; then \
          slepc_dir=$$(realpath ${SLEPC_DIR}); LOCDIR=$$(pwd | sed s"?$${slepc_dir}/??g")/; \
          echo "Generating manual example links in $${LOCDIR}" ; \
          for i in ${EXAMPLESALL} foo ; do \
            if [ "$$i" != "foo" ] ; then \
              a=`cat $$i | ${MAPNAMES} -map ${LOC}/docs/manualpages/manualpages.cit \
                   -printmatch-link -o /dev/null| cut -f 2 | cut -d '#' -f 1 |sed -e s~^../~~  | sort | uniq` ;  \
              for j in $$a ; do \
                b=`ls ${LOC}/docs/manualpages/$${j} | grep -v /all/ | cut -f9` ; \
                l=`grep -e "^<A HREF=\"\.\./\.\./\.\..*/tutorials/" -e "^<A HREF=\"\.\./\.\./\.\..*/nlevp/" -e "^<A HREF=\"\.\./\.\./\.\..*/cnetwork/" $${b} | wc -l`; \
                if [ $$l -le 10 ] ; then \
                  if [ $$l -eq 0 ] ; then \
                    echo "<P><H3><FONT COLOR=\"#883300\">Examples</FONT></H3>" >> $$b; \
                  fi; \
                  echo  "<A HREF=\"../../../XX\">BB</A><BR>" | sed s?XX?$${LOCDIR}$$i.html?g | sed s?BB?$${LOCDIR}$$i?g >> $$b; \
                  grep -v /BODY $$b > ltmp; \
                  echo "</BODY></HTML>" >> ltmp; \
                  mv -f ltmp $$b; \
                fi; \
              done; \
            fi; \
          done; \
        fi

petsc_manualpages_buildcite:
	-@export petscidx_tmp=$$(mktemp -d); \
          if [ ! -d "${LOC}/docs/manualpages" ]; then ${MKDIR} ${LOC}/docs/manualpages; fi; \
          echo Generating index of PETSc man pages; \
          petscrelease=`grep '^#define PETSC_VERSION_RELEASE ' ${PETSC_DIR}/include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
          if [ $${petscrelease} = 1 ]; then petscbranch="release"; else petscbranch="main"; fi; \
          DOCTEXT_PATH=${PETSC_DIR}/doc/manualpages/doctext; export DOCTEXT_PATH; \
          TEXTFILTER_PATH=${PETSC_DIR}/doc/manualpages/doctext; export TEXTFILTER_PATH; \
          cp ${PETSC_DIR}/include/*.h $$petscidx_tmp; \
          doctext_common_def=${PETSC_DIR}/doc/manualpages/doctext/doctextcommon.txt; \
          for f in `ls $${petscidx_tmp}`; do \
            LMANSEC=`grep SUBMANSEC $${petscidx_tmp}/$${f} | sed s'?[ ]*/\*[ ]*SUBMANSEC[ ]*=[ ]*\([a-zA-Z]*\)[ ]*\*/?\1?'g`; \
            if [ "$${LMANSEC}" = "" ]; then LMANSEC="MissingSUBMANSEC"; fi; \
            if [ ! -d $${petscidx_tmp}/$${LMANSEC} ]; then ${MKDIR} $${petscidx_tmp}/$${LMANSEC}; fi; \
            ${DOCTEXT} -html -indexdir "https://petsc.org/$${petscbranch}/manualpages/$${LMANSEC}" \
              -index $${petscidx_tmp}/petscmanualpages.cit \
              -mpath $${petscidx_tmp}/$${LMANSEC} $${doctext_common_def} $${petscidx_tmp}/$${f}; \
          done; \
          if [ -d ${PETSC_DIR}/src ]; then \
            for f in `find ${PETSC_DIR}/src \( -name tutorials -o -name tests -o -name "ftn-*" \) -prune -o -name "*.c" -print`; do \
              base=`basename $$f`; \
              makef=$${f%$${base}}makefile; \
              if [ -f $${makef} ]; then \
                LMANSEC=`grep SUBMANSEC $${makef} | grep -v BFORTSUBMANSEC | cut -d= -f2 | tr -d " \t"`; \
                if [ "$${LMANSEC}" = "" ]; then LMANSEC=`grep MANSEC $${makef} | cut -d= -f2 | tr -d " \t" | head -n 1`; fi; \
                if [ "$${LMANSEC}" != "" ]; then \
                  if [ ! -d $${petscidx_tmp}/$${LMANSEC} ]; then ${MKDIR} $${petscidx_tmp}/$${LMANSEC}; fi; \
                  ${DOCTEXT} -html -indexdir "https://petsc.org/$${petscbranch}/manualpages/$${LMANSEC}" \
                    -index $${petscidx_tmp}/petscmanualpages.cit \
                    -mpath $${petscidx_tmp}/$${LMANSEC} $${doctext_common_def} $${f}; \
                fi; \
              fi; \
            done; \
          fi; \
          cp $${petscidx_tmp}/petscmanualpages.cit ${LOC}/docs/manualpages; \
          ${RM} -r $${petscidx_tmp}

slepc_manualpages_buildcite:
	@-DOCTEXT_PATH=${PETSC_DIR}/doc/manualpages/doctext; export DOCTEXT_PATH; \
          TEXTFILTER_PATH=${PETSC_DIR}/doc/manualpages/doctext; export TEXTFILTER_PATH; \
          if [ "${MANSEC}" = "" ] ; then \
            for f in ${SOURCED}; do \
              LMANSEC=`grep SUBMANSEC $${f} | sed s'?[ ]*/\*[ ]*SUBMANSEC[ ]*=[ ]*\([a-zA-Z]*\)[ ]*\*/?\1?'g`; \
              if [ "$${LMANSEC}" = "" ] ; then \
                LMANSEC="MissingSUBMANSEC"; \
                pwd | grep -e ftn-custom -e slepc/finclude -e slepc/private > /dev/null; fnd=$?; \
                if [ "$${fnd}" == "1" ]; then \
                  echo "Missing MANSEC or SUBMANSEC definition in " `pwd`/$${f} ; \
                fi; \
              fi; \
              if [ ! -d "${LOC}/docs/manualpages/$${LMANSEC}" ] && ([ "$${LMANSEC}" != "MissingSUBMANSEC" ] || [ "$${fnd}" == "1" ]); then \
                echo Making directory ${LOC}/docs/manualpages/$${LMANSEC} for manual pages from $${f}; ${MKDIR} ${LOC}/docs/manualpages/$${LMANSEC}; \
              fi; \
              ${DOCTEXT} -html -indexdir ../$${LMANSEC} \
                -index ${LOC}/docs/manualpages/manualpages.cit \
                -mpath ${LOC}/docs/manualpages/$${LMANSEC} $${f}; \
            done; \
          else \
            if [ "${SUBMANSEC}" = "" ] ; then LMANSEC=${MANSEC}; else LMANSEC=${SUBMANSEC}; fi; \
            if [ ! -d "${LOC}/docs/manualpages/$${LMANSEC}" ]; then \
              echo Making directory ${LOC}/docs/manualpages/$${LMANSEC} for manual pages; ${MKDIR} ${LOC}/docs/manualpages/$${LMANSEC}; \
            fi; \
            ${DOCTEXT} -html -indexdir ../$${LMANSEC} \
              -index ${LOC}/docs/manualpages/manualpages.cit \
              -mpath ${LOC}/docs/manualpages/$${LMANSEC} ${SOURCED}; \
          fi

slepc_html:
	-@export htmlmap_tmp=$$(mktemp) ;\
          slepc_dir=$$(realpath ${SLEPC_DIR}); LOCDIR=$$(pwd | sed s"?$${slepc_dir}/??g")/; \
          sed -e s?man+../?man+ROOT/docs/manualpages/? ${LOC}/docs/manualpages/manualpages.cit > $$htmlmap_tmp ;\
          cat ${LOC}/docs/manualpages/petscmanualpages.cit >> $$htmlmap_tmp ;\
          cat ${PETSC_DIR}/doc/manualpages/mpi.www.index >> $$htmlmap_tmp ;\
          ROOT=`echo $${LOCDIR} | sed -e s?/[a-z0-9-]*?/..?g -e s?src/??g -e s?include/??g` ;\
          loc=`pwd | sed -e s?\$${SLEPC_DIR}?$${LOC}/?g -e s?/disks??g`;  \
          ${MKDIR} -p $${loc} ;\
          for i in ${SOURCEALL} ${EXAMPLESALL} foo ; do\
            if [ -f $$i ]; then \
              idir=`dirname $$i`;\
              if [ ! -d $${loc}/$${idir} ]; then ${MKDIR} -p $${loc}/$${idir}; fi ; \
              iroot=`echo $$i | sed -e "s?[a-z.]*/??g"`;\
              IROOT=`echo $${i} | sed -e s?[.][.]??g` ;\
              if [ "$${IROOT}" != "$${i}" ] ; then \
                IROOT=".."; \
              else \
                IROOT=$${ROOT};\
              fi;\
              ${RM} $${loc}/$$i.html; \
              echo "<center><a href=\"$${iroot}\">Actual source code: $${iroot}</a></center><br>" > $${loc}/$$i.html; \
              sed -E "s/PETSC[A-Z]*_DLLEXPORT//g" $$i | \
              ${C2HTML} -n | \
              awk '{ sub(/<pre width="80">/,"<pre width=\"80\">\n"); print }' | \
              ${PYTHON} ${SLEPC_DIR}/lib/slepc/bin/maint/fixinclude.py $$i $${SLEPC_DIR} | \
              grep -E -v '(PetscValid|#if !defined\(__|#define __|#undef __|EXTERN_C )' | \
              ${MAPNAMES} -map $$htmlmap_tmp -inhtml | sed -e s?ROOT?$${IROOT}?g >> $${loc}/$$i.html ; \
            fi; \
          done ;\
          loc=`pwd | sed -e s?\$${SLEPC_DIR}?$${LOC}/?g -e s?/disks??g`; ${RM} $${loc}/index.html; \
          if [ -f ${SLEPC_DIR}/src/docs/manualpages-sec/${MANSEC} ] ; then \
            cat ${SLEPC_DIR}/src/docs/manualpages-sec/${MANSEC} | sed -e "s?<A HREF=\"SLEPC_DIR[a-z/]*\">Examples</A>?<A HREF=\"$${ROOT}/docs/manualpages/${MANSEC}\">Manual pages</A>?g" -e "s?SLEPC_DIR?$${ROOT}/?g"> $${loc}/index.html; \
          else \
            touch $${loc}/index.html; \
          fi; \
          echo "<p>" >> $${loc}/index.html ;\
          loc=`pwd | sed -e s?\$${SLEPC_DIR}?$${LOC}/?g -e s?/disks??g`;\
          base=`basename $${LOCDIR}`; \
          if [ "$${base}" = "tutorials" -o "$${base}" = "nlevp" -o "$${base}" = "cnetwork" ] ; then \
            for file in ${EXAMPLESALL} foo ; do \
              if [ -f $$file ]; then \
                if [ "$${file%.F}.F" = "$${file}" -o "$${file%.F90}.F90" = "$${file}" ]; then \
                  cmess=`grep "Description:" $${file} | cut -d: -f2`; \
                else \
                  cmess=`grep "static\( const\)\? char help" $${file} | cut -d\" -f2 | cut -d\. -f1`; \
                fi; \
                echo "<a href=\"$${file}.html\">$${file}: $${cmess}</a><br>" >> $${loc}/index.html;\
              fi; \
            done ;\
          else \
            for file in `ls -d */ 2> /dev/null` foo; do \
              if [ -d $$file ]; then \
                echo "<a href=\"$${file}/\">$${file}/</a><br>" >> $${loc}/index.html; \
              fi; \
            done; \
            echo " " >> $${loc}/index.html; \
            for file in ${SOURCEALL} foo ; do \
              if [ -f $$file ]; then \
                echo "<a href=\"$${file}.html\">$${file}</a><br>" >> $${loc}/index.html; \
              fi; \
            done; \
            echo " " >> $${loc}/index.html; \
          fi ;\
          ${RM} $$htmlmap_tmp

cleanhtml:
	-@${RM} index.html *.{c,cxx,cu,F,F90,h,h90,m}.html *.{c,cxx,cu}.gcov.html
