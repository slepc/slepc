# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= -T -E -W -j 4
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     := source
PETSC_ARCH    ?= arch-doc
BUILDDIR      = ../${PETSC_ARCH}/doc

# Put it first so that "make" without argument is like "make help".
help:
	-@echo "Preparing python virtual environment (.venv-sphinx)"
	@python3 -m venv .venv-sphinx ;\
	. .venv-sphinx/bin/activate ;\
	python3 -m pip install -r requirements.txt ;\
	$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) ;\
	deactivate \
	-@echo "--"
	-@echo ""
	-@echo "SLEPc targets:"
	-@echo ""
	-@echo "  website     to build the HTML website"

.PHONY: help makefile

.PHONY: clean
clean:
	-@${RM} -rf ${BUILDDIR}/*
	-@${RM} -rf source/_static
	-@${RM} -rf source/manualpages/DS
	-@${RM} -rf source/manualpages/PEP
	-@${RM} -rf source/manualpages/LME
	-@${RM} -rf source/manualpages/FN
	-@${RM} -rf source/manualpages/SVD
	-@${RM} -rf source/manualpages/MFN
	-@${RM} -rf source/manualpages/NEP
	-@${RM} -rf source/manualpages/BV
	-@${RM} -rf source/manualpages/EPS
	-@${RM} -rf source/manualpages/Sys
	-@${RM} -rf source/manualpages/RG
	-@${RM} -rf source/manualpages/ST
	-@${RM} -rf logo-upv
	-@${RM} -rf .venv-sphinx
	-@${RM} -f source/petsc_objects.inv
	-@${MAKE} --no-print-directory -C $(SLEPC_DIR)/src/binding/slepc4py docsclean

.PHONY: webstatic
webstatic:
	@if [ ! -d source/_static ]; then \
		git clone -b main --depth=1 https://gitlab.com/slepc/webstatic.git source/_static; \
	else \
		cd source/_static; \
		git pull; \
		if [ $$? -ne 0 ]; then \
			cd - ; \
			rm -rf source/_static; \
			git clone -b main --depth=1 https://gitlab.com/slepc/webstatic.git source/_static; \
		fi \
	fi

# use this rule to build the documentation on your machine when working on SLEPc documentation
.PHONY: website
website: latexpdf
	@python3 -m venv .venv-sphinx && \
	. .venv-sphinx/bin/activate && \
	python3 -m pip install -r requirements.txt && \
	${RM} -rf ${BUILDDIR}/html && \
	$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) && \
	deactivate

.PHONY: logo_upv
logo_upv:
	@echo "Preparing UPV images";
	@if [ ! -f "source/_static/images/manual/pdf/logo-upv.pdf" ]; then \
		if [ ! -f "logo-upv/MARCA_UPV/MarcaUPV_PRINCIPAL/marca_UPV_principal_negro.svg" ]; then \
			if [ ! -f "logo-upv/archivos_MARCA_UPV.zip" ]; then \
				mkdir logo-upv; \
				wget -O logo-upv/archivos_MARCA_UPV.zip https://www.upv.es/comu/nueva-imagen/recursos/archivos_MARCA_UPV.zip ;\
				if [ $$? -ne 0 ]; then \
					printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
					echo "Could not download UPV images"; \
					printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; \
				fi; \
			fi; \
			cd logo-upv; \
			unzip archivos_MARCA_UPV.zip MARCA_UPV/MarcaUPV_PRINCIPAL/marca_UPV_principal_negro.svg; \
			if [ $$? -ne 0 ]; then \
				printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
				echo "Could not extract UPV logo from zip"; \
				printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; \
			fi; \
			cd -; \
		fi; \
		rsvg-convert -f pdf -o source/_static/images/manual/pdf/logo-upv.pdf logo-upv/MARCA_UPV/MarcaUPV_PRINCIPAL/marca_UPV_principal_negro.svg; \
		if [ $$? -ne 0 ]; then \
			printf ${PETSC_TEXT_HILIGHT}"*********************** ERROR ************************\n" ; \
			echo "Could not convert UPV logo to pdf"; \
			printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ;false; \
		fi; \
	fi

.PHONY: latexpdf
latexpdf: webstatic logo_upv
	-@echo "Preparing python virtual environment (.venv-sphinx)"
	@python3 -m venv .venv-sphinx && \
	. .venv-sphinx/bin/activate && \
	python3 -m pip install -r requirements.txt && \
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" -T -E -j 4 $(O) && \
	if [ ! -d $(SOURCEDIR)/_static/manual ]; then \
		mkdir $(SOURCEDIR)/_static/manual ;\
	fi && \
	cp $(BUILDDIR)/latex/slepc-manual.pdf $(SOURCEDIR)/_static/manual/ && \
	deactivate \

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: webstatic
	-@echo "Running Sphinx"
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
