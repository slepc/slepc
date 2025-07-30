# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     := source
PETSC_ARCH    ?= arch-doc
BUILDDIR      = ../${PETSC_ARCH}/doc

# Put it first so that "make" without argument is like "make help".
help:
	-@echo "Preparing python virtual environment (.venv-sphinx)"
	-@python3 -m venv .venv-sphinx ;\
	. .venv-sphinx/bin/activate ;\
	python3 -m pip install --quiet -r requirements.txt ;\
	$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) ;\
	deactivate \
	-@echo "--"
	-@echo "SLEPc targets:"
	-@echo "  sphinxhtml  to build the HTML website"

.PHONY: help Makefile

clean:
	-@${RM} -rf ${BUILDDIR}/*
	-@${RM} -rf source/_static

webstatic:
	git clone https://gitlab.com/slepc/webstatic.git -b main --depth=1 source/_static

# use this rule to build the documentation on your machine when working on SLEPc documentation
website: latexpdf
	-@python3 -m venv .venv-sphinx ;\
	. .venv-sphinx/bin/activate ;\
	python3 -m pip install --quiet -r requirements.txt ;\
	${RM} -rf ${BUILDDIR}/* ;\
	make SPHINXOPTS="-T -E -W -j 4" html ;\
	deactivate

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

latexpdf: webstatic logo_upv
	-@echo "Preparing python virtual environment (.venv-sphinx)"
	-@python3 -m venv .venv-sphinx ;\
	. .venv-sphinx/bin/activate ;\
	python3 -m pip install -r requirements.txt ;\
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) ;\
	if [ ! -d $(SOURCEDIR)/_static/manual ]; then \
		mkdir $(SOURCEDIR)/_static/manual ;\
	fi; \
	cp $(BUILDDIR)/latex/slepc-manual.pdf $(SOURCEDIR)/_static/manual/; \
	deactivate \

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%:
	-@echo "Preparing python virtual environment (.venv-sphinx)"
	-@python3 -m venv .venv-sphinx ;\
	. .venv-sphinx/bin/activate ;\
	python3 -m pip install --quiet -r requirements.txt ;\
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) ;\
	deactivate \
