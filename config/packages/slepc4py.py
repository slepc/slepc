#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

from __future__ import print_function
import sys, os, log, package

class Slepc4py(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename     = 'slepc4py'
    self.installable     = True
    self.ProcessArgs(argdb)

  def ProcessArgs(self,argdb,petscpackages=''):
    url,flag,found = argdb.PopUrl('download-slepc4py')
    if found:
      self.log.Exit('--download-slepc4py has been renamed to --with-slepc4py\nUse -h for help')
    value,found = argdb.PopBool('with-'+self.packagename)
    if found:
      self.requested = value
    have_petsc4py, have_petsc4py_cnt = argdb.PopBool('have-petsc4py')
    self.have_petsc4py = have_petsc4py if have_petsc4py_cnt > 0 else None

  def ShowHelp(self):
    wd = package.Package.wd
    print('  --with-slepc4py=<bool>'.ljust(wd)+': Build Python bindings (default: no)')
    print('  --have-petsc4py=<bool>'.ljust(wd)+': Whether petsc4py is installed (default: autodetect)')

  def ShowInfo(self):
    if self.havepackage:
      self.log.Println('Python bindings (slepc4py) will be built after SLEPc')

  def Process(self,slepcconf,slepcvars,slepcrules,slepc,petsc,archdir=''):
    if not self.requested:
      self.SkipInstall(slepcrules)
      return
    self.log.NewSection('Processing slepc4py...')

    pythonpath = get_python_path(petsc)
    sys.path = pythonpath.split(':') + sys.path

    # Check for pestc4py unless user suppressed this
    if self.have_petsc4py is None:
      try:
        from petsc4py import PETSc
      except ImportError:
        self.log.Exit('Cannot import petsc4py, if you have a PYTHONPATH variable make sure it contains the path where petsc4py can be found, e.g., $PETSC_DIR/$PETSC_ARCH/lib')
    elif not self.have_petsc4py:
      self.log.Exit('petsc4py is required but had been marked as not installed')

    builddir = os.path.join(slepc.dir,'src','binding','slepc4py')
    destdir  = os.path.join(slepc.prefixdir,'lib')

    # add makefile rules
    envvars = 'PYTHONPATH=%s ' % pythonpath
    if slepc.isinstall:
      envvars += 'PETSC_ARCH="" SLEPC_DIR=${DESTDIR}${SLEPC_INSTALLDIR} '
    if petsc.isIntel():
      envvars += 'CFLAGS="" '
    confdir = os.path.join(destdir,'slepc','conf')
    rule =  'slepc4pybuild:\n'
    rule += '\t@echo "*** Building slepc4py ***"\n'
    rule += '\t@${RM} -f '+os.path.join(confdir,'slepc4py.errorflg')+'\n'
    rule += '\t@cd '+builddir+' && ${RM} -rf build && \\\n'
    rule += '\t   %s ${PYTHON} setup.py build 2>&1 || \\\n' % envvars
    rule += '\t   (echo "**************************ERROR*************************************" && \\\n'
    rule += '\t   echo "Error building slepc4py." && \\\n'
    rule += '\t   echo "********************************************************************" && \\\n'
    rule += '\t   touch '+os.path.join(confdir,'slepc4py.errorflg')+' && \\\n'
    rule += '\t   exit 1)\n\n'
    slepcrules.write(rule)

    rule =  'slepc4pyinstall:\n'
    rule += '\t@echo "*** Installing slepc4py ***"\n'
    rule += '\t@cd '+builddir+' && \\\n'
    rule += '\t   %s ${PYTHON} setup.py install --install-lib=%s \\\n' % (envvars,destdir)
    rule += '\t      $(if $(DESTDIR),--root=\'$(DESTDIR)\')'
    rule += '\t   2>&1 || \\\n'
    rule += '\t   (echo "**************************ERROR*************************************" && \\\n'
    rule += '\t   echo "Error installing slepc4py" && \\\n'
    rule += '\t   echo "********************************************************************" && \\\n'
    rule += '\t   exit 1)\n'
    rule += '\t@echo "====================================="\n'
    rule += '\t@echo "To use slepc4py, add '+destdir+' to PYTHONPATH"\n'
    rule += '\t@echo "====================================="\n\n'
    slepcrules.write(rule)

    if slepc.isinstall:
      slepcrules.write('slepc4py-build:\n')
      slepcrules.write('slepc4py-install: slepc4pybuild slepc4pyinstall\n')
    else:
      slepcrules.write('slepc4py-build: slepc4pybuild slepc4pyinstall\n')
      slepcrules.write('slepc4py-install:\n')

    rule =  'slepc4pytest:\n'
    rule += '\t@echo "*** Testing slepc4py on ${PETSC4PY_NP} processes ***"\n'
    rule += '\t@PYTHONPATH=%s:%s PETSC_OPTIONS="{PETSC_OPTIONS} -check_pointer_intensity 0 -error_output_stdout -malloc_dump ${PETSC_TEST_OPTIONS}" ${MPIEXEC} -n ${PETSC4PY_NP} ${PYTHON} %s --verbose\n' % (destdir,pythonpath,os.path.join('src','binding','slepc4py','test','runtests.py'))
    rule += '\t@echo "====================================="\n\n'
    slepcrules.write(rule)

    slepcconf.write('#define SLEPC_HAVE_SLEPC4PY 1\n')
    slepcconf.write('#define SLEPC4PY_INSTALL_PATH %s\n' % destdir)
    self.havepackage = True


  def SkipInstall(self,slepcrules):
    # add empty rules
    slepcrules.write('slepc4py-build:\n')
    slepcrules.write('slepc4py-install:\n')
    slepcrules.write('slepc4pytest:\n')

def get_python_path(petsc):
  """Return the path to python packages from environment or PETSc"""
  if 'PYTHONPATH' in os.environ:
    return os.environ['PYTHONPATH']
  else:
    if petsc.isinstall:
      return os.path.join(petsc.dir,'lib')
    else:
      return os.path.join(petsc.dir,petsc.arch,'lib')
