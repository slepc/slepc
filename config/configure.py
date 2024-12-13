#!/usr/bin/env python3
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
import os, sys, time, shutil

def WriteModulesFile(modules,version,sdir):
  ''' Write the contents of the Modules file '''
  modules.write('#%Module\n\n')
  modules.write('proc ModulesHelp { } {\n')
  modules.write('    puts stderr "This module sets the path and environment variables for slepc-%s"\n' % version)
  modules.write('    puts stderr "     see https://slepc.upv.es/ for more information"\n')
  modules.write('    puts stderr ""\n}\n')
  modules.write('module-whatis "SLEPc - Scalable Library for Eigenvalue Problem Computations"\n\n')
  modules.write('module load petsc\n')
  modules.write('set slepc_dir "%s"\n' % sdir)
  modules.write('setenv SLEPC_DIR "$slepc_dir"\n')

def WritePkgconfigFile(pkgconfig,version,pversion,sdir,isinstall,prefixdir,singlelib,suffix):
  ''' Write the contents of the pkg-config file '''
  pkgconfig.write('prefix=%s\n' % prefixdir)
  pkgconfig.write('exec_prefix=${prefix}\n')
  pkgconfig.write('includedir=${prefix}/include\n')
  pkgconfig.write('libdir=${prefix}/lib\n\n')
  pkgconfig.write('Name: SLEPc\n')
  pkgconfig.write('Description: the Scalable Library for Eigenvalue Problem Computations\n')
  pkgconfig.write('Version: %s\n' % version)
  pkgconfig.write('Requires: PETSc >= %s\n' % pversion)
  pkgconfig.write('Cflags: -I${includedir}')
  if not isinstall:
    pkgconfig.write(' -I'+os.path.join(sdir,'include'))
  pkgconfig.write('\nLibs:')
  if singlelib:
    pkgconfig.write(' -L${{libdir}} -lslepc{0}\n'.format(suffix))
  else:
    pkgconfig.write(' -L${{libdir}} -lslepcnep{0} -lslepcpep{0} -lslepcsvd{0} -lslepceps{0} -lslepcmfn{0} -lslepclme{0} -lslepcsys{0}\n'.format(suffix))

def WriteReconfigScript(reconfig,slepcdir,usedargs):
  ''' Write the contents of the reconfigure script '''
  reconfig.write('#!/usr/bin/env python3\n\n')
  reconfig.write('import os, sys\n')
  if usedargs:
    reconfig.write('sys.argv.extend(\''+usedargs+'\'.split())\n')
  reconfig.write('exec(open(os.path.join(\''+slepcdir+'\',\'config\',\'configure.py\')).read())\n')

def ResetConfigureHash(hashfile,log):
  ''' Removes the configure hash file '''
  try:
    log.write('Deleting configure hash file: '+hashfile)
    os.remove(hashfile)
  except:
    log.write('Unable to delete configure hash file: '+hashfile)

def Epilog(slepc,petsc):
  print()
  print('xxx'+'='*74+'xxx')
  print(' Configure stage complete. Now build the SLEPc library with:')
  if petsc.isinstall:
    print('   make SLEPC_DIR='+slepc.dir+' PETSC_DIR='+petsc.dir)
  else:
    print('   make SLEPC_DIR='+slepc.dir+' PETSC_DIR='+petsc.dir+' PETSC_ARCH='+petsc.archname)
  print('xxx'+'='*74+'xxx')
  print()

# Use en_US as language so that compiler messages are in English
def fixLang(lang):
  if lang in os.environ and os.environ[lang] != '':
    lv = os.environ[lang]
    enc = ''
    try: lv,enc = lv.split('.')
    except: pass
    if lv not in ['en_US','C']: lv = 'en_US'
    if enc: lv = lv+'.'+enc
    os.environ[lang] = lv

fixLang('LC_LOCAL')
fixLang('LANG')

# Set python path
configdir = os.path.abspath('config')
if not os.path.isdir(configdir):
  sys.exit('ERROR: Run configure from $SLEPC_DIR, not '+os.path.abspath('.'))
sys.path.insert(0,configdir)
sys.path.insert(0,os.path.join(configdir,'packages'))

# Load auxiliary classes
import argdb, log
argdb = argdb.ArgDB(sys.argv)
log   = log.Log()

showhelp = argdb.PopHelp()

# Load main classes, process the corresponding command-line options
import slepc, petsc
slepc = slepc.SLEPc(argdb,log)
petsc = petsc.PETSc(argdb,log)

# Check environment and PETSc version
if not showhelp:
  log.Print('Checking environment...')
  petsc.InitDir(slepc.prefixdir)
  petsc.LoadVersion()
slepc.InitDir()
slepc.LoadVersion()

# Load PETSc configuration
if not showhelp:
  petsc.LoadConf()
  packagesinpetsc = petsc.packages
else:
  packagesinpetsc = ''

# Load classes for packages and process their command-line options
import arpack, blopex, chase, elemental, elpa, evsl, feast, hpddm, ksvd, polar, primme, scalapack, slepc4py, slicot, trlan, sowing, lapack
arpack    = arpack.Arpack(argdb,log)
blopex    = blopex.Blopex(argdb,log)
chase     = chase.Chase(argdb,log)
elemental = elemental.Elemental(argdb,log,packagesinpetsc)
elpa      = elpa.Elpa(argdb,log)
evsl      = evsl.Evsl(argdb,log)
feast     = feast.Feast(argdb,log,packagesinpetsc)
ksvd      = ksvd.Ksvd(argdb,log)
polar     = polar.Polar(argdb,log)
primme    = primme.Primme(argdb,log)
trlan     = trlan.Trlan(argdb,log)
sowing    = sowing.Sowing(argdb,log)
lapack    = lapack.Lapack(argdb,log)
scalapack = scalapack.Scalapack(argdb,log,packagesinpetsc)
slepc4py  = slepc4py.Slepc4py(argdb,log)
slicot    = slicot.Slicot(argdb,log)
hpddm     = hpddm.HPDDM(argdb,log)

# The next list sorts the packages in a way that dependencies of X appear before X.
# SLEPc's configure does not build a graph of package dependencies, every dependency is searched linearly
externalwithdeps = [arpack, blopex, chase, elpa, evsl, hpddm, polar, ksvd, primme, slicot, trlan]
# List of packages in alphabetical order
externalpackages = sorted(externalwithdeps, key=lambda p: p.packagename.upper())

petscpackages    = [lapack, elemental, feast, scalapack]
specialpackages  = [slepc, petsc, slepc4py, sowing]
checkpackages    = specialpackages + petscpackages + externalwithdeps

# Print help if requested and check for wrong command-line options
if showhelp:
  print('\nConfiguration script for SLEPc '+slepc.version)
  print('\nUsage: ./configure [OPTION]...\n')
  print('  Brackets indicate an optional part')
  print('  <bool> means a boolean, use either 0 or 1')
  print('  <dir> means a directory')
  print('  <fname> means a file name, can also include the full path or url')
  print('  <libraries> means a quoted list of libraries, e.g., --with-arpack-lib="-lparpack -larpack"')
  print('  <flags> means a string of flags, e.g., --download-primme-cflags="-std=c99 -g"')
  for pkg in specialpackages:
    pkg.ShowHelp()
  print('\nOptional packages via PETSc (these are tested by default if present in PETSc\'s configuration):\n')
  for pkg in petscpackages:
    pkg.ShowHelp()
  print('\nOptional packages (external):\n')
  for pkg in externalpackages:
    pkg.ShowHelp()
  print('')
  sys.exit(0)
argdb.ErrorPetscOptions()
argdb.ErrorIfNotEmpty()

# Check if packages-download directory contains requested packages
if slepc.downloaddir:
  l = list(filter(None, [pkg.MissingTarball(slepc.downloaddir) for pkg in externalpackages + [sowing]]))
  if l:
    log.Println('\n\nDownload the following packages and run the script again:')
    for pkg in l: log.Println(pkg)
    log.Exit('Missing files in packages-download directory')

# Create directories for configuration files
archdir, archdirexisted = slepc.CreateDirTest(slepc.dir,petsc.archname)
libdir  = slepc.CreateDir(archdir,'lib')
confdir = slepc.CreateDirTwo(libdir,'slepc','conf')

# Open log file
log.Open(slepc.dir,confdir,'configure.log')
log.write('='*80)
log.write('Starting Configure Run at '+time.ctime(time.time()))
log.write('Configure Options: '+' '.join(sys.argv[1:]))
log.write('Working directory: '+os.getcwd())
log.write('Python version:\n'+sys.version)
log.write('make: '+petsc.make)

# Some checks related to PETSc configuration
if petsc.nversion < slepc.nversion:
  log.Exit('This SLEPc version is not compatible with PETSc version '+petsc.version)
if not petsc.precision in ['double','single','__float128']:
  log.Exit('This SLEPc version does not work with '+petsc.precision+' precision')

# Display versions and paths
log.write('PETSc source directory: '+petsc.dir)
log.write('PETSc install directory: '+petsc.prefixdir)
log.write('PETSc version: '+petsc.lversion)
if not petsc.isinstall:
  log.write('PETSc architecture: '+petsc.arch)
log.write('SLEPc source directory: '+slepc.dir)
if slepc.isinstall:
  log.write('SLEPc install directory: '+slepc.prefixdir)
log.write('SLEPc version: '+slepc.lversion)

# Clean previous configuration if needed
if archdirexisted:
  if slepc.isinstall and not slepc.clean:
    log.Warn('You are requesting a prefix install but the arch directory '+archdir+' already exists and may contain files from previous builds; consider adding option --with-clean')
  try:
    os.unlink(os.path.join(confdir,'files'))
  except: pass
  if slepc.clean:
    log.Println('\nCleaning arch dir '+archdir+'...')
    try:
      for root, dirs, files in os.walk(archdir,topdown=False):
        for name in files:
          if name!='configure.log':
            os.remove(os.path.join(root,name))
    except:
      log.Exit('Cannot remove existing files in '+archdir)
    for rdir in ['obj','externalpackages']:
      try:
        shutil.rmtree(os.path.join(archdir,rdir))
      except: pass

# Generate/check configure hash file
configurehash = slepc.GetConfigureHash(argdb,petsc)
hashfile = os.path.join(confdir,'configure-hash')
log.write('Checking configure hashfile')
if slepc.force:
  ResetConfigureHash(hashfile,log)
else:
  # compare hash with existing configure-hash file
  a = ''
  try:
    with open(hashfile,'r') as f:
      a = f.read()
  except:
    log.write('No previous hashfile found')
  if a == configurehash:
    log.Println('\n')
    log.write('configure hash file '+hashfile+' matches; no need to run configure.')
    log.Println('Your configure options and state have not changed; no need to rerun configure.')
    log.Println('However you can force a configure run using the option: --force')
    Epilog(slepc,petsc)
    sys.exit(0)

# Write main configuration files
if not slepc.prefixdir:
  slepc.prefixdir = archdir
includedir = slepc.CreateDir(archdir,'include')
with slepc.CreateFile(confdir,'slepcvariables') as slepcvars:
  with slepc.CreateFile(confdir,'slepcrules') as slepcrules:
    with slepc.CreateFile(includedir,'slepcconf.h') as slepcconf:
      for pkg in checkpackages:
        pkg.Process(slepcconf,slepcvars,slepcrules,slepc,petsc,archdir)
      slepcconf.write('#define SLEPC_HAVE_PACKAGES ":')
      for pkg in petscpackages + externalpackages:
        if hasattr(pkg,'havepackage') and pkg.havepackage: slepcconf.write(pkg.packagename+':')
      slepcconf.write('"\n#endif\n')
      libflags = []
      includeflags = []
      for pkg in externalwithdeps:
        if hasattr(pkg,'havepackage') and pkg.havepackage:
          for entry in pkg.libflags.split():
            if entry not in libflags:
               libflags.append(entry)
          if hasattr(pkg,'includeflags'):
            for entry in pkg.includeflags.split():
              if entry not in includeflags:
                 includeflags.append(entry)
      slepcvars.write('SLEPC_EXTERNAL_LIB = '+' '.join(libflags)+'\n')
      slepcvars.write('SLEPC_EXTERNAL_INCLUDES = '+' '.join(includeflags)+'\n')
      slepcvars.write('SLEPC_EXTERNAL_LIB_BASIC = '+' '.join(list(set(libflags).difference(set(slepc.libflags.split()))))+'\n')
      slepcvars.write('SLEPC_EXTERNAL_INCLUDES_BASIC = '+' '.join(list(set(includeflags).difference(set(slepc.includeflags.split()))))+'\n')

log.NewSection('Writing various configuration files...')

# Create global configuration file for the case of empty PETSC_ARCH
if petsc.isinstall:
  with slepc.CreateFile(os.path.join(slepc.dir,'lib','slepc','conf'),'slepcvariables') as globconf:
    globconf.write('SLEPC_DIR = '+slepc.dir+'\n')
    globconf.write('PETSC_ARCH = '+petsc.archname+'\n')

# Write Modules configuration file
modulesdir = slepc.CreateDirTwo(confdir,'modules','slepc')
log.write('Modules file in '+modulesdir)
with slepc.CreateFile(modulesdir,slepc.lversion) as modules:
  WriteModulesFile(modules,slepc.lversion,slepc.prefixdir if slepc.isinstall else slepc.dir)

# Write pkg-config configuration file
pkgconfdir = slepc.CreateDir(libdir,'pkgconfig')
log.write('pkg-config file in '+pkgconfdir)
for pkfile in ['SLEPc.pc','slepc.pc']:
  with slepc.CreateFile(pkgconfdir,pkfile) as pkgconfig:
    WritePkgconfigFile(pkgconfig,slepc.lversion,petsc.version,slepc.dir,slepc.isinstall,slepc.prefixdir,petsc.singlelib,petsc.lib_name_suffix)

# Write reconfigure file
if not slepc.isinstall:
  log.write('Reconfigure file in '+confdir)
  with slepc.CreateFile(confdir,'reconfigure-'+petsc.archname+'.py') as reconfig:
    WriteReconfigScript(reconfig,slepc.dir,argdb.UsedArgs())
  try:
    os.chmod(os.path.join(confdir,'reconfigure-'+petsc.archname+'.py'),0o775)
  except OSError as e:
    log.Exit('Unable to make reconfigure script executable:\n'+str(e))

# Configure successful, write configure-hash file
ResetConfigureHash(hashfile,log)
try:
  with open(hashfile,'w') as f:
    f.write(configurehash)
except:
  log.write('Error when writing configure-hash file')

# Print summary
log.NewSection('')
log.Println('')
log.Println('='*80)
log.Println('SLEPc Configuration')
log.Println('='*80)
for pkg in checkpackages:
  pkg.ShowInfo()
log.write('\nFinishing Configure Run at '+time.ctime(time.time()))
log.write('='*80)
log.Close()
Epilog(slepc,petsc)
