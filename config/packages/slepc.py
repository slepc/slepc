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
import argdb, os, sys, package

class SLEPc(package.Package):

  def __init__(self,argdb,log):
    self.log         = log
    self.ProcessArgs(argdb)
    self.isinstall   = not self.prefixdir==''

  def ShowHelp(self):
    wd = package.Package.wd
    print('\nConfiguration:')
    print('  --help, -h'.ljust(wd)+': Display this help and exit')
    print('  --with-clean=<bool>'.ljust(wd)+': Delete prior build files including externalpackages')
    print('  --force=<bool>'.ljust(wd)+': Bypass configure hash caching, and run to completion')
    print('  --with-packages-download-dir=<dir>'.ljust(wd)+': Skip network download of tarballs and locate them in specified dir')
    print('  --with-packages-build-dir=<dir>'.ljust(wd)+': Location to unpack and run the build process for downloaded packages')
    print('\nSLEPc:')
    print('  --prefix=<dir>'.ljust(wd)+': Specify location to install SLEPc (e.g., /usr/local)')
    print('  --DATAFILESPATH=<dir>'.ljust(wd)+': Location of datafiles (available at https://slepc.upv.es/datafiles)')

  def ProcessArgs(self,argdb):
    self.clean       = argdb.PopBool('with-clean')[0]
    self.force       = argdb.PopBool('force')[0]
    self.datadir     = argdb.PopPath('DATAFILESPATH',exist=True)[0]
    self.downloaddir = argdb.PopPath('with-packages-download-dir',exist=True)[0]
    self.pkgbuilddir = argdb.PopPath('with-packages-build-dir',exist=True)[0]
    self.prefixdir   = argdb.PopPath('prefix')[0]
    if self.prefixdir:
      self.prefixdir = os.path.realpath(os.path.normpath(self.prefixdir))

  def Process(self,slepcconf,slepcvars,slepcrules,slepc,petsc,archdir=''):
    slepcvars.write('include '+petsc.petscvariables+'\n')
    slepcvars.write('SLEPC_CONFIGURE_OPTIONS = '+' '.join(sys.argv[1:])+'\n')
    slepcvars.write('SLEPC_INSTALLDIR = '+slepc.prefixdir+'\n')
    if petsc.isinstall:
      slepcvars.write('INSTALLED_PETSC = 1\n')
    if slepc.datadir:
      slepcvars.write('DATAFILESPATH = '+slepc.datadir+'\n')
    # This should be the first write to file slepcconf.h
    slepcconf.write('#if !defined(INCLUDED_SLEPCCONF_H)\n')
    slepcconf.write('#define INCLUDED_SLEPCCONF_H\n\n')
    self.AddDefine(slepcconf,'PETSC_DIR',petsc.dir)
    self.AddDefine(slepcconf,'PETSC_ARCH',petsc.arch)
    self.AddDefine(slepcconf,'DIR',slepc.dir)
    self.AddDefine(slepcconf,'LIB_DIR',os.path.join(slepc.prefixdir,'lib'))
    libdir = os.path.join(slepc.prefixdir,'lib')
    incdir = os.path.join(slepc.prefixdir,'include')
    if petsc.buildsharedlib:
      self.libflags = petsc.cc_linker_slflag + libdir + ' -L' + libdir
    else:
      self.libflags = '-L' + libdir
    self.includeflags = '-I' + incdir
    if slepc.isrepo and slepc.gitrev != 'unknown':
      self.AddDefine(slepcconf,'VERSION_GIT',slepc.gitrev)
      self.AddDefine(slepcconf,'VERSION_DATE_GIT',slepc.gitdate)
      self.AddDefine(slepcconf,'VERSION_BRANCH_GIT',slepc.branch)
    # Single library installation
    if petsc.singlelib:
      slepcvars.write('SHLIBS = libslepc${LIB_NAME_SUFFIX}\n')
      slepcvars.write('LIBNAME = '+os.path.join('${INSTALL_LIB_DIR}','libslepc${LIB_NAME_SUFFIX}.${AR_LIB_SUFFIX}')+'\n')
      for module in ['SYS','EPS','SVD','PEP','NEP','MFN','LME']:
        slepcvars.write('SLEPC_'+module+'_LIB_NOPETSC = ${CC_LINKER_SLFLAG}${SLEPC_LIB_DIR} -L${SLEPC_LIB_DIR} -lslepc${LIB_NAME_SUFFIX} ${SLEPC_EXTERNAL_LIB_BASIC}\n')
      slepcvars.write('SLEPC_LIB_NOPETSC = ${CC_LINKER_SLFLAG}${SLEPC_LIB_DIR} -L${SLEPC_LIB_DIR} -lslepc${LIB_NAME_SUFFIX} ${SLEPC_EXTERNAL_LIB_BASIC}\n')

  def ShowInfo(self):
    self.log.Println('\nSLEPc directory:\n  '+self.dir)
    if self.isrepo and self.gitrev != 'unknown':
      self.log.Println('  It is a git repository on branch: '+self.branch)
    if self.isinstall:
      self.log.Println('SLEPc prefix directory:\n  '+self.prefixdir)

  def InitDir(self):
    if 'SLEPC_DIR' in os.environ:
      self.dir = os.path.normpath(os.environ['SLEPC_DIR'])
      if not os.path.exists(self.dir) or not os.path.exists(os.path.join(self.dir,'config')):
        self.log.Exit('SLEPC_DIR environment variable is not valid')
      if os.path.realpath(os.getcwd()) != os.path.realpath(self.dir):
        self.log.Exit('SLEPC_DIR is not the current directory')
    else:
      self.dir = os.getcwd()
      if not os.path.exists(os.path.join(self.dir,'config')):
        self.log.Exit('Current directory is not valid')

  def LoadVersion(self):
    try:
      with open(os.path.join(self.dir,'include','slepcversion.h')) as f:
        for l in f.readlines():
          l = l.split()
          if len(l) == 3:
            if l[1] == 'SLEPC_VERSION_RELEASE':
              self.release = l[2]
            if l[1] == 'SLEPC_VERSION_MAJOR':
              major = l[2]
            elif l[1] == 'SLEPC_VERSION_MINOR':
              minor = l[2]
            elif l[1] == 'SLEPC_VERSION_SUBMINOR':
              subminor = l[2]
      if self.release=='0': subminor = '99'
      self.version = major + '.' + minor
      self.lversion = major + '.' + minor + '.' + subminor
      self.nversion = int(major)*100 + int(minor)
    except:
      self.log.Exit('File error while reading SLEPc version')

    # Check whether this is a working copy of the repository
    self.isrepo = False
    if os.path.exists(os.path.join(self.dir,'src','docs')):
      self.log.write('This appears to be a repository clone - src/docs exists')
      self.isrepo = True
      if os.path.exists(os.path.join(self.dir,'.git')):
        self.log.write('.git directory exists')
        (status, output) = self.RunCommand('git help')
        if status:
          self.log.Warn('SLEPC_DIR appears to be a git working copy, but git is not found in PATH')
          self.gitrev = 'unknown'
          self.gitdate = 'unknown'
          self.branch = 'unknown'
        else:
          (status, self.gitrev) = self.RunCommand('git describe --match=v*')
          if not self.gitrev:
            (status, self.gitrev) = self.RunCommand('git log -1 --pretty=format:%H')
          (status, self.gitdate) = self.RunCommand('git log -1 --pretty=format:%ci')
          (status, self.branch) = self.RunCommand('git describe --contains --all HEAD')
          if status:
            (status, self.branch) = self.RunCommand('git rev-parse --abbrev-ref HEAD')
            if status:
              self.branch = 'unknown'
      else:
        self.log.write('This repository clone is obtained as a tarball as no .git dirs exist')
        self.gitrev = 'unknown'
        self.gitdate = 'unknown'
        self.branch = 'unknown'
    else:
      if os.path.exists(os.path.join(self.dir,'.git')):
        self.log.Exit('Your SLEPc source tree is broken. Use "git status" to check, or remove the entire directory and start all over again')
      else:
        self.log.write('This is a tarball installation')

  def CreateFile(self,basedir,fname):
    ''' Create file basedir/fname and return path string '''
    newfile = os.path.join(basedir,fname)
    try:
      newfile = open(newfile,'w')
    except:
      self.log.Exit('Cannot create '+fname+' file in '+basedir)
    return newfile

  def CreateDir(self,basedir,dirname):
    ''' Create directory basedir/dirname and return path string '''
    newdir = os.path.join(basedir,dirname)
    if not os.path.exists(newdir):
      try:
        os.mkdir(newdir)
      except:
        self.log.Exit('Cannot create '+dirname+' directory: '+newdir)
    return newdir

  def CreateDirTwo(self,basedir,dir1,dir2):
    ''' Create directory basedir/dir1/dir2 and return path string '''
    newbasedir = os.path.join(basedir,dir1)
    if not os.path.exists(newbasedir):
      try:
        os.mkdir(newbasedir)
      except:
        self.log.Exit('Cannot create '+dir1+' directory: '+newbasedir)
    newdir = os.path.join(newbasedir,dir2)
    if not os.path.exists(newdir):
      try:
        os.mkdir(newdir)
      except:
        self.log.Exit('Cannot create '+dir2+' directory: '+newdir)
    return newdir

  def CreateDirTest(self,basedir,dirname):
    ''' Create directory, return path string and flag indicating if already existed '''
    newdir = os.path.join(basedir,dirname)
    if not os.path.exists(newdir):
      existed = False
      try:
        os.mkdir(newdir)
      except:
        self.log.Exit('Cannot create '+dirname+' directory: '+newdir)
    else:
      existed = True
    return newdir, existed

  def CreatePrefixDirs(self,prefixdir):
    ''' Create directories include and lib under prefixdir, and return path strings '''
    if not os.path.exists(prefixdir):
      try:
        os.mkdir(prefixdir)
      except:
        self.log.Exit('Cannot create prefix directory: '+prefixdir)
    incdir = os.path.join(prefixdir,'include')
    if not os.path.exists(incdir):
      try:
        os.mkdir(incdir)
      except:
        self.log.Exit('Cannot create include directory: '+incdir)
    libdir = os.path.join(prefixdir,'lib')
    if not os.path.exists(libdir):
      try:
        os.mkdir(libdir)
      except:
        self.log.Exit('Cannot create lib directory: '+libdir)
    return incdir,libdir

  def GetExternalPackagesDir(self,archdir):
    ''' Create externalpackages if needed, unless user specified --with-packages-build-dir '''
    externdir = self.pkgbuilddir
    if not externdir:
      externdir = os.path.join(archdir,'externalpackages')
      if not os.path.exists(externdir):
        try:
          os.mkdir(externdir)
        except:
          self.log.Exit('Cannot create directory: '+externdir)
    return externdir

  def AddDefine(self,conffile,name,value,prefix='SLEPC_'):
    conffile.write('#define '+prefix+name+' "'+value+'"\n')

  def GetConfigureHash(self,argdb,petsc):
    ''' Generate string to be saved as configure-hash to uniquely identify the current configure run '''
    import hashlib
    import platform
    hash = 'Uname: '+platform.uname().system+' '+platform.uname().processor+'\n'
    hash += 'PATH=' + os.environ.get('PATH','') + '\n'
    hash += str(petsc) + '\n'
    args = sorted(set(filter(lambda x: (x != '--force'),argdb.UsedArgsList())))
    hash += 'SLEPc configure options:\n' + '\n'.join('    '+a for a in args) + '\n'
    chash = ''
    try:
      for root, dirs, files in os.walk('config'):
        for f in files:
          if not f.endswith('.py') or f.startswith('.') or f.startswith('#'):
            continue
          fname = os.path.join(root,f)
          with open(fname,'rb') as f:
            chash += hashlib.sha256(f.read()).hexdigest() + '  ' + fname + '\n'
    except:
      self.log.Warn('Error generating file list/hash from config directory for configure hash, forcing new configuration')
      return ''
    hash += '\n'.join(sorted(chash.splitlines()))
    return hash

