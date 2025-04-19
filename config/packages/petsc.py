#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import package, os, sys

class PETSc(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename = 'petsc'

  def __str__(self):
    ''' String conversion, used when writing the configure-hash file '''
    if not hasattr(self,'dir'): return ''
    conf = 'PETSC_DIR=' + self.dir + '\n'
    if not hasattr(self,'lversion'): return conf
    conf += 'PETSc version: ' + self.lversion + '\n'
    if not hasattr(self,'arch'): return conf
    conf += 'PETSC_ARCH=' + self.arch + '\n'
    conf += 'PETSC_SCALAR=' + self.scalar + '\n'
    conf += 'PETSC_PRECISION=' + self.precision + '\n'
    conf += 'BLASLAPACK_LIB=' + self.blaslapack_lib + '\n'
    conf += 'CC=' + self.cc + '\n'
    conf += 'CC_FLAGS=' + self.cc_flags + '\n'
    if hasattr(self,'fc'):
      conf += 'FC=' + self.fc + '\n'
      conf += 'FC_FLAGS=' + (self.fc_flags if hasattr(self,'fc_flags') else '') + '\n'
    if hasattr(self,'cxx'):
      conf += 'CXX=' + self.cxx + '\n'
      conf += 'CXX_FLAGS=' + (self.cxx_flags if hasattr(self,'cxx_flags') else '') + '\n'
    conf += 'PETSc configure options:\n'
    args = sorted(set(self.configure_options.split()))
    conf += '\n'.join('    '+a for a in args) + '\n'
    return conf

  def ShowHelp(self):
    pass

  def Check(self):
    (result, output) = self.Link([],[],'')
    self.havepackage = result

  def InitDir(self,prefixdir):
    if 'PETSC_DIR' in os.environ:
      self.dir = os.path.normpath(os.environ['PETSC_DIR'])
      if not os.path.exists(self.dir):
        self.log.Exit('PETSC_DIR environment variable is not valid')
    else:
      if prefixdir:
        self.dir = prefixdir
        os.environ['PETSC_DIR'] = self.dir
      else:
        self.log.Exit('PETSC_DIR environment variable is not set')

  def LoadVersion(self):
    try:
      f = open(os.path.join(self.dir,'include','petscversion.h'))
      for l in f.readlines():
        l = l.split()
        if len(l) == 3:
          if l[1] == 'PETSC_VERSION_RELEASE':
            self.release = l[2]
          if l[1] == 'PETSC_VERSION_MAJOR':
            major = l[2]
          elif l[1] == 'PETSC_VERSION_MINOR':
            minor = l[2]
          elif l[1] == 'PETSC_VERSION_SUBMINOR':
            subminor = l[2]
      f.close()
      if self.release=='0': subminor = '99'
      self.version = major + '.' + minor
      self.lversion = major + '.' + minor + '.' + subminor
      self.nversion = int(major)*100 + int(minor)
    except:
      self.log.Exit('File error while reading PETSc version')

    # Check whether this is a working copy of the repository
    self.isrepo = False
    if os.path.exists(os.path.join(self.dir,'.git')):
      (status, output) = self.RunCommand('cd '+self.dir+';git rev-parse')
      if not status:
        self.isrepo = True
        (status, self.gitrev)  = self.RunCommand('cd '+self.dir+';git log -1 --pretty=format:%H')
        (status, self.gitdate) = self.RunCommand('cd '+self.dir+';git log -1 --pretty=format:%ci')
        (status, self.branch)  = self.RunCommand('cd '+self.dir+';git describe --contains --all HEAD')

  def LoadConf(self):
    if 'PETSC_ARCH' in os.environ and os.environ['PETSC_ARCH']:
      self.isinstall = False
      self.arch = os.environ['PETSC_ARCH']
      if os.path.basename(self.arch) != self.arch:
        suggest = os.path.basename(self.arch)
        if not suggest: suggest = os.path.basename(self.arch[0:-1])
        self.log.Exit('Variable PETSC_ARCH must not be a full path\nYou set PETSC_ARCH=%s, maybe you meant PETSC_ARCH=%s'% (self.arch,suggest))
      self.petscvariables = os.path.join(self.dir,self.arch,'lib','petsc','conf','petscvariables')
      petscconf_h = os.path.join(self.dir,self.arch,'include','petscconf.h')
    else:
      self.isinstall = True
      self.petscvariables = os.path.join(self.dir,'lib','petsc','conf','petscvariables')
      petscconf_h = os.path.join(self.dir,'include','petscconf.h')

    self.buildsharedlib = False
    self.lib_name_suffix = ''
    try:
      with open(self.petscvariables) as f:
        for l in f.readlines():
          r = l.split('=',1)
          if len(r)!=2: continue
          k = r[0].strip()
          v = r[1].strip()
          v = v.replace('${PETSC_DIR}',self.dir)  # needed in some Cray installations
          if k == 'PETSC_SCALAR':
            self.scalar = v
          elif k == 'PETSC_PRECISION':
            self.precision = v
          elif k == 'FC' and not v=='':
            self.fc = v
          elif k == 'BUILDSHAREDLIB' and v=='yes':
            self.buildsharedlib = True
          else:
            if k in ['AR','AR_FLAGS','AR_LIB_SUFFIX','BFORT','BLASLAPACK_LIB','CC','CC_FLAGS','CC_LINKER_SLFLAG','CMAKE','CONFIGURE_OPTIONS','CPP','CXX','CXX_FLAGS','FC_FLAGS','FC_VERSION','LIB_NAME_SUFFIX','MAKE','MAKE_NP','PREFIXDIR','RANLIB','SCALAPACK_LIB','SEDINPLACE','SL_LINKER_SUFFIX']:
              setattr(self,k.lower(),v)
    except:
      self.log.Exit('Cannot process file ' + self.petscvariables)

    self.ind64 = False
    self.mpiuni = False
    self.msmpi = False
    self.debug = False
    self.singlelib = False
    self.blaslapackmangling = ''
    self.blaslapackint64 = False
    self.fortran = False
    self.language = 'c'
    self.maxcxxdialect = ''
    self.packages = []
    try:
      with open(petscconf_h) as f:
        for l in f.readlines():
          l = l.split()
          if len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_64BIT_INDICES' and l[2]=='1':
            self.ind64 = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_MPIUNI' and l[2]=='1':
            self.mpiuni = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_MSMPI' and l[2]=='1':
            self.msmpi = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_DEBUG' and l[2]=='1':
            self.debug = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_SINGLE_LIBRARY' and l[2]=='1':
            self.singlelib = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_BLASLAPACK_UNDERSCORE' and l[2]=='1':
            self.blaslapackmangling = 'underscore'
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_BLASLAPACK_CAPS' and l[2]=='1':
            self.blaslapackmangling = 'caps'
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_64BIT_BLAS_INDICES' and l[2]=='1':
            self.blaslapackint64 = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_USE_FORTRAN_BINDINGS' and l[2]=='1':
            self.fortran = True
          elif len(l)==3 and l[0]=='#define' and l[1]=='PETSC_CLANGUAGE_CXX' and l[2]=='1':
            self.language = 'c++'
          elif self.isinstall and len(l)==3 and l[0]=='#define' and l[1]=='PETSC_ARCH':
            self.arch = l[2].strip('"')
          else:
            for p in ['elemental','hpddm','mkl_libs','mkl_includes','mkl_pardiso','scalapack','slepc']:
              if len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_'+p.upper() and l[2]=='1':
                self.packages.append(p)
            for p in ['20','17','14','11']:
              if len(l)==3 and l[0]=='#define' and l[1]=='PETSC_HAVE_CXX_DIALECT_CXX'+p and l[2]=='1' and (self.maxcxxdialect=='' or p>self.maxcxxdialect):
                self.maxcxxdialect = p
                break
      if 'mkl_libs' in self.packages or 'mkl_includes' in self.packages:
        if 'mkl_libs' in self.packages: self.packages.remove('mkl_libs')
        if 'mkl_includes' in self.packages: self.packages.remove('mkl_includes')
        self.packages.append('mkl')
    except:
      if self.isinstall:
        self.log.Exit('Cannot process file ' + petscconf_h + ', maybe you forgot to set PETSC_ARCH')
      else:
        self.log.Exit('Cannot process file ' + petscconf_h)

    if self.isinstall:
      pseudoarch = 'arch-' + sys.platform.replace('cygwin','mswin')+ '-' + self.language.lower().replace('+','x')
      if self.debug:
        pseudoarch += '-debug'
      else:
        pseudoarch += '-opt'
      if not 'real' in self.scalar:
        pseudoarch += '-' + self.scalar
      self.archname = 'installed-'+pseudoarch.replace('linux-','linux2-')
    else:
      self.archname = self.arch

  def Process(self,slepcconf,slepcvars,slepcrules,slepc,petsc,archdir=''):
    self.log.NewSection('Checking PETSc installation...')
    if petsc.nversion > slepc.nversion:
      self.log.Warn('PETSc version '+petsc.version+' is newer than SLEPc version '+slepc.version)
    if slepc.release=='1' and not petsc.release=='1':
      errmsg = 'A release version of SLEPc requires a release version of PETSc, not a development version'
      if self.isrepo and slepc.isrepo:
        errmsg += '\nType "git checkout release" (or "git checkout main") in both PETSc and SLEPc repositories'
      self.log.Exit(errmsg)
    if slepc.release=='0' and petsc.release=='1':
      if not 'slepc' in petsc.packages:
        errmsg = 'A development version of SLEPc cannot be built with a release version of PETSc'
        if self.isrepo and slepc.isrepo:
          errmsg += '\nType "git checkout release" (or "git checkout main") in both PETSc and SLEPc repositories'
        self.log.Exit(errmsg)
    if petsc.isinstall:
      if os.path.realpath(petsc.prefixdir) != os.path.realpath(petsc.dir):
        self.log.Warn('PETSC_DIR does not point to PETSc installation path')
    petsc.Check()
    if not petsc.havepackage:
      self.log.Exit('Unable to link with PETSc')
    if slepc.isrepo and petsc.isrepo and not petsc.isinstall and petsc.branch!='release' and slepc.branch!='release':
      try:
        import dateutil.parser, datetime
        petscdate = dateutil.parser.parse(petsc.gitdate)
        slepcdate = dateutil.parser.parse(slepc.gitdate)
        if abs(petscdate-slepcdate)>datetime.timedelta(days=30):
          self.log.Warn('Your PETSc and SLEPc repos may not be in sync (more than 30 days apart)')
      except ImportError: pass

  def ShowInfo(self):
    self.log.Println('PETSc directory:\n  '+self.dir)
    if self.isrepo:
      self.log.Println('  It is a git repository on branch: '+self.branch)
    if self.isinstall:
      self.log.Println('Prefix install with '+self.precision+' precision '+self.scalar+' numbers')
    else:
      self.log.Println('Architecture "'+self.archname+'" with '+self.precision+' precision '+self.scalar+' numbers')

  def isGfortran100plus(self):
    '''returns true if the compiler is gfortran-10.0.x or later'''
    try:
      (result, output) = self.RunCommand(self.fc+' --version',)
      import re
      strmatch = re.match(r'GNU Fortran\s+\(.*\)\s+(\d+)\.(\d+)',output)
      if strmatch:
        VMAJOR,VMINOR = strmatch.groups()
        if (int(VMAJOR),int(VMINOR)) >= (10,0):
          return 1
    except RuntimeError:
      pass

  def isIntel(self):
    '''returns true if the compiler is an Intel compiler'''
    try:
      (result, output) = self.RunCommand(self.cc+' --help | head -n 80',)
      if 'Intel(R)' in output:
        return 1
    except RuntimeError:
      pass

  def removeWarningFlags(self,flags):
    outflags = []
    for flag in flags:
      if not flag in ['-Werror','-Wall','-Wwrite-strings','-Wno-strict-aliasing','-Wno-unknown-pragmas','-Wno-unused-variable','-Wno-unused-dummy-argument','-fvisibility=hidden','-std=c89','-pedantic','--coverage','-Mfree','-fdefault-integer-8']:
        outflags.append(flag)
    return outflags

  def getCFlags(self):
    outflags = self.removeWarningFlags(self.cc_flags.split())
    return ' '.join(outflags)

  def getCXXFlags(self):
    if hasattr(self,'cxx_flags'):
      outflags = self.removeWarningFlags(self.cxx_flags.split())
      return ' '.join(outflags)
    else:
      return ''

  def getFFlags(self):
    if hasattr(self,'fc_flags'):
      outflags = self.removeWarningFlags(self.fc_flags.split())
      if self.isGfortran100plus():
        outflags.append('-fallow-argument-mismatch')
      return ' '.join(outflags)
    else:
      return ''

