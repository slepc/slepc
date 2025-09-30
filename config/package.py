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
import os, sys, tempfile, shutil, tarfile
import log, argdb
from urllib.request import urlretrieve
from urllib import parse as urlparse_local
import subprocess
import socket
from shutil import which  # just to break compatibility with python2

# Fix parsing for nonstandard schemes
urlparse_local.uses_netloc.extend(['bk', 'ssh', 'svn'])

class Package:

  packages = []    # list of packages processed and added so far

  def __init__(self,argdb,log):
    self.installable     = False  # an already installed package can be picked --with-xxx-dir
    self.downloadable    = False  # package can be downloaded and installed with --download-xxx
    self.downloadpackage = 0
    self.packagetype     = ''     # can be 'gnu', 'cmake', 'source_c', or empty
    self.packagedir      = ''
    self.packagelibs     = ''
    self.packageincludes = ''
    self.packageurl      = ''
    self.buildflags      = ''
    self.log             = log
    self.supportsscalar  = ['real', 'complex']
    self.supportssingle  = False
    self.supports64bint  = False
    self.fortran         = False
    self.hasheaders      = False
    self.requested       = False
    self.havepackage     = False

  def RunCommand(self,instr):
    try:
      self.log.write('- '*35+'\nRunning command:\n'+instr+'\n'+'- '*35)
    except AttributeError: pass
    try:
      output = subprocess.check_output(instr,shell=True,stderr=subprocess.STDOUT)
      result = 0
    except subprocess.CalledProcessError as ex:
      output = ex.output
      result = ex.returncode
    output = output.decode(encoding='UTF-8',errors='replace').rstrip()
    try:
      self.log.write('Output:\n'+output+'\n'+'- '*35)
    except AttributeError: pass
    return (result,output)

  def ProcessArgs(self,argdb,petscpackages=''):
    if hasattr(self,'petscdepend') and self.petscdepend in petscpackages:
      self.requested = True
    if self.installable and not hasattr(self,'petscdepend'):
      string,found = argdb.PopPath('with-'+self.packagename+'-dir',exist=True)
      if found:
        self.requested = True
        self.packagedir = string
      string,found = argdb.PopString('with-'+self.packagename+'-lib')
      if found:
        if self.packagedir:
          self.log.Exit('Specify either "--with-'+self.packagename+'-dir" or "--with-'+self.packagename+'-lib%s", but not both!' % (' --with-'+self.packagename+'-include' if self.hasheaders else ''))
        self.requested = True
        self.packagelibs = string
      if self.hasheaders:
        string,found = argdb.PopString('with-'+self.packagename+'-include')
        if found:
          self.requested = True
          self.packageincludes = string
    if self.installable:
      value,found = argdb.PopBool('with-'+self.packagename)
      if found:
        self.requested = value
    if self.downloadable:
      flagsfound = False
      if self.packagetype == 'gnu':
        string,flagsfound = argdb.PopString('download-'+self.packagename+'-configure-arguments')
      elif self.packagetype == 'cmake':
        string,flagsfound = argdb.PopString('download-'+self.packagename+'-cmake-arguments')
      elif self.packagetype == 'source_c':
        string,flagsfound = argdb.PopString('download-'+self.packagename+'-cflags')
      url,flag,found = argdb.PopUrl('download-'+self.packagename)
      if found:
        if self.requested:
          self.log.Exit('Cannot request both download and install simultaneously')
        self.requested = True
        self.download = True
        self.packageurl = url
        self.downloadpackage = flag
      if flagsfound:
        if not getattr(self,'download',False):
          if self.packagetype == 'gnu':
            self.log.Exit('--download-'+self.packagename+'-configure-arguments must be used together with --download-'+self.packagename)
          elif self.packagetype == 'cmake':
            self.log.Exit('--download-'+self.packagename+'-cmake-arguments must be used together with --download-'+self.packagename)
          elif self.packagetype == 'source_c':
            self.log.Exit('--download-'+self.packagename+'-cflags must be used together with --download-'+self.packagename)
        self.buildflags = string

  def Process(self,slepcconf,slepcvars,slepcrules,slepc,petsc,archdir=''):
    self.make = petsc.make
    if petsc.buildsharedlib:
      self.slflag = petsc.cc_linker_slflag
    if self.requested:
      name = self.packagename.upper()
      if self.downloadpackage:
        if hasattr(self,'version') and self.packageurl=='':
          self.log.NewSection('Installing '+name+' version '+self.version+'...')
        else:
          self.log.NewSection('Installing '+name+'...')
        self.Precondition(slepc,petsc)
        self.DownloadAndInstall(slepcconf,slepcvars,slepc,petsc,archdir,slepc.prefixdir)
      elif self.installable:
        self.log.NewSection('Checking '+name+'...')
        self.Precondition(slepc,petsc)
        if petsc.buildsharedlib:
          self.packagelibs = self.DistilLibList(self.packagelibs,petsc)
        self.Check(slepcconf,slepcvars,petsc,archdir)
        if not self.havepackage: self.log.setLastFailed()
      try:
        self.LoadVersion(slepcconf)
        self.log.write('Version number for '+name+' is '+self.iversion)
      except AttributeError:
        pass
      if self.havepackage: Package.packages.append(self)
    else: # not requested
      if hasattr(self,'SkipInstall'):
        self.SkipInstall(slepcrules)

  def Precondition(self,slepc,petsc):
    package = self.packagename.upper()
    if petsc.scalar == 'complex':
      if 'complex' not in self.supportsscalar:
        self.log.Exit(package+' does not support complex scalars')
    elif petsc.scalar == 'real':
      if 'real' not in self.supportsscalar:
        self.log.Exit(package+' is supported only with complex scalars')
    if petsc.precision == 'single':
      if not self.supportssingle:
        self.log.Exit(package+' is supported only in double precision')
    elif petsc.precision != 'double':
      self.log.Exit('Precision '+petsc.precision+' is not supported for external packages')
    if petsc.ind64 and not self.supports64bint:
      self.log.Exit(package+' cannot be used with 64-bit integers')
    if self.downloadpackage and self.fortran and not hasattr(petsc,'fc'):
      self.log.Exit('Option --download-'+self.packagename+' requires a Fortran compiler')

  def Require(self,packagename):
    for p in Package.packages:
      if p.packagename.upper() == packagename.upper():
        return p
    self.log.Exit('The package '+self.packagename.upper()+' requires configuring also with '+packagename.upper()+". Run configure --help for details on how to install it")

  def DistilLibList(self,packagelibs,petsc):
    libs = []
    for l in packagelibs.split():
      if l.endswith(petsc.sl_linker_suffix):
        filename = os.path.basename(l)
        libs.append('-L'+l.rstrip(filename))
        libs.append(self.slflag+l.rstrip(filename))
        libs.append('-l'+filename.lstrip('lib').rstrip('.'+petsc.sl_linker_suffix))
      else:
        libs.append(l)
    newldflags = []
    newlibs = []
    dupflags = ['-L',self.slflag]
    for j in libs:
      # remove duplicate -L, -Wl,-rpath options - and only consecutive -l options
      if j in newldflags and any([j.startswith(flg) for flg in dupflags]): continue
      if newlibs and j == newlibs[-1]: continue
      if list(filter(j.startswith,['-l'])) or list(filter(j.endswith,['.lib','.a','.so','.o'])) or j in ['-Wl,-Bstatic','-Wl,-Bdynamic','-Wl,--start-group','-Wl,--end-group']:
        newlibs.append(j)
      else:
        newldflags.append(j)
    liblist = ' '.join(newldflags + newlibs)
    return liblist

  def GetArchiveName(self):
    '''Return name of archive after downloading'''
    if self.packageurl=='':
      archivename = self.archive
    else:
      parsed = urlparse_local.urlparse(self.packageurl)
      archivename = os.path.basename(parsed[2])
    if archivename[0] == 'v':
      archivename = archivename[1:]
    try:
      if archivename[0].isdigit() or int(archivename.split('.')[0],16):
        archivename = self.packagename+'-'+archivename
    except: pass
    return archivename

  def GetDirectoryName(self):
    '''Return name of the directory after extracting the tarball'''
    dirname = self.GetArchiveName()
    for suffix in ('.tar.gz','.tgz'):
      if dirname.endswith(suffix):
        dirname = dirname[:-len(suffix)]
    return dirname

  def MissingTarball(self,downloaddir):
    '''Check if tarball is missing in downloaddir'''
    if self.downloadable and getattr(self,'download',False):
      localFile = os.path.join(downloaddir,self.GetArchiveName())
      if not os.path.exists(localFile):
        url = self.packageurl
        if url=='':
          url = self.url
        return self.packagename+': '+url+' --> '+localFile

  def Download(self,externdir,downloaddir):
    # Quick return: check if source is already available
    if os.path.exists(os.path.join(externdir,self.GetDirectoryName())):
      self.log.write('Using '+os.path.join(externdir,self.GetDirectoryName()))
      return os.path.join(externdir,self.GetDirectoryName())

    if downloaddir:
      # Get tarball from download dir
      localFile = os.path.join(downloaddir,self.GetArchiveName())
      if not os.path.exists(localFile):
        self.log.Exit('Could not find file '+self.GetArchiveName()+' under '+downloaddir)
      url = localFile
      filename = os.path.basename(url)
    else:
      # Download tarball
      url = self.packageurl
      if url=='':
        url = self.url
      if os.path.exists(url):
        url = 'file:'+url
      filename = os.path.basename(urlparse_local.urlparse(url)[2])
      localFile = os.path.join(externdir,self.GetArchiveName())
      self.log.write('Downloading '+url+' to '+localFile)

      if os.path.exists(localFile):
        os.remove(localFile)
      try:
        sav_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(30)
        urlretrieve(url, localFile)
        socket.setdefaulttimeout(sav_timeout)
      except Exception as e:
        socket.setdefaulttimeout(sav_timeout)
        failureMessage = '''\
Unable to download package %s from: %s
* If URL specified manually - perhaps there is a typo?
* If your network is disconnected - please reconnect and rerun ./configure
* Or perhaps you have a firewall blocking the download
* You can run with --with-packages-download-dir=/adirectory and ./configure will instruct you what packages to download manually
* or you can download the above URL manually, to /yourselectedlocation/%s
  and use the configure option:
  --download-%s=/yourselectedlocation/%s
''' % (self.packagename.upper(), url, filename, self.packagename, filename)
        self.log.Exit(failureMessage)

    # Uncompress tarball
    extractdir = os.path.join(externdir,self.GetDirectoryName())
    self.log.write('Uncompressing '+localFile+' to directory '+extractdir)
    if os.path.exists(extractdir):
      for root, dirs, files in os.walk(extractdir, topdown=False):
        for name in files:
          os.remove(os.path.join(root,name))
        for name in dirs:
          os.rmdir(os.path.join(root,name))
    failureMessage = '''\
Downloaded package %s from: %s is not a tarball.
[or installed python cannot process compressed files]
* If you are behind a firewall - please fix your proxy and rerun ./configure
  For example at LANL you may need to set the environmental variable http_proxy (or HTTP_PROXY?) to  http://proxyout.lanl.gov
* You can run with --with-packages-download-dir=/adirectory and ./configure will instruct you what packages to download manually
* or you can download the above URL manually, to /yourselectedlocation/%s
  and use the configure option:
  --download-%s=/yourselectedlocation/%s
''' % (self.packagename.upper(), url, filename, self.packagename, filename)
    try:
      tf = tarfile.open(localFile)
    except tarfile.ReadError as e:
      self.log.Exit(str(e)+'\n'+failureMessage)
    else:
      if not tf: self.log.Exit(failureMessage)
      with tf:
        #git puts 'pax_global_header' as the first entry and some tar utils process this as a file
        firstname = tf.getnames()[0]
        if firstname == 'pax_global_header':
          firstmember = tf.getmembers()[1]
        else:
          firstmember = tf.getmembers()[0]
        # some tarfiles list packagename/ but some list packagename/filename in the first entry
        if firstmember.isdir():
          dirname = firstmember.name
        else:
          dirname = os.path.dirname(firstmember.name)
        tf.extractall(path=externdir)

    # fix file permissions for the untared tarballs
    try:
      # check if 'dirname' is set'
      if dirname:
        (result,output) = self.RunCommand('cd '+externdir+'; chmod -R a+r '+dirname+'; find '+dirname+r' -type d -name "*" -exec chmod a+rx {} \;')
        if dirname != self.GetDirectoryName():
          self.log.write('The directory name '+dirname+' is different from the expected one, renaming to '+self.GetDirectoryName())
          os.rename(os.path.join(externdir,dirname),os.path.join(externdir,self.GetDirectoryName()))
          dirname = self.GetDirectoryName()
      else:
        self.log.Warn('Could not determine dirname extracted by '+localFile+' to fix file permissions')
    except RuntimeError as e:
      self.log.Exit('Error changing permissions for '+dirname+' obtained from '+localFile+ ' : '+str(e))
    if not downloaddir:
      os.remove(localFile)
    return os.path.join(externdir,dirname)

  wd = 36

  def ShowHelp(self):
    wd = Package.wd
    if self.downloadable or self.installable:
      print(self.packagename.upper()+':')
    if self.downloadable:
      print(('  --download-'+self.packagename+'[=<fname>]').ljust(wd)+': Download and install '+self.packagename.upper())
      if self.packagetype == 'gnu':
        print(('  --download-'+self.packagename+'-configure-arguments=<flags>').ljust(wd)+': Indicate extra flags to configure '+self.packagename.upper())
      elif self.packagetype == 'cmake':
        print(('  --download-'+self.packagename+'-cmake-arguments=<flags>').ljust(wd)+': Indicate extra flags to build '+self.packagename.upper()+' with CMake')
      elif self.packagetype == 'source_c':
        print(('  --download-'+self.packagename+'-cflags=<flags>').ljust(wd)+': Indicate extra flags to compile '+self.packagename.upper())
    if self.installable:
      print(('  --with-'+self.packagename+'=<bool>').ljust(wd)+': Test for '+self.packagename.upper()+' (requires PETSc with %s)'%getattr(self,'petscdepend','').upper())
    if self.installable and not hasattr(self,'petscdepend'):
      print(('  --with-'+self.packagename+'-dir=<dir>').ljust(wd)+': Indicate the root directory of the '+self.packagename.upper()+' installation')
      print(('  --with-'+self.packagename+'-lib=<libraries>').ljust(wd)+': Indicate quoted list of libraries and link flags for '+self.packagename.upper())
      if self.hasheaders:
        print(('  --with-'+self.packagename+'-include=<dirs>').ljust(wd)+': Indicate the directory of the '+self.packagename.upper()+' include files')

  def ShowInfo(self):
    if self.havepackage:
      if hasattr(self,'version') and self.downloadpackage and self.packageurl=='':
        packagename = self.packagename.upper()+' version '+self.version
      else:
        packagename = self.packagename.upper()
      if hasattr(self,'petscdepend'):
        self.log.Println(packagename+' from %s linked by PETSc' % self.petscdepend.upper())
      elif hasattr(self,'packageflags'):
        self.log.Println(packagename+' library flags:')
        self.log.Println('  '+self.packageflags)
      else:
        self.log.Println(packagename+' installed')

  def Link(self,functions,callbacks,flags,givencode='',cflags='',clanguage='c',logdump=True):

    # Create temporary directory and makefile
    try:
      tmpdir = tempfile.mkdtemp(prefix='slepc-')
      if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    except:
      self.log.Exit('Cannot create temporary directory')
    try:
      with open(os.path.join(tmpdir,'makefile'),'w') as makefile:
        makefile.write('checklink: checklink.o\n')
        makefile.write('\t${CLINKER} -o checklink checklink.o ${LINKFLAGS} ${PETSC_SNES_LIB}\n')
        makefile.write('\t@${RM} -f checklink checklink.o\n')
        makefile.write('include '+os.path.join('${PETSC_DIR}','lib','petsc','conf','variables')+'\n')
        makefile.write('include '+os.path.join('${PETSC_DIR}','lib','petsc','conf','rules')+'\n')
        if cflags:
          if clanguage=='c++': makefile.write('CXXFLAGS='+cflags+'\n')
          else: makefile.write('CFLAGS='+cflags+'\n')
    except:
      self.log.Exit('Cannot create makefile in temporary directory')

    # Create source file
    if givencode == '':
      code = '#include "petscsnes.h"\n'
      for f in functions:
        code += 'PETSC_EXTERN int\n' + f + '();\n'

      for c in callbacks:
        code += 'int '+ c + '() { return 0; } \n'

      code += 'int main() {\n'
      code += 'PetscErrorCode ierr; Vec v; Mat m; KSP k;\n'
      code += 'ierr = PetscInitializeNoArguments();\n'
      code += 'ierr = VecCreate(PETSC_COMM_WORLD,&v);\n'
      code += 'ierr = MatCreate(PETSC_COMM_WORLD,&m);\n'
      code += 'ierr = KSPCreate(PETSC_COMM_WORLD,&k);\n'
      code += '(void)ierr;\n'
      for f in functions:
        code += f + '();\n'
      code += 'return 0;\n}\n'
    else:
      code = givencode

    with open(os.path.join(tmpdir,'checklink.cxx' if clanguage=='c++' else 'checklink.c'),'w') as cfile:
      cfile.write(code)
    if logdump:
      try:
        self.log.write('- '*35+'\nChecking link with code:\n')
        self.log.write(code)
      except AttributeError: pass

    # Try to compile test program
    (result, output) = self.RunCommand('cd ' + tmpdir + ';' + self.make + ' checklink LINKFLAGS="'+flags+'"')
    shutil.rmtree(tmpdir)

    if result:
      return (0,code + output)
    else:
      return (1,code + output)

  def FortranLink(self,functions,callbacks,flags):
    f = []
    for i in functions:
      f.append(i+'_')
    c = []
    for i in callbacks:
      c.append(i+'_')
    (result, output1) = self.Link(f,c,flags,logdump=False)
    output1 = '\n====== With underscore Fortran names\n' + output1
    if result: return ('UNDERSCORE',output1)

    f = []
    for i in functions:
      f.append(i.upper())
    c = []
    for i in callbacks:
      c.append(i.upper())
    (result, output2) = self.Link(f,c,flags,logdump=False)
    output2 = '\n====== With capital Fortran names\n' + output2
    if result: return ('CAPS',output2)

    output = '\n=== With linker flags: '+flags
    return ('',output + output1 + output2)

  def GenerateGuesses(self,name,archdir,word='lib'):
    installdirs = [os.path.join(os.path.sep,'usr','local'),os.path.join(os.path.sep,'opt')]
    if 'HOME' in os.environ:
      installdirs.insert(0,os.environ['HOME'])

    dirs = []
    for i in installdirs:
      dirs = dirs + [os.path.join(i,word)]
      for d in [name,name.upper(),name.lower()]:
        dirs = dirs + [os.path.join(i,d)]
        dirs = dirs + [os.path.join(i,d,word)]
        dirs = dirs + [os.path.join(i,word,d)]

    for d in dirs[:]:
      if not os.path.exists(d):
        dirs.remove(d)
    dirs = [''] + dirs + [os.path.join(archdir,word)]
    return dirs

  def FortranLib(self,slepcconf,slepcvars,dirs,libs,functions,callbacks = []):
    name = self.packagename.upper()
    mangling = ''
    if isinstance(libs, str): # user-provided string with link options
      flags = libs
      (mangling, output) = self.FortranLink(functions,callbacks,flags)
      error = output
    else:
      error = ''
      flags = ''
      for d in dirs:
        for l in libs:
          if d:
            if hasattr(self,'slflag'):
              flags = ' '.join([self.slflag + d] + ['-L' + d] + l)
            else:
              flags = ' '.join(['-L' + d] + l)
          else:
            flags = ' '.join(l)
          (mangling, output) = self.FortranLink(functions,callbacks,flags)
          error += output
          if mangling: break
        if mangling: break

    if mangling:
      self.log.write(output)
    else:
      self.log.write(error)
      self.log.Exit('Unable to link with '+name+' library in directories '+' '.join(dirs)+' with libraries and link flags '+flags)

    slepcconf.write('#define SLEPC_HAVE_' + name + ' 1\n#define SLEPC_' + name + '_HAVE_'+mangling+' 1\n')
    if not hasattr(self,'petscdepend'):
      slepcvars.write(name + '_LIB = '+flags+'\n')
      self.libflags = flags
      self.packageflags = flags
    self.havepackage = True

  def WriteMakefile(self,fname,builddir,cont):
    self.log.write('Using makefile:\n')
    self.log.write(cont)
    with open(os.path.join(builddir,fname),'w') as mfile:
      mfile.write(cont)

  def DefaultIncludePath(self,petsc,file):
    (result,output) = self.RunCommand('echo | '+petsc.cpp+' -Wp,-v -')
    if not result:
      import re
      dirs = re.findall('^ .*',output,re.MULTILINE)
      for s in dirs:
        d = s[1:]
        if os.path.isfile(os.path.join(d,file)):
          self.log.write('Found '+os.path.join(d,file))
          return d
    return '/usr/include'

