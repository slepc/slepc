#!/usr/bin/env python3
from __future__ import print_function
import re, os, sys, shutil
import subprocess

try:
  WindowsError
except NameError:
  WindowsError = None

class Installer:
  def __init__(self, args = None):
    if len(args)<6:
      print('********************************************************************')
      print('Installation script error - not enough arguments:')
      print('./config/install.py SLEPC_DIR PETSC_DIR SLEPC_INSTALLDIR -destDir=DESTDIR [-no-examples] PETSC_ARCH LIB_SUFFIX RANLIB')
      print('********************************************************************')
      sys.exit(1)
    self.rootDir     = args[0]
    self.petscDir    = args[1]
    self.destDir     = os.path.abspath(args[2])
    if args[3].startswith('-destDir='):
      if args[3][9:]:
        # prepend to prefix dir
        self.destDir   = os.path.join(args[3][9:], os.path.relpath(self.destDir, os.sep))
    else:
      print('********************************************************************')
      print('Error in argument, should start with -destDir=')
      print('********************************************************************')
      sys.exit(1)
    self.copyexamples = True
    if args[4] == '-no-examples':
      self.copyexamples = False
      del args[4]
    self.arch        = args[4]
    self.arLibSuffix = args[5]
    self.ranlib      = ' '.join(args[6:])
    self.copies = []
    return

  def readInstallDir(self, src):
    try:
      f = open(src)
      for l in f.readlines():
        r = l.split('=',1)
        if len(r)!=2: continue
        if r[0].strip() == 'SLEPC_INSTALLDIR':
          break
      f.close()
    except:
      print('********************************************************************')
      print('Error reading SLEPC_INSTALLDIR from slepcvariables')
      print('********************************************************************')
      sys.exit(1)
    return r[1].strip()

  def readPetscCC(self, src):
    try:
      f = open(src)
      for l in f.readlines():
        r = l.split('=',1)
        if len(r)!=2: continue
        if r[0].strip() == 'CC':
          break
      f.close()
    except:
      print('********************************************************************')
      print('Error reading CC from petscvariables')
      print('********************************************************************')
      sys.exit(1)
    return r[1].strip()

  def setupDirectories(self):
    self.archDir           = os.path.join(self.rootDir, self.arch)
    self.rootIncludeDir    = os.path.join(self.rootDir, 'include')
    self.archIncludeDir    = os.path.join(self.rootDir, self.arch, 'include')
    self.rootConfDir       = os.path.join(self.rootDir, 'lib','slepc','conf')
    self.archConfDir       = os.path.join(self.rootDir, self.arch, 'lib','slepc','conf')
    self.rootBinDir        = os.path.join(self.rootDir, 'lib','slepc','bin')
    self.archBinDir        = os.path.join(self.rootDir, self.arch, 'bin')
    self.archLibDir        = os.path.join(self.rootDir, self.arch, 'lib')
    self.destIncludeDir    = os.path.join(self.destDir, 'include')
    self.destConfDir       = os.path.join(self.destDir, 'lib','slepc','conf')
    self.destLibDir        = os.path.join(self.destDir, 'lib')
    self.destBinDir        = os.path.join(self.destDir, 'lib','slepc','bin')
    self.installDir        = self.readInstallDir(os.path.join(self.archConfDir,'slepcvariables'))
    self.installIncludeDir = os.path.join(self.installDir, 'include')
    self.installBinDir     = os.path.join(self.installDir, 'lib','slepc','bin')
    self.rootShareDir      = os.path.join(self.rootDir, 'share')
    self.destShareDir      = os.path.join(self.destDir, 'share')
    self.rootSrcDir        = os.path.join(self.rootDir, 'src')
    arch = '' if self.arch.startswith('installed-') else self.arch
    self.petscConfDir      = os.path.join(self.petscDir, arch, 'lib', 'petsc', 'conf')
    self.petscCC           = self.readPetscCC(os.path.join(self.petscConfDir,'petscvariables'))
    return

  def checkDestdir(self):
    if os.path.exists(self.destDir):
      if os.path.samefile(self.destDir, self.rootDir):
        print('********************************************************************')
        print('Incorrect prefix usage. Specified destDir same as current SLEPC_DIR')
        print('********************************************************************')
        sys.exit(1)
      if os.path.samefile(self.destDir, os.path.join(self.rootDir,self.arch)):
        print('********************************************************************')
        print('Incorrect prefix usage. Specified destDir same as current SLEPC_DIR/PETSC_ARCH')
        print('********************************************************************')
        sys.exit(1)
      if not os.path.isdir(os.path.realpath(self.destDir)):
        print('********************************************************************')
        print('Specified destDir', self.destDir, 'is not a directory. Cannot proceed!')
        print('********************************************************************')
        sys.exit(1)
      if not os.access(self.destDir, os.W_OK):
        print('********************************************************************')
        print('Unable to write to ', self.destDir, 'Perhaps you need to do "sudo make install"')
        print('********************************************************************')
        sys.exit(1)
    return

  def copyfile(self, src, dst, symlinks = False, copyFunc = shutil.copy2):
    """Copies a single file    """
    copies = []
    errors = []
    if not os.path.exists(dst):
      os.makedirs(dst)
    elif not os.path.isdir(dst):
      raise shutil.Error('Destination is not a directory')
    srcname = src
    dstname = os.path.join(dst, os.path.basename(src))
    try:
      if symlinks and os.path.islink(srcname):
        linkto = os.readlink(srcname)
        os.symlink(linkto, dstname)
      else:
        copyFunc(srcname, dstname)
        copies.append((srcname, dstname))
    except (IOError, os.error) as why:
      errors.append((srcname, dstname, str(why)))
    except shutil.Error as err:
      errors.append((srcname,dstname,str(err.args[0])))
    if errors:
      raise shutil.Error(errors)
    return copies

  def copyexamplefiles(self, src, dst, copyFunc = shutil.copy2):
    """Copies all files, but not directories in a single file    """
    names  = os.listdir(src)
    for name in names:
      if not name.endswith('.html'):
        srcname = os.path.join(src, name)
        if os.path.isfile(srcname):
           self.copyfile(srcname,dst)

  def fixExamplesMakefile(self, src):
    '''Change ./${PETSC_ARCH} in makefile in root slepc directory with ${SLEPC_DIR}'''
    lines   = []
    oldFile = open(src, 'r')
    alllines=oldFile.read()
    oldFile.close()
    newlines=alllines.split('\n')[0]+'\n'  # Firstline
    # Hardcode SLEPC_DIR, PETSC_DIR and PETSC_ARCH to avoid users doing the wrong thing
    newlines+='SLEPC_DIR='+self.installDir+'\n'
    newlines+='PETSC_DIR='+self.petscDir+'\n'
    newlines+='PETSC_ARCH=\n'
    for line in alllines.split('\n')[1:]:
      if line.startswith('TESTLOGTAPFILE'):
        newlines+='TESTLOGTAPFILE = $(TESTDIR)/test_install_tap.log\n'
      elif line.startswith('TESTLOGERRFILE'):
        newlines+='TESTLOGERRFILE = $(TESTDIR)/test_install_err.log\n'
      elif line.startswith('$(generatedtest)') and 'slepcvariables' in line:
        newlines+='all: test\n\n'+line+'\n'
      else:
        newlines+=line+'\n'
    newFile = open(src, 'w')
    newFile.write(newlines)
    newFile.close()
    return

  def copyExamples(self, src, dst):
    """Recursively copy the examples directories
    """
    if not os.path.isdir(dst):
      raise shutil.Error('Destination is not a directory')

    self.copyfile(os.path.join(src,'makefile'),dst)
    names  = os.listdir(src)
    nret2 = 0
    for name in names:
      srcname = os.path.join(src, name)
      dstname = os.path.join(dst, name)
      if not name.startswith('arch') and os.path.isdir(srcname) and os.path.isfile(os.path.join(srcname,'makefile')):
        os.mkdir(dstname)
        nret = self.copyExamples(srcname,dstname)
        if name == 'tests' or name == 'tutorials' or name == 'nlevp' or name == 'cnetwork':
          self.copyexamplefiles(srcname,dstname)
          if os.path.isdir(os.path.join(srcname,'output')):
            os.mkdir(os.path.join(dstname,'output'))
            self.copyexamplefiles(os.path.join(srcname,'output'),os.path.join(dstname,'output'))
          nret = 1
        if not nret:
          # prune directory branches that don't have examples under them
          os.unlink(os.path.join(dstname,'makefile'))
          os.rmdir(dstname)
        nret2 = nret + nret2
    return nret2

  def copyConfig(self, src, dst):
    """Recursively copy the examples directories
    """
    if not os.path.isdir(dst):
      raise shutil.Error('Destination is not a directory')

    self.copies.extend(self.copyfile('gmakefile.test',dst))
    #newConfigDir=os.path.join(dst,'config')  # Am not renaming at present
    #if not os.path.isdir(newConfigDir): os.mkdir(newConfigDir)
    #testConfFiles="".split()
    #for tf in testConfFiles:
    #  self.copies.extend(self.copyfile(os.path.join('config',tf),newConfigDir))
    #return

  def copytree(self, src, dst, symlinks = False, copyFunc = shutil.copy2, exclude = [], exclude_ext= ['.DSYM','.o','.pyc','.html'], recurse = 1):
    """Recursively copy a directory tree using copyFunc, which defaults to shutil.copy2().

       The copyFunc() you provide is only used on the top level, lower levels always use shutil.copy2

    The destination directory must not already exist.
    If exception(s) occur, an shutil.Error is raised with a list of reasons.

    If the optional symlinks flag is true, symbolic links in the
    source tree result in symbolic links in the destination tree; if
    it is false, the contents of the files pointed to by symbolic
    links are copied.
    """
    copies = []
    names  = os.listdir(src)
    if not os.path.exists(dst):
      os.makedirs(dst)
    elif not os.path.isdir(dst):
      raise shutil.Error('Destination is not a directory')
    errors = []
    for name in names:
      srcname = os.path.join(src, name)
      dstname = os.path.join(dst, name)
      try:
        if symlinks and os.path.islink(srcname):
          linkto = os.readlink(srcname)
          os.symlink(linkto, dstname)
        elif os.path.isdir(srcname) and recurse and not os.path.basename(srcname) in exclude:
          copies.extend(self.copytree(srcname, dstname, symlinks,exclude = exclude, exclude_ext = exclude_ext))
        elif os.path.isfile(srcname) and not os.path.basename(srcname) in exclude and os.path.splitext(name)[1] not in exclude_ext:
          copyFunc(srcname, dstname)
          copies.append((srcname, dstname))
        # XXX What about devices, sockets etc.?
      except (IOError, os.error) as why:
        errors.append((srcname, dstname, str(why)))
      # catch the Error from the recursive copytree so that we can
      # continue with other files
      except shutil.Error as err:
        errors.append((srcname,dstname,str(err.args[0])))
    try:
      shutil.copystat(src, dst)
    except OSError as e:
      if WindowsError is not None and isinstance(e, WindowsError):
        # Copying file access times may fail on Windows
        pass
      else:
        errors.append((src, dst, str(e)))
    if errors:
      raise shutil.Error(errors)
    return copies


  def fixConfFile(self, src):
    lines   = []
    oldFile = open(src, 'r')
    for line in oldFile.readlines():
      # paths generated by configure could be different link-path than what's used by user, so fix both
      line = line.replace(os.path.join(self.rootDir, self.arch), self.installDir)
      line = line.replace(os.path.realpath(os.path.join(self.rootDir, self.arch)), self.installDir)
      line = line.replace(os.path.join(self.rootDir, 'bin'), self.installBinDir)
      line = line.replace(os.path.realpath(os.path.join(self.rootDir, 'bin')), self.installBinDir)
      line = line.replace(os.path.join(self.rootDir, 'include'), self.installIncludeDir)
      line = line.replace(os.path.realpath(os.path.join(self.rootDir, 'include')), self.installIncludeDir)
      line = line.replace(self.rootDir, self.installDir)
      # remove SLEPC_DIR/PETSC_ARCH variables from conf-makefiles. They are no longer necessary
      line = line.replace('${SLEPC_DIR}/${PETSC_ARCH}', self.installDir)
      line = line.replace('PETSC_ARCH=${PETSC_ARCH}', '')
      line = line.replace('${SLEPC_DIR}', self.installDir)
      lines.append(line)
    oldFile.close()
    newFile = open(src, 'w')
    newFile.write(''.join(lines))
    newFile.close()
    return

  def fixConf(self):
    import shutil
    for file in ['slepc_rules', 'slepc_rules_doc.mk', 'slepc_rules_util.mk', 'slepc_variables', 'slepcrules', 'slepcvariables']:
      self.fixConfFile(os.path.join(self.destConfDir,file))
    self.fixConfFile(os.path.join(self.destLibDir,'pkgconfig','slepc.pc'))
    self.fixConfFile(os.path.join(self.destIncludeDir,'slepcconf.h'))
    return

  def fixPythonWheel(self):
    import glob
    import shutil
    #
    for pattern in (
        self.destLibDir  + '/*.a',
        self.destLibDir  + '/*.la',
        self.destLibDir  + '/pkgconfig',  # TODO: keep?
        self.destConfDir + '/configure-hash',
        self.destConfDir + '/uninstall.py',
        self.destConfDir + '/reconfigure-*.py',
        self.destConfDir + '/pkg.conf.*',
        self.destConfDir + '/pkg.git*.*',
        self.destConfDir + '/modules',  # TODO: keep?
        self.destShareDir + '/*/examples/src/*',
        self.destShareDir + '/*/datafiles',
    ):
      for pathname in glob.glob(pattern):
        if os.path.isdir(pathname):
          shutil.rmtree(pathname)
        elif os.path.exists(pathname):
          os.remove(pathname)
    #
    for filename in (
      self.destIncludeDir + '/slepcconf.h',
      self.destShareDir + '/slepc/examples/gmakefile.test',
      self.destConfDir + '/slepc_rules_doc.mk',
      self.destConfDir + '/slepc_rules_util.mk',
      self.destConfDir + '/slepc_rules',
      self.destConfDir + '/slepcrules',
      self.destConfDir + '/slepc_variables',
      self.destConfDir + '/slepcvariables',
    ):
      with open(filename, 'r') as oldFile:
        contents = oldFile.read()
      contents = contents.replace(self.installDir, '${SLEPC_DIR}')
      contents = contents.replace(self.rootDir, '${SLEPC_DIR}')
      contents = contents.replace(self.petscDir, '${PETSC_DIR}')
      with open(filename, 'w') as newFile:
        newFile.write(contents)
    #
    def lsdir(dirname, *patterns):
      return glob.glob(os.path.join(dirname, *patterns))
    def shell(*args):
      out  = subprocess.check_output(list(args), universal_newlines=True)
      return out[:-1] if out[-1:] == '\n' else out
    plibdir = os.path.join(self.petscDir, 'lib')
    slibdir = os.path.join(self.installDir, 'lib')
    if sys.platform == 'linux':
      libraries = [
        lib for lib in lsdir(self.destLibDir, 'lib*.so*')
        if not os.path.islink(lib)
      ]
      for shlib in libraries:
        # fix shared library rpath
        rpath = shell('patchelf', '--print-rpath', shlib)
        rpath = rpath.split(os.path.pathsep)
        if plibdir in rpath:
          rpath.insert(0, '$ORIGIN/../../petsc/lib')
          while plibdir in rpath:
            rpath.remove(plibdir)
        if slibdir in rpath:
          rpath.insert(0, '$ORIGIN')
          while slibdir in rpath:
            rpath.remove(slibdir)
        if rpath:
          rpath = os.path.pathsep.join(rpath)
          shell('patchelf', '--set-rpath', rpath, shlib)
        # fix shared library file and symlink
        basename = os.path.basename(shlib)
        libname, ext, _ = basename.partition('.so')
        liblink = libname + ext
        soname = shell('patchelf', '--print-soname', shlib)
        for symlink in lsdir(self.destLibDir, liblink + '*'):
          if os.path.islink(symlink):
            os.unlink(symlink)
        curdir = os.getcwd()
        try:
          os.chdir(os.path.dirname(shlib))
          if soname != basename:
            os.rename(basename, soname)
          if soname != liblink:
            os.symlink(soname, liblink)
        finally:
          os.chdir(curdir)
    if sys.platform == 'darwin':
      def otool(cmd, dylib):
        pattern = r'''
          ^\s+ cmd \s %s$\n
          ^\s+ cmdsize \s \d+$\n
          ^\s+ (?:name|path) \s (.*) \s \(offset \s \d+\)$
        ''' % cmd
        return re.findall(
          pattern, shell('otool', '-l', dylib),
          flags=re.VERBOSE | re.MULTILINE,
        )
      libraries = [
        lib for lib in lsdir(self.destLibDir, 'lib*.dylib')
        if not os.path.islink(lib)
      ]
      for dylib in libraries:
        install_name = otool('LC_ID_DYLIB', dylib)[0]
        dependencies = otool('LC_LOAD_DYLIB', dylib)
        runtime_path = otool('LC_RPATH', dylib)
        # fix shared library install name and rpath
        install_name = '@rpath/' + os.path.basename(install_name)
        shell('install_name_tool', '-id', install_name, dylib)
        for libdir in (plibdir, slibdir):
          if libdir in runtime_path:
            shell('install_name_tool', '-delete_rpath', libdir, dylib)
        for rpath in ('@loader_path', '@loader_path/../../petsc/lib'):
          if rpath not in runtime_path:
            shell('install_name_tool', '-add_rpath', rpath, dylib)
        for dep in dependencies:
          if os.path.dirname(dep) in (plibdir, slibdir):
            newid = '@rpath/' + os.path.basename(dep)
            shell('install_name_tool', '-change', dep, newid, dylib)
        # fix shared library file and symlink
        basename = os.path.basename(dylib)
        libname, ext = os.path.splitext(basename)
        libname = libname.partition('.')[0]
        liblink = libname + ext
        dyname = os.path.basename(install_name)
        for symlink in lsdir(self.destLibDir, libname + '*' + ext):
          if os.path.islink(symlink):
            os.unlink(symlink)
        curdir = os.getcwd()
        try:
          os.chdir(os.path.dirname(dylib))
          if dyname != basename:
            os.rename(basename, dyname)
          if dyname != liblink:
            os.symlink(dyname, liblink)
        finally:
          os.chdir(curdir)
    #
    return

  def createUninstaller(self):
    uninstallscript = os.path.join(self.destConfDir, 'uninstall.py')
    f = open(uninstallscript, 'w')
    # Could use the Python AST to do this
    f.write('#!'+sys.executable+'\n')
    f.write('import os\n')

    f.write('copies = '+repr(self.copies).replace(self.destDir,self.installDir))
    f.write('''
for src, dst in copies:
  try:
    os.remove(dst)
  except:
    pass
''')
    #TODO: need to delete libXXX.YYY.dylib.dSYM directory on Mac
    dirs = [os.path.join('include','slepc','finclude'),os.path.join('include','slepc','private'),os.path.join('lib','slepc','conf')]
    newdirs = []
    for dir in dirs: newdirs.append(os.path.join(self.installDir,dir))
    f.write('dirs = '+str(newdirs))
    f.write('''
for dir in dirs:
  import shutil
  try:
    shutil.rmtree(dir)
  except:
    pass
''')
    f.close()
    os.chmod(uninstallscript,0o744)
    return

  def installIncludes(self):
    # TODO: should exclude slepc/finclude except for fortran builds
    self.copies.extend(self.copytree(self.rootIncludeDir, self.destIncludeDir,exclude = ['makefile']))
    self.copies.extend(self.copytree(self.archIncludeDir, self.destIncludeDir))
    return

  def installConf(self):
    self.copies.extend(self.copytree(self.rootConfDir, self.destConfDir, exclude = ['install.py','bfort-base.txt','bfort-mpi.txt','bfort-petsc.txt','bfort-slepc.txt']))
    self.copies.extend(self.copytree(self.archConfDir, self.destConfDir, exclude = ['configure.log','error.log','files','gmake.log','make.log','test.log','memoryerror.log']))
    return

  def installBin(self):
    #if os.path.exists(self.rootBinDir):
    #  self.copies.extend(self.copytree(self.rootBinDir, self.destBinDir))
    #if os.path.exists(self.archBinDir):
    #  self.copies.extend(self.copytree(self.archBinDir, self.destBinDir))
    return

  def installShare(self):
    self.copies.extend(self.copytree(self.rootShareDir, self.destShareDir))
    if self.copyexamples:
      examplesdir=os.path.join(self.destShareDir,'slepc','examples')
      if os.path.exists(examplesdir):
        shutil.rmtree(examplesdir)
      os.mkdir(examplesdir)
      os.mkdir(os.path.join(examplesdir,'src'))
      self.copyExamples(self.rootSrcDir,os.path.join(examplesdir,'src'))
      self.copyConfig(self.rootDir,examplesdir)
      self.fixExamplesMakefile(os.path.join(examplesdir,'gmakefile.test'))
    return

  def copyLib(self, src, dst):
    '''Run ranlib on the destination library if it is an archive. Also run install_name_tool on dylib on Mac'''
    # Symlinks (assumed local) are recreated at dst
    if os.path.islink(src):
      linkto = os.readlink(src)
      try:
        os.remove(dst)            # In case it already exists
      except OSError:
        pass
      os.symlink(linkto, dst)
      return
    # Do not install object files
    if not os.path.splitext(src)[1] == '.o':
      shutil.copy2(src, dst)
    if os.path.splitext(dst)[1] == '.'+self.arLibSuffix:
      if not 'win32fe' in self.petscCC:
        (result, output) = subprocess.getstatusoutput(self.ranlib+' '+dst)
    if os.path.splitext(dst)[1] == '.dylib' and shutil.which('otool') and shutil.which('install_name_tool'):
      output = subprocess.check_output(['otool', '-D', src], universal_newlines=True)
      oldname = output.splitlines()[1]
      installName = oldname.replace(os.path.realpath(self.archDir), self.installDir)
      subprocess.check_output(['install_name_tool', '-id', installName, dst])
    # preserve the original timestamps - so that the .a vs .so time order is preserved
    shutil.copystat(src,dst)
    return

  def installLib(self):
    self.copies.extend(self.copytree(self.archLibDir, self.destLibDir, copyFunc = self.copyLib, exclude = ['.DIR'],recurse = 0))
    self.copies.extend(self.copytree(os.path.join(self.archLibDir,'pkgconfig'), os.path.join(self.destLibDir,'pkgconfig'), copyFunc = self.copyLib, exclude = ['.DIR'],recurse = 0))
    return


  def outputInstallDone(self):
    arch=self.arch
    if arch.startswith('installed-'): arch='""'
    print('''\
====================================
Install complete.
Now to check if the libraries are working do (in current directory):
make SLEPC_DIR=%s PETSC_DIR=%s PETSC_ARCH=%s check
====================================\
''' % (self.installDir,self.petscDir,arch))
    return

  def outputDestDirDone(self):
    print('''\
====================================
Copy to DESTDIR %s is now complete.
Before use - please copy/install over to specified prefix: %s
====================================\
''' % (self.destDir,self.installDir))
    return

  def runsetup(self):
    self.setupDirectories()
    self.checkDestdir()
    return

  def runcopy(self):
    if self.destDir == self.installDir:
      print('*** Installing SLEPc at prefix location:',self.destDir, ' ***')
    else:
      print('*** Copying SLEPc to DESTDIR location:',self.destDir, ' ***')
    if not os.path.exists(self.destDir):
      try:
        os.makedirs(self.destDir)
      except:
        print('********************************************************************')
        print('Unable to create', self.destDir, 'Perhaps you need to do "sudo make install"')
        print('********************************************************************')
        sys.exit(1)
    self.installIncludes()
    self.installConf()
    self.installBin()
    self.installLib()
    self.installShare()
    self.createUninstaller()
    return

  def runfix(self):
    self.fixConf()
    using_build_backend = any(
      os.environ.get(prefix + '_BUILD_BACKEND')
      for prefix in ('_PYPROJECT_HOOKS', 'PEP517')
    )
    if using_build_backend:
      self.fixPythonWheel()
    return

  def rundone(self):
    if self.destDir == self.installDir:
      self.outputInstallDone()
    else:
      self.outputDestDirDone()
    return

  def run(self):
    self.runsetup()
    self.runcopy()
    self.runfix()
    self.rundone()
    return

if __name__ == '__main__':
  Installer(sys.argv[1:]).run()
