#!/usr/bin/python
#
#
# Check ABI/API compatibility of two SLEPc versions (the old and the new). Old and new SLEPc should have already been built using GCC and with -g
#
# Usage:
#   ./abicheck.py -old_dir        <old SLEPC_DIR>
#                 -old_arch       <old PETSC_ARCH>
#                 -old_petsc_dir  <old PETSC_DIR>
#                 -new_dir        <new SLEPC_DIR>
#                 -new_arch       <new PETSC_ARCH>
#                 -new_petsc_dir  <new PETSC_DIR>
#                 -report_format  <format of the report file, either xml or html. Optional with default=html>
#

import os, sys
import subprocess
import argparse

# Copy from config/BuildSystem/config/setCompilers.py
def isGNU(compiler):
  '''Returns true if the compiler is a GNU compiler'''
  try:
    output = subprocess.check_output([compiler,'--help'], stderr=subprocess.STDOUT)
    return (any([s in output for s in ['www.gnu.org',
                                       'bugzilla.redhat.com',
                                       'gcc.gnu.org',
                                       'gcc version',
                                       '-print-libgcc-file-name',
                                       'passed on to the various sub-processes invoked by gcc',
                                       'passed on to the various sub-processes invoked by cc',
                                       'passed on to the various sub-processes invoked by gfortran',
                                       'passed on to the various sub-processes invoked by g++',
                                       'passed on to the various sub-processes invoked by c++',
                                       ]])
            and not any([s in output for s in ['Intel(R)',
                                               'Unrecognised option --help passed to ld', # NAG f95 compiler
                                               ' clang '
                                               ]]))
  except:
    pass
  return 0

def gen_xml_desc(slepc_dir, petsc_arch, petsc_dir, libs, xmlfile):
  '''Generate an ABI descriptor file for the given SLEPC_DIR and PETSC_ARCH'''

  # Get branch and hash of the source used to build the library. We use them as an informative version number of the library
  slepcconf = os.path.join(slepc_dir, petsc_arch, 'include', 'slepcconf.h')

  branch = subprocess.check_output(['grep', '^#define SLEPC_VERSION_BRANCH_GIT', slepcconf])
  branch = branch.split("\"")[1]
  hash   = subprocess.check_output(['grep', '^#define SLEPC_VERSION_GIT', slepcconf])
  hash   = hash.split("\"")[1]

  slepclibs = ''
  slepcheaders = ''

  for name in libs:
    name = name.lower()
    if not name.startswith('slepc'): # Support slepc, slepcsys, or sys etc
      name = 'slepc'+name
    libname       = 'lib'+name
    libname_so    = os.path.join(slepc_dir, petsc_arch, 'lib', libname+'.so')
    libname_dylib = os.path.join(slepc_dir, petsc_arch, 'lib', libname+'.dylib')
    libname_a     = os.path.join(slepc_dir, petsc_arch, 'lib', libname+'.a')
    lib = ''
    if os.path.isfile(libname_so):
      lib = libname_so
    elif os.path.isfile(libname_dylib):
      lib = libname_dylib
    elif os.path.isfile(libname_a):
      lib = libname_a

    if lib == '':
      raise RuntimeError('Could not find library %s for SLEPC_DIR=%s PETSC_ARCH=%s. Please configure and build it before doing ABI/API checking\n'% (name,slepc_dir,petsc_arch))
    slepclibs += lib+'\n'

    header = os.path.join(slepc_dir,'include',name+'.h')
    if not os.path.isfile(header):
      raise RuntimeError('File %s does not exist.\n'% (header))
    slepcheaders += header+'\n'

  with open(xmlfile, "w") as file:
    file.write("<version>\n")
    file.write(branch + ' (' + hash +')')
    file.write("</version>\n\n")

    file.write("<headers>\n")
    file.write(slepcheaders)
    file.write("\n")
    file.write("</headers>\n\n")

    file.write("<libs>\n")
    file.write(slepclibs)
    file.write("\n")
    file.write("</libs>\n\n")

    file.write('<include_paths>\n')
    file.write(os.path.join(slepc_dir,'include'))
    file.write("\n")
    file.write(os.path.join(slepc_dir,petsc_arch,'include'))
    file.write("\n")
    file.write(os.path.join(petsc_dir,'include'))
    file.write("\n")
    file.write(os.path.join(petsc_dir,petsc_arch,'include'))
    file.write("\n")
    file.write("</include_paths>\n")

def run_abi_checker(slepc_dir, petsc_arch, petsc_dir, abi_dir, oldxml, newxml, cc, rformat):
  if rformat == 'html':
    report = 'report.html'
  else:
    report = 'report.xml'
  reportfile   = os.path.join(abi_dir, report)

  abichecker     = os.path.join(petsc_dir, 'lib', 'petsc', 'bin', 'maint', 'abi-compliance-checker', 'abi-compliance-checker.pl')
  ierr = subprocess.call([abichecker, '-l', 'slepc', '--lang=C', '-old', oldxml, '-new', newxml, '--gcc-path', cc, '--report-path', reportfile, '--report-format', rformat])
  return ierr

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-old_dir",  help="Old SLEPC_DIR, which defines the SLEPc library to compare with", required=True)
  parser.add_argument("-old_arch", help="Old PETSC_ARCH, which defines the SLEPc library to compare with", required=True)
  parser.add_argument("-old_petsc_dir",  help="Old PETSC_DIR", required=True)
  parser.add_argument("-new_dir",  help="New SLEPC_DIR", required=True)
  parser.add_argument("-new_arch", help="New PETSC_ARCH", required=True)
  parser.add_argument("-new_petsc_dir",  help="New PETSC_DIR", required=True)
  parser.add_argument("-report_format", help="Format of the report file", default='html', required=False)

  args      = parser.parse_args()
  old_dir   = args.old_dir
  old_arch  = args.old_arch
  old_petsc = args.new_petsc_dir
  new_dir   = args.new_dir
  new_arch  = args.new_arch
  new_petsc = args.new_petsc_dir
  rformat   = args.report_format

  # Get the compiler and raise error if it is not GNU
  cc = ''
  for [slepc_dir, petsc_arch, petsc_dir] in [[old_dir, old_arch, old_petsc],[new_dir, new_arch, new_petsc]]:
    petscvariables = os.path.join(petsc_dir, petsc_arch, 'lib', 'petsc','conf', 'petscvariables')
    if not os.path.isfile(petscvariables):
      raise RuntimeError('File %s does not exist.\n'% (petscvariables))
    ccstring = subprocess.check_output(['grep', '^CC =', petscvariables])
    cc       = ccstring.split('=')[1].strip()
    if not isGNU(cc):
      raise RuntimeError('The compiler %s is not a GNU compiler.\n'% (cc))

  # Check the PETSC_USE_SINGLE_LIBRARY status
  petscconf_h = os.path.join(new_petsc, new_arch, 'include', 'petscconf.h')
  if not os.path.isfile(petscconf_h):
    raise RuntimeError('File %s does not exist.\n'% (petscconf_h))

  if 'PETSC_USE_SINGLE_LIBRARY' in open(petscconf_h).read():
    libs=['slepc']
  else:
    libs=['sys', 'eps', 'svd', 'pep', 'nep', 'mfn', 'lme']

  # Check report format
  if rformat != 'html' and rformat != 'xml':
     raise RuntimeError('Unsupported report format "%s". Only html and xml are supported.\n'% (rformat))

  # We generate report under new_dir/abi
  abi_dir = os.path.join(new_dir, new_arch, 'abi')
  if not os.path.isdir(abi_dir):
    os.makedirs (abi_dir)

  print("\nChecking libraries %s ..." %(libs))
  oldxml = os.path.join(abi_dir, 'old.xml')
  newxml = os.path.join(abi_dir, 'new.xml')
  gen_xml_desc(old_dir,old_arch,old_petsc, libs, oldxml)
  gen_xml_desc(new_dir,new_arch,new_petsc, libs, newxml)
  ierr = run_abi_checker(new_dir, new_arch, new_petsc, abi_dir, oldxml, newxml, cc, rformat)

  print("=========================================================================================")
  if (ierr):
    print("Error: ABI/API compatibility check failed. Open the compatibility report file to see details.")
  else:
    print("ABI/API of the two versions are compatible.")

if __name__ == '__main__':
    main()
