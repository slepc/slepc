#!/usr/bin/env python3
#
#    Generates Fortran function stubs and module interface definitions for PETSc
#
#    Note 1:
#      const char *title[] gets mapped to character(*) title and the string gets copied into the space provided by the caller
#
#    This tool looks for the values MANSEC and [BFORT]SUBMANSEC (where BFORTSUBMANSEC has priority over SUBMANSEC)
#    defined in the makefile
#
#    The F90 generated interface files are stored in $PETSC_ARCH/ftn/MANSEC/slepc[BFORT]SUBMANSEC.h90
#    The Fortran stub files are stored in $PETSC_ARCH/ftn/MANSEC/**/ where ** is the directory under MANSEC of the original source
#
#    Stubs/interfaces generated from include can only involve sys files
#
#    These are then included by the src/MANSEC/ftn-mod/MANSECmod.F90 files to create the Fortran module files
#
#    SUBMANSEC (but not BFORTSUBMANSEC) is also used (in the documentation generating part of PETSc) to determine what
#    directory in doc/manualpages/ the manual pages are deposited.
#
#    An example of when the BFORTSUBMANSEC may be different than SUBMANSEC is src/dm/label/impls/ephemeral/plex where we would like
#    the documentation to be stored under DMLabel but the function interfaces need to go into the DMPLEX Fortran module
#    (not the DM Fortran module) since they depend on DMPlexTransform.
#
from __future__ import print_function
import os
import pathlib
import shutil
import sys
import subprocess
from subprocess import check_output

CToFortranTypes = {'int':'integer4', 'ptrdiff_t':'PetscInt64', 'float':'PetscFortranFloat', 'int32_t':'integer4',
                   'double':'PetscFortranDouble', 'short':None, 'size_t':'PetscSizeT', 'rocblas_status':None, 'PetscBT':None,
                   'PetscEnum':None, 'PetscDLHandle':None}

Letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','w','x','y']
verbose = False
skipinc = ['slepcversion.h']

def verbosePrint(text):
  '''Prints the text if run with verbose option'''
  if verbose: print(text)

def cross(L):
  '''Given a list of the form 'aaOOaO' generates a list of lists where the lists are all combinations of the original list with O replaced with a'''
  A = L[0]
  L = L[1:]
  if not L:
    if A == 'O':
      return [['a'],['O']]
    else:
      return [['a']]
  nL = []
  sL = cross(L)
  if A == 'O':
    for l in sL:
      k = l.copy()
      q = l.copy()
      k.insert(0,'a')
      q.insert(0,'O')
      nL.append(k)
      nL.append(q)
  else:
    for l in sL:
      l.insert(0, 'a')
      nL.append(l)
  return nL

def crossCreate(L):
  '''Given a list of arguments to a function, generates a list of combinations of the form 'aaOOaO' indicating potential optional arguments'''
  '''Currently not used because it takes too long to build the Fortran modules'''
  opts = []
  if len(L.arguments) == 0: return ['']
  for i in L.arguments:
    if i.optional: opts.append('O')
    else: opts.append('a')
  return cross(opts)

def generateFortranInterface(petscarch, classes, enums, structs, senums, funname, fun):
  '''Generates the interface definition for a function'''
  '''This is used both by class functions and standalone functions'''
  # check for functions for which we cannot build intefaces
  if fun.opaque:
    return
  if fun.name.find('_') > -1: return
  for k in fun.arguments:
    ktypename = k.typename
    if ktypename in CToFortranTypes and not CToFortranTypes[ktypename]:
      fun.opaque = True
      return
    if ktypename.find('_') and not ktypename.startswith('MPI') > -1:
      fun.opaque = True
      return
    if ktypename.endswith('Func') or ktypename.endswith('Fn') :
      fun.opaque = True
      return
    if ktypename == 'void' or ktypename == 'PeCtx':
      return

  mansec = fun.mansec
  if not mansec: mansec = fun.submansec
  file = fun.includefile + '90'
  if not file.startswith('slepc'): file = 'slepc' + file
  mansecpath = (os.path.join('sys','classes',mansec) if mansec in ['bv','ds','fn','rg','st'] else mansec)
  with open(os.path.join(petscarch,'ftn', mansecpath,file),"a") as fd:
    # Currently not used because it takes too long to build the Fortran modules
    #opts = crossCreate(fun)
    if False: # len(opts) > 1:
      for k in fun.arguments:
        if k.array and k.stars and not k.typename == 'char': return
        if k.stars and k.typename == 'MPI_Fint': return   # TODO add support for returning MPI_Fint
        if k.stars == 2 and k.typename == 'void': return

      # there are multiple interfaces for the function, they are generated in $PETSC_ARCH/ftn/MANSEC/*.hf90
      fd.write('  interface ' + funname + '\n')
      fd.write('  module procedure ' + funname)
      cnt = 0
      for fi in opts:
        if cnt > 0: fd.write(',' + funname)
        fd.write(''.join(fi))
        cnt = cnt + 1
      fd.write('\n')
    # generate the single needed interface for the function
    else:
      fd.write('  interface ' + funname + '\n')

      fi = fun
      func = ''
      dims = ['']
      petscclasses = ['PetscObject','PetscViewer','PetscRandom','PetscSubcomm','IS','Vec','Mat','KSP','PC','SNES']
      petscenums = ['NormType','MatStructure']
      # if ((funname).startswith('MatDenseGetArray') or (funname).startswith('MatDenseRestoreArray')) and fi[-1].endswith('[]'): dims = ['1d','2d']
      for dim in dims:
        fd.write('  subroutine ' + funname + func + dim + '(')
        simportset = set()
        simport = ''
        cnt = 0
        for k in fi.arguments:
          if k.stringlen: continue
          ktypename = k.typename
          if cnt: fd.write(',')
          fd.write(Letters[cnt])
          if not ktypename in simportset:
            if ktypename in classes or ktypename in petscclasses or ktypename == 'VecScatter':
              if simport: simport = simport + ','
              simport = simport + 't' + ktypename
            if ktypename in enums or ktypename in petscenums:
              if simport: simport = simport + ','
              simport = simport + 'e' + ktypename
            if ktypename in structs and not structs[ktypename].opaque:
              if simport: simport = simport + ','
              simport = simport + 's' + ktypename
          simportset.add(ktypename)
          cnt = cnt + 1
        if cnt: fd.write(',')
        fd.write(' z)\n')
        if simport: fd.write('  import ' + simport + '\n')

        cnt = 0
        for k in fun.arguments:
          if k.stringlen: continue
          ktypename = k.typename
          if ktypename in CToFortranTypes:
            ktypename =CToFortranTypes[ktypename]
          if ktypename in senums or ktypename == 'char':
            fd.write('  character(*) :: ' + Letters[cnt] + '\n')
          elif k.array and k.stars:
            if not dim or dim == '1d': fd.write('  ' + ktypename + ', pointer :: ' +  Letters[cnt]  + '(:)\n')
            else: fd.write('  ' + ktypename + ', pointer :: ' +  Letters[cnt]  + '(:,:)\n')
          elif k.array:
            fd.write('  ' + ktypename + ' :: ' +  Letters[cnt]  + '(*)\n')
          else:
            fd.write('  ' + ktypename + ' :: ' + Letters[cnt] + '\n')
          cnt = cnt + 1
        fd.write('  PetscErrorCode z\n')
        fd.write('  end subroutine\n')
    fd.write('  end interface\n')
    fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
    fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + funname + func + dim  + '\n')
    fd.write('#endif\n')

def generateCStub(petscarch,senums,classes,funname,fun):
  '''Generates the C stub that is callable from Fortran for a function'''

  #  PETSc returns strings in two ways either
  #     with a pointer to an array: char *[]
  #     or by copying the string into a given array with a given length: char [], size_t len
  #  in both cases the Fortran API passes in a character(*) with enough room to hold the result
  #  it does not pass in the string length as a separate argument
  #  The following variables are used to indicate which case is being generated
  #      self.stringlen    - True indicates the argument is the length of the previous character string
  #      self.const        - indicates the string argument is an input, not an output
  #      self.stars        - indicates the string is (in C) returned by a pointer to a string array
  if fun.opaque or fun.opaquestub: return
  # temporary
  for k in fun.arguments:
    # no C stub if function returns an array, except if it is a string
    # TODO: generate fillible stub for functions that return arrays
    if k.array and k.stars and not k.typename == 'char': return
    if k.stars and k.typename == 'MPI_Fint': return   # TODO add support for returning MPI_Fint
    if k.stars == 2 and k.typename == 'void': return

  if fun.file.endswith('.h'): return # currently ignore stubs generated from include/*.h
  dir = os.path.join(petscarch,fun.dir.replace('src/','ftn/'))
  if not os.path.isdir(dir): os.makedirs(dir)
  with open(os.path.join(dir,fun.file.replace('.c','f.c')),'a') as fd:
    fd.write('#include "petscsys.h"\n')
    fd.write('#include "petscfix.h"\n')
    fd.write('#include "petsc/private/ftnimpl.h"\n')
    fd.write('#include <slepc' + fun.mansec + '.h>\n')
    fd.write('#include <slepc' + fun.includefile.replace('slepc','') + '>\n')

    suffix = ''
    # not used because generating the Fortran modules takes too long
    #for k in fun.arguments:
    #  if k.optional:
    #    suffix = 'raw'

    fd.write('#if defined(PETSC_HAVE_FORTRAN_CAPS)\n')
    fd.write('  #define ' + (funname + suffix).lower() + '_ ' + (funname + suffix).upper() + '\n')
    fd.write('#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)\n')
    fd.write('  #define ' + (funname + suffix).lower() + '_ ' + (funname + suffix).lower() + '\n')
    fd.write('#endif\n')

    # output function declaration prototype
    fd.write('SLEPC_EXTERN void ' + (funname + suffix).lower() + '_(')
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      ktypename = k.typename
      if cnt: fd.write(', ')
      if k.const and not ((k.typename == 'char' or k.typename in senums) and k.array):
        # see note one at the top of the file why const is not added for this case
        fd.write('const ')
      if k.typename in senums:
        fd.write('char *')
      else:
        fd.write(ktypename)
      fd.write(' ')
      if not (k.typename == 'char' or k.typename in senums or k.array or k.typename == 'PeCtx'):
        fd.write('*')
      fd.write(Letters[cnt])
      if k.array: fd.write('[]')
      cnt = cnt + 1
    if cnt: fd.write(', ')
    fd.write('PetscErrorCode *ierr')
    # add the lengths of the string arguments put in by the Fortran compiler
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.typename in senums or k.typename == 'char':
        fd.write(', PETSC_FORTRAN_CHARLEN_T l_'  + Letters[cnt])
      cnt = cnt + 1
    fd.write(')\n{\n')

    # functions that destroy objects should return immediately if null, -2, -3
    if fun.arguments and fun.arguments[0].typename in classes and fun.name.endswith('Destroy'):
      fd.write('  PETSC_FORTRAN_OBJECT_F_DESTROYED_TO_C_NULL(a);\n')

    # handle arguments that may return a null object
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.stars and k.typename  in classes:
        fd.write('  PetscBool null_' + Letters[cnt] + ' = !*(void**) ' + Letters[cnt] + ' ? PETSC_TRUE : PETSC_FALSE;\n')
      cnt = cnt + 1

    # prevent an existing object from being overwritten by a new create
    if fun.arguments and fun.arguments[-1].typename in classes and fun.name.endswith('Create'):
      fd.write('  PETSC_FORTRAN_OBJECT_CREATE(' + Letters[len(fun.arguments)-1] + ');\n')

    # handle string argument fixes
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.typename == 'char' or k.typename in senums:
        if not k.stars and (k.const or (k.typename in senums)):
          fd.write('  char* c_' + Letters[cnt] + ';\n')
          fd.write('  FIXCHAR(' + Letters[cnt] + ', l_' + Letters[cnt] + ', c_' + Letters[cnt] + ');\n')
        elif k.stars:
          fd.write('  char* c_' + Letters[cnt] + ' = PETSC_NULLPTR;\n')
      cnt = cnt + 1

    # handle viewer argument fixes
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.typename == 'PetscViewer' and not k.stars and not k.array:
        fd.write('  PetscViewer v_' + Letters[cnt] + ' = PetscPatchDefaultViewers(' + Letters[cnt] + ');\n')
      cnt = cnt + 1

    # handle any arguments that may be null
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.typename in classes and k.stars:
        fd.write('  CHKFORTRANNULLOBJECT(' + Letters[cnt] + ');\n')
      if k.typename == 'PetscInt' and (k.stars or k.array):
        fd.write('  CHKFORTRANNULLINTEGER(' + Letters[cnt] + ');\n')
      elif k.typename == 'PetscReal' and (k.stars or k.array):
        fd.write('  CHKFORTRANNULLREAL(' + Letters[cnt] + ');\n')
      elif k.typename == 'PetscScalar' and (k.stars or k.array):
        fd.write('  CHKFORTRANNULLSCALAR(' + Letters[cnt] + ');\n')
      elif k.typename == 'PetscBool' and (k.stars or k.array):
        fd.write('  CHKFORTRANNULLBOOL(' + Letters[cnt] + ');\n')
      cnt = cnt + 1

    # call C function
    fd.write('  *ierr = ' + funname + '(')
    cnt = 0
    for k in fun.arguments:
      if cnt: fd.write(', ')
      if k.typename in senums or k.typename == 'char':
        if k.stars:
          fd.write('(const char **)&')
        if k.const or (k.typename in senums):
          fd.write('c_')
      elif k.typename == 'MPI_Fint':
         fd.write('MPI_Comm_f2c(*(')
      elif not k.stars and not k.array and not k.stringlen and not k.typename == 'PetscViewer' and not k.typename == 'PeCtx':
        fd.write('*')
#      if k.typename == 'void' and k.stars == 2:
#        fd.write('&')
      if k.stringlen:
        fd.write('l_' + Letters[cnt - 1])
        continue
      if k.typename == 'PetscViewer' and not k.stars and not k.array:
        fd.write('v_')
      fd.write(Letters[cnt])
      if k.typename == 'PetscBool' and not k.stars and not k.array:
        # handle bool argument fixes (-1 needs to be corrected to 1 for Intel compilers)
        fd.write(' ? PETSC_TRUE : PETSC_FALSE')
      if k.typename == 'MPI_Fint':
        fd.write('))')
      cnt = cnt + 1
    fd.write(');\n')
    fd.write('  if (*ierr) return;\n');

    # cleanup any string arguments fixes
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.typename == 'char' or k.typename in senums:
        if not k.stars and (k.const or (k.typename in senums)):
          fd.write('  FREECHAR(' + Letters[cnt] + ', c_' + Letters[cnt] + ');\n')
        else:
          if k.stars:
            fd.write('  *ierr = PetscStrncpy((char *)' + Letters[cnt] + ', c_' + Letters[cnt] + ', l_' + Letters[cnt] + ');\n')
            fd.write('  if (*ierr) return;\n');
          fd.write('  FIXRETURNCHAR(PETSC_TRUE, ' + Letters[cnt] + ', l_' + Letters[cnt] + ');\n')
      cnt = cnt + 1

    # handle arguments that may return a null PETSc object
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if k.stars and k.typename in classes:
        fd.write('  if (! null_' + Letters[cnt] + ' && !*(void**) ' + Letters[cnt] + ') *(void **) ' + Letters[cnt] + ' = (void *)-2;\n')
      cnt = cnt + 1

    fd.write('}\n')
  shutil.copy(os.path.join(fun.dir,'makefile'), os.path.join(dir,'makefile'))

def generateFortranStub(senums, funname, fun, fd, opts):
  '''For functions with optional arguments generate the Fortran stub that calls the C stub'''
  for k in fun.arguments:
    # no C stub if function returns an array, except if it is a string
    # TODO: generate fillible stub for functions that return arrays
    if k.array and k.stars and not k.typename == 'char': return
    if k.stars and k.typename == 'MPI_Fint': return   # TODO add support for returning MPI_Fint
    if k.stars == 2 and k.typename == 'void': return
  for fi in opts:
    fd.write('  subroutine ' + funname + ''.join(fi) + '(')
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      if cnt: fd.write(',')
      fd.write(Letters[cnt])
      cnt = cnt + 1
    fd.write(',z)\n')
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      ktypename = k.typename
      if fi[cnt] == 'O':
        fd.write('  PetscNull :: ' + Letters[cnt] + '\n')
      elif ktypename in senums or ktypename == 'char':
        fd.write('  character(*) :: ' + Letters[cnt] + '\n')
      elif k.array and k.stars:
        fd.write('  ' + ktypename + ', pointer :: ' +  Letters[cnt] + '(:)\n')
      elif k.array:
        fd.write('  ' + ktypename + ' :: '+ Letters[cnt] + '(*)\n')
      else:
        fd.write('  '+ ktypename + ' :: ' + Letters[cnt] + '\n')
      cnt = cnt + 1
    fd.write('  PetscErrorCode z\n')
    fd.write('  call ' + funname + 'Raw(')
    cnt = 0
    for k in fun.arguments:
      if k.stringlen: continue
      ktypename = k.typename
      if cnt: fd.write(',')
      if fi[cnt] == 'a':
        fd.write(Letters[cnt])
      else:
        typename = ktypename.upper().replace('PETSC','')
        if typename == 'INT': typename = 'INTEGER'
        fd.write('PETSC_NULL_' + typename)
        if k.array and k.stars: fd.write('_POINTER')
        elif k.array: fd.write('_ARRAY')
      cnt = cnt + 1
    fd.write(',z)\n')
    fd.write('  end subroutine\n')

##########  main

def main(petscdir,slepcdir,petscarch):
  '''Generates all the Fortran include and C stub files needed for the Fortran API'''
  sys.path.insert(0,os.path.realpath(os.path.dirname(__file__)))
  import getAPI
  del sys.path[0]

  classes, enums, senums, typedefs, structs, funcs, files, mansecs, submansecs = getAPI.getAPI()

  # hardwired PetscObject functions
  petscobjectfunctions = {}

  argpetscobj = getAPI.Argument()
  argpetscobj.typename = 'PetscObject'
  argpetscobjstar = getAPI.Argument()
  argpetscobjstar.typename = 'PetscObject'
  argpetscobjstar.stars = 1
  argpetscobjid = getAPI.Argument()
  argpetscobjid.typename = 'PetscObjectId'
  argpetscobjidstar = getAPI.Argument()
  argpetscobjidstar.typename = 'PetscObjectId'
  argpetscobjidstar.stars = 1
  argpetscclassidstar = getAPI.Argument()
  argpetscclassidstar.typename = 'PetscClassId'
  argpetscclassidstar.stars = 1
  argpetscviewer = getAPI.Argument()
  argpetscviewer.typename = 'PetscViewer'
  argpetscoptions = getAPI.Argument()
  argpetscoptions.typename = 'PetscOptions'
  argpetscoptionsstar = getAPI.Argument()
  argpetscoptionsstar.typename = 'PetscOptions'
  argpetscoptionsstar.stars = 1
  argpetscoptitem = getAPI.Argument()
  argpetscoptitem.typename = 'PetscOptionItems'
  argstring = getAPI.Argument()
  argstring.typename = 'char'
  argstring.array = True
  argstringstar = getAPI.Argument()
  argstringstar.typename = 'char'
  argstringstar.array = True
  argstringstar.stars = 1
  argint = getAPI.Argument()
  argint.typename = 'PetscInt'
  argintstar = getAPI.Argument()
  argintstar.typename = 'PetscInt'
  argintstar.stars = 1
  argmpiintstar = getAPI.Argument()
  argmpiintstar.typename = 'PetscMPIInt'
  argmpiintstar.stars = 1
  argboolstar = getAPI.Argument()
  argboolstar.typename = 'PetscBool'
  argboolstar.stars = 1
  argcommstar = getAPI.Argument()
  argcommstar.typename = 'MPI_Comm'
  argcommstar.stars = 1

  fname = 'PetscObjectGetType'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstringstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetPrintedOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectInheritPrintedOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectProcessOptionsHandlers'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscoptitem]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectDestroyOptionsHandlers'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectReference'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetReference'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj,argintstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectDereference'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectRemoveReference'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectHasFunction'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring, argboolstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetFromOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetUp'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetName'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstringstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectDestroy'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobjstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectView'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscviewer]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectViewFromOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectTypeCompare'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring, argboolstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectObjectTypeCompare'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobj, argboolstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectBaseTypeCompare'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring, argboolstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectRegisterDestroy'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetId'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobjidstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectCompareId'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobjid, argboolstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscoptionsstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetOptions'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscoptions]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetOptionsPrefix'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectAppendOptionsPrefix'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetOptionsPrefix'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstringstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectPrependOptionsPrefix'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetComm'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argcommstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetTabLevel'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argintstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetTabLevel'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argint]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectIncrementTabLevel'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscobj, argint]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetNewTag'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argmpiintstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetClassId'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscclassidstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectGetClassName'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstringstar]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectSetName'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectPrintClassNamePrefixType'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argpetscviewer]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectName'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj]
  petscobjectfunctions[fname] = fun

  fname = 'PetscObjectChangeTypeName'
  fun = getAPI.Function(fname)
  fun.arguments = [argpetscobj, argstring]
  petscobjectfunctions[fname] = fun

##########  $PETSC_ARCH/include/petsc/finclude/*.h

  dir = os.path.join(petscarch,'include', 'slepc', 'finclude')
  if os.path.isdir(dir): shutil.rmtree(dir)
  os.makedirs(dir)

  for i in files.keys():
    if i.endswith('types.h') or i in skipinc: continue
    with open(os.path.join(dir, i),'w') as fd:
      dname = 'SLEPC' + i.upper()[0:-2] + 'DEF_H'
      fd.write('#if !defined(' + dname + ')\n#define ' + dname + '\n\n')
      fb = os.path.join('include', 'slepc', 'finclude',i.replace('.h','base.h'))
      if os.path.isfile(fb):
        fd.write('#include "' + os.path.join('slepc','finclude',i.replace('.h','base.h')) + '"\n')
      for j in files[i].included:
        if j in skipinc: continue
        j = j.replace('types.h','.h')
        if i == j: continue
        fd.write('#include "' + os.path.join(('petsc' if j.startswith('petsc') else 'slepc'),'finclude',j) + '"\n')
      fd.write('\n')

  for i in enums.keys():
    with open(os.path.join(dir, enums[i].includefile),"a") as fd:
      fd.write('#define ' + i + ' type(e' + i + ')\n')

  for i in typedefs.keys():
    if not typedefs[i].name: continue
    value = typedefs[i].value
    if value in CToFortranTypes:
      if not CToFortranTypes[value]: continue
      value = CToFortranTypes[value]
    with open(os.path.join(dir, typedefs[i].includefile),"a") as fd:
      fd.write('#define ' + i + ' ' + value + '\n')

  for i in structs.keys():
    with open(os.path.join(dir, structs[i].includefile),"a") as fd:
      if not structs[i].opaque:
        fd.write('#define ' +  i  + ' type(s' + i + ')\n')
      else:
        fd.write('#define ' + i + ' PetscFortranAddr\n')

  for i in files.keys():
    if i.endswith('types.h') or i in skipinc: continue
    with open(os.path.join(dir, i),'a') as fd:
      fd.write('\n')

  for i in senums.keys():
    with open(os.path.join(dir, senums[i].includefile),"a") as fd:
      fd.write('#define ' + i + ' CHARACTER(80)\n')

  for i in files.keys():
    if i.endswith('types.h') or i in skipinc: continue
    with open(os.path.join(dir, i),'a') as fd:
      fd.write('\n')

  for i in classes.keys():
    with open(os.path.join(dir, classes[i].includefile),"a") as fd:
      fd.write('#define ' + i + ' type(t' + i + ')\n')

  for i in files.keys():
    if i.endswith('types.h') or i in skipinc: continue
    with open(os.path.join(dir, i),'a') as fd:
      fd.write('\n#endif\n')

###########  $PETSC_ARCH/ftn/MANSEC/*.h

  dir = os.path.join(petscarch,'ftn')
  if os.path.isdir(dir): shutil.rmtree(dir)
  os.makedirs(dir)

  for i in mansecs.keys():
    if i!='sys':
      mansecpath = (os.path.join('sys','classes',i) if i in ['bv','ds','fn','rg','st'] else i)
      dir = os.path.join(petscarch,'ftn', mansecpath)
      if os.path.isdir(dir): shutil.rmtree(dir)
      os.makedirs(dir)

  for i in classes.keys():
    mansecpath = (os.path.join('sys','classes',classes[i].mansec) if classes[i].mansec in ['bv','ds','fn','rg','st'] else classes[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath, classes[i].includefile),"a") as fd:
      if not classes[i].petscobject:
        fd.write('  type t' + i + '\n')
        fd.write('    PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE\n')
        fd.write('  end type t' + i + '\n')
      else:
        fd.write('  type, extends(tPetscObject) ::  t' + i + '\n')
        fd.write('  end type t' + i + '\n')
      v = ('SLEPC_NULL_' + i.upper().replace('SLEPC','').strip('_')).replace('_NULL_NULL','_NULL')
      fd.write('  ' + i + ', parameter :: ' + v + ' = t' + i + '(0)\n')
      fd.write('  ' + i + ', target :: ' + v + '_ARRAY(1) = [t' + i + '(0)]\n')
      fd.write('  ' + i + ', pointer :: ' + v + '_POINTER(:) => ' + v + '_ARRAY\n')
      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + v + '\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + v + '_ARRAY\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + v + '_POINTER\n')
      fd.write('#endif\n')
      fd.write('\n')

  for i in enums.keys():
    mansecpath = (os.path.join('sys','classes',enums[i].mansec) if enums[i].mansec in ['bv','ds','fn','rg','st'] else enums[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,enums[i].includefile),"a") as fd:
      fd.write('  type e' + i + '\n')
      fd.write('    PetscEnum:: v PETSC_FORTRAN_TYPE_INITIALIZE\n')
      fd.write('  end type e' + i + '\n\n')
      v = ('SLEPC_NULL_' + i.upper().replace('SLEPC','').replace('NULL','')).strip('_').replace('_NULL_NULL','_NULL')
      fd.write('  ' + i + ', parameter :: ' + v + ' = e' + i + '(-50)\n')
      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + v + '\n')
      fd.write('#endif\n')
      cnt = 0
      givenvalue = 0
      for j in enums[i].values:
        if j.find('=') > -1:
          if givenvalue == -1:
            print('Some enum values for ' + i + ' are set but others are not set')
          v = j.replace(' = ',' = e' + i + '(') + ')'
          givenvalue = 1
        else:
          if givenvalue == 1:
            print('Some enum values for ' + i + ' are set but others are not set')
          v = j + ' = e' + i + '(' + str(cnt) + ')'
          givenvalue = -1
        fd.write('    ' + i + ', parameter :: ' + v + '\n')
        cnt = cnt + 1
      fd.write('\n')

      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      for j in enums[i].values:
         if j.count('='): v = j[0:j.find('=')]
         else: v = j
         fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + v + '\n')
      fd.write('#endif\n')
      fd.write('\n')

  for i in senums.keys():
    mansecpath = (os.path.join('sys','classes',senums[i].mansec) if senums[i].mansec in ['bv','ds','fn','rg','st'] else senums[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,senums[i].includefile),"a") as fd:
      for j in senums[i].values:
        fd.write('  CHARACTER(LEN=*), PARAMETER :: ' + j + ' = \'' + senums[i].values[j].replace('"','') + '\'\n')
      fd.write('\n')

      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      for j in senums[i].values:
        fd.write('!DEC$ ATTRIBUTES DLLEXPORT::' + j + '\n')
      fd.write('#endif\n')
      fd.write('\n')

  for i in structs.keys():
    if structs[i].opaque: continue
    mansecpath = (os.path.join('sys','classes',structs[i].mansec) if structs[i].mansec in ['bv','ds','fn','rg','st'] else structs[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,structs[i].includefile),"a") as fd:
      fd.write('  type s' + i + '\n')
      for j in structs[i].records:
        fd.write('    ' + j.type.replace('[','(').replace(']',')') + '\n')
      fd.write('  end type s' + i + '\n')
      fd.write('\n')

###########  $PETSC_ARCH/ftn/MANSEC/*.h90

  for i in classes.keys():
    # generate interface definitions for all objects' methods
    for j in classes[i].functions: # loop over functions in class
      generateFortranInterface(petscarch,classes,enums,structs,senums,j,classes[i].functions[j])

    if i in ['SlepcConvMon']: continue
    file = classes[i].includefile + '90'
    if not file.startswith('slepc'): file = 'slepc' + file
    mansecpath = (os.path.join('sys','classes',classes[i].mansec) if classes[i].mansec in ['bv','ds','fn','rg','st'] else classes[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,file),"a") as fd:
      fd.write('  interface operator(.ne.)\n')
      fd.write('    module procedure ' + i + 'notequals\n')
      fd.write('  end interface operator (.ne.)\n')
      fd.write('  interface operator(.eq.)\n')
      fd.write('    module procedure ' + i + 'equals\n')
      fd.write('  end interface operator (.eq.)\n\n')

    # generate interface definitions for PetscObject methods for each PetscObject subclass (KSP etc)
    if not classes[i].petscobject: continue
    mansecpath = (os.path.join('sys','classes',classes[i].mansec) if classes[i].mansec in ['bv','ds','fn','rg','st'] else classes[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,file),"a") as fd:
      ii = i.replace('Petsc','')
      fd.write('  interface PetscObjectCast\n')
      fd.write('    module procedure PetscObjectCast' + ii + '\n')
      fd.write('  end interface\n')
      fd.write('  interface PetscBarrier\n')
      fd.write('  module procedure PetscBarrier' + ii + '\n')
      fd.write('  end interface\n')

      for funname in petscobjectfunctions:
        fi = petscobjectfunctions[funname]
        fd.write('  interface ' + funname + '\n')
        fd.write('    module procedure ' + funname  + ii + '\n')
        fd.write('  end interface\n')

  # generate interface definitions for all standalone functions
  for j in funcs.keys():
    generateFortranInterface(petscarch,classes,enums,structs,senums,funcs[j].name,funcs[j])

  # generate .eq. and .neq. for enums
  for i in enums.keys():
    file = enums[i].includefile + '90'
    if not file.startswith('slepc'): file = 'slepc' + file
    mansecpath = (os.path.join('sys','classes',enums[i].mansec) if enums[i].mansec in ['bv','ds','fn','rg','st'] else enums[i].mansec)
    with open(os.path.join(petscarch,'ftn',mansecpath,file),"a") as fd:
      fd.write('  interface operator(.ne.)\n')
      fd.write('    module procedure ' + i + 'notequals\n')
      fd.write('  end interface operator (.ne.)\n')
      fd.write('  interface operator(.eq.)\n')
      fd.write('    module procedure ' + i + 'equals\n')
      fd.write('  end interface operator (.eq.)\n\n')

##########  $PETSC_ARCH/ftn/MANSEC/*.hf90

  for i in classes.keys():
    if i in ['SlepcConvMon']: continue
    mansecpath = (os.path.join('sys','classes',classes[i].mansec) if classes[i].mansec in ['bv','ds','fn','rg','st'] else classes[i].mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,classes[i].includefile + 'f90'),"a") as fd:

      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT:: ' + i + 'notequals\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT:: ' + i + 'equals\n')
      fd.write('#endif\n\n')
      fd.write('  function ' + i + 'notequals(A,B)\n')
      fd.write('    logical ' + i + 'notequals\n')
      fd.write('    type(t' + i + '), intent(in) :: A,B\n')
      fd.write('    ' + i + 'notequals = (A%v .ne. B%v)\n')
      fd.write('  end function\n')
      fd.write('  function ' + i + 'equals(A,B)\n')
      fd.write('    logical ' + i + 'equals\n')
      fd.write('    type(t' + i + '), intent(in) :: A,B\n')
      fd.write('    ' + i + 'equals = (A%v .eq. B%v)\n')
      fd.write('  end function\n')

       # generate Fortran subroutines for PetscObject methods for each PetscObject subclass (KSP etc)
      if classes[i].petscobject:
        ii = i.replace('Petsc', '')
        fd.write('  function PetscObjectCast' + ii + '(a)\n')
        fd.write('    ' + i + ' a\n')
        fd.write('    PetscObject PetscObjectCast' + ii + '\n')
        fd.write('    PetscObjectCast' + ii + '%v = a%v\n')
        fd.write('  end function \n')
        fd.write('  subroutine PetscBarrier' + ii + '(a,z)\n')
        fd.write('    ' + i + ' a\n')
        fd.write('    PetscErrorCode z\n')
        fd.write('    call PetscBarrier(PetscObjectCast(a),z)\n')
        fd.write('  end subroutine \n')

        for funname in petscobjectfunctions:
          fi = petscobjectfunctions[funname]
          fd.write('  subroutine ' + funname + ii + '(')
          cnt = 0
          for k in fi.arguments:
            if k.stringlen: continue
            if cnt: fd.write(', ')
            fd.write(Letters[cnt])
            cnt = cnt + 1
          fd.write(' , z)\n')
          fd.write('  ' + i + ' a\n')
          cnt = 1
          for k in fi.arguments[1:]:
            if k.stringlen: continue
            ktypename = k.typename
            if ktypename in CToFortranTypes:
              ktypename = CToFortranTypes[ktypename]
            if ktypename in senums or ktypename == 'char':
              fd.write('  character(*) :: ' + Letters[cnt] + '\n')
            elif k.array and k.stars:
              fd.write('  ' + ktypename + ', pointer :: ' +  Letters[cnt]  + '(:)\n')
            elif k.array:
              fd.write('  ' + ktypename + ' :: ' +  Letters[cnt]  + '(*)\n')
            else:
              fd.write('  ' + ktypename + ' :: ' + Letters[cnt] + '\n')
            cnt = cnt + 1
          fd.write('  PetscErrorCode z\n')
          fd.write('  call ' + funname  + '(PetscObjectCast(')
          cnt = 0
          for k in fi.arguments:
            if k.stringlen: continue
            if cnt: fd.write(', ')
            fd.write(Letters[cnt])
            if cnt == 0: fd.write(')')
            cnt = cnt + 1
          fd.write(', z)\n')
          fd.write('  end subroutine \n')

  # generate all the polymorphic Fortran subroutines needed for class methods with optional arguments
  # not used because it takes too long to build the Fortran modules
  for i in []: #classes.keys():
    if i in ['PetscIntStack']: continue
    for j in classes[i].functions: # loop over functions in class
      # check for functions for which we cannot build interfaces
      if classes[i].functions[j].opaque: continue
      opts = crossCreate(classes[i].functions[j])
      if len(opts) == 1: continue
      mansec = classes[i].mansec
      file = classes[i].functions[j].includefile + 'f90'
      if not file.startswith('petsc'): file = 'slepc' + file
      mansecpath = (os.path.join('sys','classes',mansec) if mansec in ['bv','ds','fn','rg','st'] else mansec)
      with open(os.path.join(petscarch,'ftn', mansecpath,file),'a') as fd:
        generateFortranStub(senums, j, classes[i].functions[j], fd, opts)

  # generate all the polymorphic Fortran subroutines needed for class-less routines with optional arguments
  # not used because it takes too long to build the Fortran modules
  for j in []: # funcs.keys():
    if funcs[j].opaque: continue
    opts = crossCreate(funcs[j])
    if len(opts) == 1: continue
    mansec = funcs[j].mansec
    file = funcs[j].includefile + 'f90'
    if not file.startswith('slepc'): file = 'slepc' + file
    mansecpath = (os.path.join('sys','classes',mansec) if mansec in ['bv','ds','fn','rg','st'] else mansec)
    with open(os.path.join(petscarch,'ftn', mansecpath,file),'a') as fd:
      generateFortranStub(senums,funcs[j].name,funcs[j], fd, opts)

  # generate .eq. and .neq. for enums
  for i in enums.keys():
    mansecpath = (os.path.join('sys','classes',enums[i].mansec) if enums[i].mansec in ['bv','ds','fn','rg','st'] else enums[i].mansec)
    with open(os.path.join(petscarch,'ftn',mansecpath,enums[i].includefile + 'f90'),"a") as fd:

      fd.write('#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT:: ' + i + 'notequals\n')
      fd.write('!DEC$ ATTRIBUTES DLLEXPORT:: ' + i + 'equals\n')
      fd.write('#endif\n\n')
      fd.write('  function ' + i + 'notequals(A,B)\n')
      fd.write('    logical ' + i + 'notequals\n')
      fd.write('    type(e' + i + '), intent(in) :: A,B\n')
      fd.write('    ' + i + 'notequals = (A%v .ne. B%v)\n')
      fd.write('  end function\n')
      fd.write('  function ' + i + 'equals(A,B)\n')
      fd.write('    logical ' + i + 'equals\n')
      fd.write('    type(e' + i + '), intent(in) :: A,B\n')
      fd.write('    ' + i + 'equals = (A%v .eq. B%v)\n')
      fd.write('  end function\n')

##########  $PETSC_ARCH/ftn/MANSEC/**/*f.c

  # convert function arguments from MPI_Comm to MPI_Fint
  for i in funcs:
    for j in funcs[i].arguments:
      j.typename = j.typename.replace('MPI_Comm','MPI_Fint')

  for i in classes:
    for j in classes[i].functions:
      for k in classes[i].functions[j].arguments:
        k.typename = k.typename.replace('MPI_Comm','MPI_Fint')

  for i in classes.keys():
    if i in ['PetscIntStack']: continue
    for j in classes[i].functions: # loop over functions in class
      generateCStub(petscarch,senums,classes,j,classes[i].functions[j])

  for j in funcs.keys():
    if funcs[j].name in ['SlepcDebugViewMatrix']: continue
    generateCStub(petscarch,senums,classes,funcs[j].name,funcs[j])

##########  $PETSC_ARCH/ftn/MANSEC/petscall.*

  # slepcall.* contains all the include files associated with C slepcMANSEC.h
  # these are used by src/MANSEC/ftn-mod/slepcMANSECmod.F to generate the module for C slepcMANSEC.h
  # src/MANSEC/ftn-mod/slepcMANSECmod.F may also define additional modules that use slepcMANSEC
  for i in mansecs.keys():
    mansecpath = (os.path.join('sys','classes',i) if i in ['bv','ds','fn','rg','st'] else i)
    d = os.path.join(petscarch,'ftn', mansecpath)
    dd = os.path.join('../','ftn', mansecpath)
    args = [os.path.join(d,i) for i in os.listdir(d) if i.endswith('.h')]
    for j in args:
      if not os.path.getsize(j): os.path.remove(j)
    with open(os.path.join(d,'slepcall.h'),'w') as fd, open(os.path.join(d,'slepcall.h90'),'w') as fd90, open(os.path.join(d,'slepcall.hf90'),'w') as fdf90:
      if not i.startswith('slepc'): f = 'slepc' + i + '.h'
      else: f = i + '.h'
      includes = set()
      for j in files[f].included:
        if j in skipinc: continue
        j = j.replace('types.h','.h')
        includes.add(j)
        fd.write('#include <' + os.path.join(('petsc' if j.startswith('petsc') else 'slepc'),'finclude',j) + '>\n')
        if os.path.isfile(os.path.join(d,j)) :
          fd.write('#include <' + os.path.join(dd,j) + '>\n')
        if os.path.isfile(os.path.join(d,j + '90')):
          fd90.write('#include <' + os.path.join(dd,j + '90') + '>\n')
        if os.path.isfile(os.path.join(d,j + 'f90')):
          fdf90.write('#include <' + os.path.join(dd,j + 'f90') + '>\n')
      if f not in includes:
        fd.write('#include <' + os.path.join('slepc','finclude',f) + '>\n')
        fd.write('#include <' + os.path.join(dd,f) + '>\n')
        if os.path.isfile(os.path.join(d,f + '90')):
          fd90.write('#include <' + os.path.join(dd,f + '90') + '>\n')
        if os.path.isfile(os.path.join(d,f + 'f90')):
          fdf90.write('#include <' + os.path.join(dd,f + 'f90') + '>\n')

#
if __name__ ==  '__main__':
  import sys
  import argparse

  parser = argparse.ArgumentParser(description='generate SLEPc FORTRAN stubs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--petsc-dir', metavar='path', required=True, help='PETSc root directory')
  parser.add_argument('--slepc-dir', metavar='path', required=True, help='SLEPc root directory')
  parser.add_argument('--petsc-arch', metavar='string', required=True, help='PETSc arch name')
  args = parser.parse_args()

  ret = main(args.petsc_dir, args.slepc_dir, args.petsc_arch)
  sys.exit(ret)
