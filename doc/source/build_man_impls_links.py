#!/usr/bin/env python3
""" Adds links in the manual pages to implementations of the function
    Also adds References section if any {cite} are found in the manual page
"""

import os
import re

def processfile(slepc_dir,srcdir,dir,file,implsClassAll,subimplsClassAll,implsFuncAll):
  #print('Processing '+os.path.join(dir,file))
  isclass = False
  with open(os.path.join(dir,file),'r') as f:
    text = f.read()
    isclass = text.find('typedef struct') > -1
  if text.find('{cite}') > -1 or text.find('{cite:p}') > -1 or text.find('{cite:t}') > -1:
    with open(os.path.join(dir,file),'w') as f:
      f.write(text[0:text.find('## See Also')])
      f.write('\n## References\n```{bibliography}\n:filter: docname in docnames\n```\n\n')
      f.write(text[text.find('## See Also'):])

  itemName = file[0:-3]
  if isclass:
    iclass = list(filter(lambda x: x.find('_p_'+itemName+' ') > -1, implsClassAll))
    func = None
    isubclass = list(filter(lambda x: x.find('} '+itemName+'_') > -1, subimplsClassAll))
  else:
    iclass = None
    isubclass = None
    func = list(filter(lambda x: x.find(' '+itemName+'_') > -1, implsFuncAll))
  if func or iclass:
    with open(os.path.join(dir,file),'a') as f:
      f.write('\n## Implementations\n')
      if func:
        for str in func:
          f.write(re.sub(r'(.*\.[ch]x*u*).*('+itemName+r'.*)(\(.*\))','<A HREF=\"PETSC_DOC_OUT_ROOT_PLACEHOLDER/\\1.html#\\2\">\\2() in \\1</A><BR>',str,count=1)+'\n')
      if iclass:
        for str in iclass:
          f.write(re.sub(r'(.*\.[ch]x*u*):.*struct.*(_p_'+itemName+').*{','<A HREF=\"PETSC_DOC_OUT_ROOT_PLACEHOLDER/\\1.html#\\2\">\\2 in \\1</A><BR>',str,count=1)+'\n')
      if isubclass:
        for str in isubclass:
          f.write(re.sub(r'(.*\.[ch]x*u*):} ('+itemName+'_.*);','<A HREF=\"PETSC_DOC_OUT_ROOT_PLACEHOLDER/\\1.html#\\2\">\\2 in \\1</A><BR>',str,count=1)+'\n')

def loadstructfunctions(slepc_dir):
  '''Creates the list of structs and class functions'''
  import subprocess
  implsClassAll = subprocess.check_output(['git', 'grep', '-E', r'struct[[:space:]]+_[pn]_[^[:space:]]+.*\{', '--', '*.c', '*.cpp', '*.cu', '*.c', '*.h', '*.cxx'], cwd = slepc_dir).strip().decode('utf-8')
  implsClassAll = list(filter(lambda x: not (x.find('/tests/') > -1 or x.find('/tutorials') > -1 or x.find(';') > -1), implsClassAll.split('\n')))

  subimplsClassAll = subprocess.check_output(['git', 'grep', '-E', '} [A-Z][A-Za-z]*_[A-Za-z]*;', '--', '*.c', '*.cpp', '*.cu', '*.c', '*.h', '*.cxx'], cwd = slepc_dir).strip().decode('utf-8')
  subimplsClassAll = list(filter(lambda x: not (x.find('/tests/') > -1 or x.find('/tutorials') > -1), subimplsClassAll.split('\n')))

  implsFuncAll = subprocess.check_output(['git', 'grep', '-nE', r'^(static )?(PETSC_EXTERN )?(PETSC_INTERN )?(extern )?PetscErrorCode +[^_ ]+_([^_ ]+|SuperLU_DIST|MKL_[C]{0,1}PARDISO)\(', '--', '*/impls/*.c', '*/impls/*.cpp', '*/impls/*.cu', '*/impls/*.c', '*/impls/*.h', '*/impls/*.cxx'], cwd = slepc_dir).strip().decode('utf-8')
  implsFuncAll = list(filter(lambda x: not (x.find('_Private') > -1 or x.find('_private') > -1 or x.find(';') > -1), implsFuncAll.split('\n')))
  return (implsClassAll,subimplsClassAll,implsFuncAll)

def main(slepc_dir,srcdir):
    (implsClassAll,subimplsClassAll,implsFuncAll) = loadstructfunctions(slepc_dir)
    for dirpath, dirnames, filenames in os.walk(os.path.join(srcdir,'manualpages'),topdown=True):
      #print('Processing directory '+dirpath)
      for file in filenames:
        if file.endswith('.md'): processfile(slepc_dir,srcdir,dirpath,file,implsClassAll,subimplsClassAll,implsFuncAll)
