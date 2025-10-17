#!/usr/bin/env python3
""" Runs c2html and mapnames on a single source file and cleans up all text"""

from typing import Dict
import os
import re
import subprocess
import pathlib
import sys
import posixpath

from sphinx.util.inventory import InventoryFile

word_pattern = re.compile(r'\w+')

def dict_complete_links(string_to_link: Dict[str,str], prefix: str = '') -> Dict[str,str]:
  """
  Prepend a prefix to any links not starting with 'http' so Sphinx will recognize them as URLs
  """
  def link_string(name: str, link: str, prefix: str) -> str:
    url = link if link.startswith('http') else prefix + link
    return '<a href=\"' + url + '\">' + name + '</a>'
  return dict((k, link_string(k, v, prefix)) for (k, v) in string_to_link.items())

def _get_inventory() -> Dict[str,str]:
  inventory_prefix = 'https://petsc.org/release/'
  inventory_raw = {}
  inventory_filename = os.path.join('doc','source','petsc_objects.inv')
  with open(inventory_filename, 'rb') as f:
    try:
      _inventory = InventoryFile.load(f, inventory_prefix, posixpath.join)
    except ValueError as exc:
      raise ValueError('unknown or unsupported inventory version: %r' % exc) from exc
    stddoc = _inventory['std:doc']
    for k, v in stddoc.items():
      value = stddoc[k]
      path = value[2]
      inventory_raw[k] = path
  inventory = dict_complete_links(inventory_raw, inventory_prefix)
  return inventory

def main(slepc_dir,loc,git_sha,c2html,mapnames,rel_dir,file):
  with open(os.path.join(rel_dir,file), "r") as fd:
    txt = fd.read()

  # TODO change text processing parts to Python
  cmd = 'sed -E "s/PETSC[A-Z]*_DLLEXPORT//g"  | '+ \
         c2html + ' -n | \
         awk \'{ sub(/<pre width="80">/,"<pre width=\"80\">\\n"); print }\' | \
         grep -v "#if !defined(__" | grep -E -v "(PetscValid|#define __|#undef __|EXTERN_C )" | ' + \
         mapnames + ' -map htmlmap.tmp -inhtml'
  txt = subprocess.check_output(cmd, text=True, input=txt, shell = True)

  # make the links to manual pages relative
  rel_dot = '../'
  for c in rel_dir:
    if c == '/':
      rel_dot = rel_dot + '../'
  txt = txt.replace('HTML_ROOT/',rel_dot)

  # make the links to include files relative
  ntxt = ''
  for line in txt.split('\n'):
    if 'include' in line:
      ins = re.search('#include [ ]*&lt;',line)
      if ins:
        includename = line[ins.end():re.search('&gt;[a-zA-Z0-9/<>#*"=. ]*',line).start()]
        ln = re.search('<a name="line[0-9]*">[ 0-9]*: </a>',line)
        linenumber = line[ln.start():ln.end()]
        if os.path.isfile(includename):
          line = linenumber+'#include <A href="'+includename+'.html">&lt;'+includename+'&gt;</A>'
        elif os.path.isfile(os.path.join('include',includename)):
          line = linenumber+'#include <A href="'+os.path.relpath(os.path.join(rel_dot,'include',includename))+'.html">&lt;'+includename+'&gt;</A>'
        elif os.path.isfile(os.path.join(includename)):
          line = linenumber+'#include <A href="'+os.path.relpath(os.path.join(rel_dot,includename))+'.html">&lt;'+includename+'&gt;</A>'
    ntxt = ntxt + line + '\n'

  # use petsc_objects.inv to create links to PETSc's documentation
  inventory = _get_inventory()
  def replace(matchobj):
    word = matchobj.group(0)
    if word in inventory:
      return inventory[word]
    return word
  txt = ''
  for line in ntxt.split('\n'):
    txt = txt + word_pattern.sub(replace,line) + '\n'
  txt = txt + '\n'

  with open(os.path.join(loc,rel_dir,file+'.html'), "w") as fdw:
    fdw.write('<center><a href="https://gitlab.com/slepc/slepc/-/blob/'+git_sha+'/'+rel_dir+'/'+file+'">Actual source code: '+file+'</a></center><br>\n')
    fdw.write(txt)

if __name__ == "__main__":
  main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],os.path.dirname(sys.argv[6]),os.path.basename(sys.argv[6]))
