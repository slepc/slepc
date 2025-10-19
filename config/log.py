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
import os, sys

class Log:

  def __init__(self):
    self.fd = None
    self.laststatus = 'done'

  def Open(self,slepcdir,confdir,fname):
    filename = os.path.join(confdir,fname)
    self.fd = open(filename,'w')
    self.filename = os.path.relpath(filename,slepcdir)
    try: # symbolic link to log file in current directory
      if os.path.isfile(fname) or os.path.islink(fname): os.remove(fname)
      os.symlink(self.filename,fname)
    except: pass

  def Println(self,string):
    print(string)
    if self.fd:
      self.fd.write(string+'\n')

  def Print(self,string):
    print(string, end=' ')
    if self.fd:
      self.fd.write(string+' ')

  def NewSection(self,string):
    if not self.laststatus == 'done':
      colorfail = '\033[91m'
      colornorm = '\033[0m'
      print(colorfail+self.laststatus+colornorm+'\n'+string, end=' ')
    else:
      print(self.laststatus+'\n'+string, end=' ')
    sys.stdout.flush()
    if self.fd:
      self.fd.write('='*80+'\n'+string+'\n')
    self.laststatus = 'done'

  def write(self,string):
    if self.fd:
      self.fd.write(string+'\n')

  def Warn(self,string):
    msg = '\nxxx'+'='*74+'xxx\nWARNING: '+string+'\nxxx'+'='*74+'xxx'
    print(msg)
    if self.fd:
      self.fd.write(msg+'\n')

  def Exit(self,string):
    msg = '\nERROR: '+string
    print(msg)
    if self.fd:
      self.fd.write(msg+'\n')
      self.fd.close()
      msg = 'ERROR: See "' + self.filename + '" file for details'
    else:
      msg = 'ERROR during configure (log file not open yet)'
    sys.exit(msg)

  def setLastStatus(self, stat):
    if stat not in ['done', 'failed', 'skipped']:
      self.Exit('Unknown value of argument stat='+stat)
    self.laststatus = stat

  def Close(self):
    if self.fd:
      self.fd.close()
