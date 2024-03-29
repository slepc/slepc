#!/usr/bin/env python3

""" Reads in all the generated manual pages, and creates the index
for the manualpages, ordering the indices into sections based
on the 'Level of Difficulty'.

 Usage:
   wwwindex.py PETSC_DIR LOC
"""

import os
import sys
import re
import glob
import posixpath
import subprocess

# This routine reorders the entries int he list in such a way, so that
# When they are printed in a row order, the entries are sorted by columns
# In this subroutine row,col,nrow,ncol correspond to the new layout
# of a given data-value
def maketranspose(data,ncol):
      nrow = (len(data)+ncol-1)//ncol
      newdata = []
      # use the complete nrow by ncol matrix
      for i in range(nrow*ncol):
            newdata.append('')
      for i in range(len(data)):
            col           = i//nrow
            row           = i%nrow
            newi          = row*ncol+col
            newdata[newi] = data[i]
      return newdata

# Now use the level info, and print a html formatted index
# table. Can also provide a header file, whose contents are
# first copied over.
def printindex(outfilename, headfilename, levels, titles, tables):
      # Read in the header file
      headbuf = ''
      if posixpath.exists(headfilename) :
            with open(headfilename, "r") as fd:
                headbuf = fd.read()
                headbuf = headbuf.replace('PETSC_DIR', '../../../')
      else:
            print('Header file \'' + headfilename + '\' does not exist')

      # Now open the output file.
      try:
            fd = open(outfilename,'w')
      except:
            print('Error writing to file',outfilename)
            exit()

      # Add the HTML Header info here.
      fd.write(headbuf)
      # Add some HTML separators
      fd.write('\n<P>\n')
      fd.write('<TABLE>\n')
      for i in range(len(levels)):
            level = levels[i]
            title = titles[i]

            if len(tables[i]) == 0:
                  # If no functions in 'None' category, then don't print
                  # this category.
                  if level == 'none':
                        continue
                  else:
                        # If no functions in any other category, then print
                        # the header saying no functions in this cagetory.
                        fd.write('<TR><TD WIDTH=250 COLSPAN="3">')
                        fd.write('<B>' + 'No ' + level +' routines' + '</B>')
                        fd.write('</TD></TR>\n')
                        continue

            fd.write('<TR><TD WIDTH=250 COLSPAN="3">')
            #fd.write('<B>' + upper(title[0])+title[1:] + '</B>')
            fd.write('<B>' + title + '</B>')
            fd.write('</TD></TR>\n')
            # Now make the entries in the table column oriented
            tables[i] = maketranspose(tables[i],3)
            for filename in tables[i]:
                  path,name     = posixpath.split(filename)
                  func_name,ext = posixpath.splitext(name)
                  mesg          = ' <TD WIDTH=250><A HREF="'+ './' + name + '">' + \
                                  func_name + '</A></TD>\n'
                  fd.write(mesg)
                  if tables[i].index(filename) % 3 == 2 : fd.write('<TR>\n')
      fd.write('</TABLE>\n')
      # Add HTML tail info here
      fd.write('<BR><A HREF="../../../docs/manual.html">Table of Contents</A>\n')
      fd.close()

# This routine takes in as input a dictionary, which contains the
# alhabetical index to all the man page functions, and prints them all in
# a single index page
def printsingleindex(outfilename,alphabet_dict):
      # Now open the output file.
      try:
            fd = open(outfilename,'w')
      except:
            print('Error writing to file',outfilename)
            exit()

      alphabet_index = list(alphabet_dict.keys())
      alphabet_index.sort()

      # Now print each section, beginning with a title
      for key in alphabet_index:

            # Print the HTML tag for this section
            fd.write('<A NAME="' + key + '"></A>\n' )

            # Print the HTML index at the beginning of each section
            fd.write('<H3> <CENTER> | ')
            for key_tmp in alphabet_index:
                  if key == key_tmp:
                        fd.write( '<FONT COLOR="#883300">' + key_tmp.upper() + '</FONT> | \n' )
                  else:
                        fd.write('<A HREF="singleindex.html#' + key_tmp + '"> ' + \
                                 key_tmp.upper() + ' </A> | \n')
            fd.write('</CENTER></H3> \n')

            # Now write the table entries
            fd.write('<TABLE>\n')
            fd.write('<TR><TD WIDTH=250 COLSPAN="3">')
            fd.write('</TD></TR>\n')
            function_dict  = alphabet_dict[key]
            function_index = list(function_dict.keys())
            function_index.sort()
            function_index = maketranspose(function_index,3)
            for name in function_index:
                  if name:
                        path_name = function_dict[name]
                  else:
                        path_name = ''
                  mesg = '<TD WIDTH=250><A HREF="'+ './' + path_name + '">' + \
                         name + '</A></TD>\n'
                  fd.write(mesg)
                  if function_index.index(name) %3 == 2: fd.write('<TR>\n')

            fd.write('</TABLE>')

      fd.close()
      return


# Read in the filename contents, and search for the formatted
# String 'Level:' and return the level info.
# Also adds the BOLD HTML format to Level field
def modifylevel(filename,secname):
      with open(filename, "r") as fd:
          buf = fd.read()

      re_level = re.compile(r'(Level:)\s+(\w+)')
      m = re_level.search(buf)
      level = 'none'
      if m:
            level = m.group(2)
      else:
            print('Error! No level info in file:', filename)

      # Reformat level and location
      tmpbuf = re_level.sub('',buf)
      re_loc = re.compile('(<FONT COLOR="#883300">Location: </FONT>)')
      tmpbuf = re_loc.sub('</B><H3><FONT COLOR="#883300">Level</FONT></H3>' + level + r'<BR>\n<H3><FONT COLOR="#883300">Location</FONT></H3>\n',tmpbuf)

      # Modify .c#,.h#,.cu#,.cxx# to .c.html#,.h.html#,.cu.html#,.cxx.html#
      tmpbuf = re.sub('.c#', '.c.html#', tmpbuf)
      tmpbuf = re.sub('.h#', '.h.html#', tmpbuf)
      tmpbuf = re.sub('.cu#', '.cu.html#', tmpbuf)
      tmpbuf = re.sub('.cxx#', '.cxx.html#', tmpbuf)

      re_loc = re.compile('</BODY></HTML>')
      outbuf = re_loc.sub('<BR><BR><A HREF="./index.html">Index of all ' + secname + ' routines</A>\n<BR><A HREF="../../../docs/manual.html">Table of Contents for all manual pages</A>\n<BR><A HREF="../singleindex.html">Index of all manual pages</A>\n</BODY></HTML>',tmpbuf)

      re_loc = re.compile(r' (http://[A-Za-z09_\(\)\./]*)[ \n]')
      outbuf = re_loc.sub(' <a href="\\1">\\1 </a> ',outbuf)

      # write the modified manpage
      with open(filename, "w") as fd:
          fd.write(outbuf)

      return level

# Go through each manpage file, present in dirname,
# and create and return a table for it, wrt levels specified.
def createtable(dirname,levels,secname):
      htmlfiles = [os.path.join(dirname,f) for f in os.listdir(dirname) if f.endswith('.html')]
      htmlfiles.sort()
      if htmlfiles == []:
            print('Error! Empty directory:',dirname)
            return None

      table = []
      for level in levels: table.append([])

      for filename in htmlfiles:
            level = modifylevel(filename,secname)
            #if not level: continue
            if level.lower() in levels:
                  table[levels.index(level.lower())].append(filename)
            else:
                  print('Error! Unknown level \''+ level + '\' in', filename)
      return table

# This routine is called for each man dir. Each time, it
# adds the list of manpages, to the given list, and returns
# the union list.

def addtolist(dirname,singlelist):
      htmlfiles = [os.path.join(dirname,f) for f in os.listdir(dirname) if f.endswith('.html')]
      htmlfiles.sort()
      if htmlfiles == []:
            print('Error! Empty directory:',dirname)
            return None

      for filename in htmlfiles:
            singlelist.append(filename)

      return singlelist

# This routine creates a dictionary, with entries such that each
# key is the alphabet, and the vaue corresponds to this key is a dictionary
# of FunctionName/PathToFile Pair.
def createdict(singlelist):
      newdict = {}
      for filename in singlelist:
            path,name     = posixpath.split(filename)
            # grab the short path Mat from /wired/path/Mat
            junk,path     = posixpath.split(path)
            index_char    = name[0:1].lower()
            # remove the .name suffix from name
            func_name,ext = posixpath.splitext(name)
            if index_char not in newdict:
                  newdict[index_char] = {}
            newdict[index_char][func_name] = path + '/' + name

      return newdict


def getallmandirs(dirs):
      """ Gets the list of man* dirs present in the doc dir. Each dir will have an index created for it. """
      mandirs = []
      for filename in dirs:
            path,name = posixpath.split(filename)
            if name == 'RCS' or name == 'sec' or name == "concepts" or name  == "SCCS" : continue
            if posixpath.isdir(filename):
                  mandirs.append(filename)
      return mandirs


# Extracts PETSC_DIR from the command line and
# starts genrating index for all the manpages.
def main():
      arg_len = len(sys.argv)

      if arg_len < 3:
            print('Error! Insufficient arguments.')
            print('Usage:', sys.argv[0], 'PETSC_DIR','LOC')
            exit()

      PETSC_DIR = sys.argv[1]
      LOC       = sys.argv[2]
      HEADERDIR = (sys.argv[3] if arg_len > 3 else 'doc/manualpages/MANSECHeaders')
      #fd        = os.popen('/bin/ls -d '+ PETSC_DIR + '/manualpages/*')
      #buf       = fd.read()
      #dirs      = split(strip(buf),'\n')
      dirs      = glob.glob(LOC + '/docs/manualpages/*')
      mandirs   = getallmandirs(dirs)

      levels = ['beginner','intermediate','advanced','developer','deprecated','none']
      titles = ['Beginner - Basic usage',
                'Intermediate - Setting options for algorithms and data structures',
                'Advanced - Setting more advanced options and customization',
                'Developer - Interfaces intended primarily for library developers, not for typical applications programmers',
                'Deprecated - Functionality scheduled for removal in future versions',
                'None: Not yet cataloged']

      singlelist = []
      for dirname in mandirs:
            outfilename  = dirname + '/index.html'
            dname,secname  = posixpath.split(dirname)
            headfilename = PETSC_DIR + '/' + HEADERDIR + '/' + secname
            table        = createtable(dirname,levels,secname)
            if not table: continue
            singlelist   = addtolist(dirname,singlelist)
            printindex(outfilename,headfilename,levels,titles,table)

      alphabet_dict = createdict(singlelist)
      outfilename   = LOC + '/docs/manualpages/singleindex.html'
      printsingleindex (outfilename,alphabet_dict)


if __name__ == '__main__':
      main()

