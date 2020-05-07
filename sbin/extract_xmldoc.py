#!/usr/bin/env python
# -*- coding: utf-8 -*-

#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2018, Ignacio Fdez. Galván                             *
#***********************************************************************

import sys
import os
import codecs
import os.path
import re
import textwrap

sys.dont_write_bytecode = True

try:
  docdir = os.path.abspath(sys.argv[1])
except IndexError:
  print('Need doc directory location')
  sys.exit(1)

try:
  datadir = os.path.abspath(sys.argv[2])
except IndexError:
  print('Need data directory location')
  sys.exit(1)

sys.path.append(os.path.join(docdir, 'extensions'))

import xmldoc

xml_re = re.compile(r'(\s*)..\s+xmldoc::\s*')

# Find all ".. xmldoc::" blocks in the source files
xmldocs = {}
for root, subdirs, files in os.walk(os.path.join(docdir, 'source')):
  for fname in files:
    if (not fname.endswith('.rst')):
      continue
    docs = []
    with open(os.path.join(root, fname), 'rb') as f:
      inxml = False
      lines = f.readlines()
      for l in lines:
        try:
          line = l.decode('utf8').rstrip()
        except UnicodeDecodeError:
          continue
        if (inxml):
          match = re.match(r'(\s{{{0},}})'.format(inxml), line)
          if (match or (line == '')):
            try:
              docs[-1].append(line[len(match.group(1)):])
            except AttributeError:
              docs[-1].append('')
          else:
            inxml = False
        if (not inxml):
          match = xml_re.match(line)
          if match:
            inxml = len(match.group(0))
            docs.append([line[inxml:]])
            inxml = len(match.group(1))+1
    if (docs):
      xmldocs[f.name] = docs

# Create Help file
with codecs.open(os.path.join(datadir, 'keyword.db'), 'w', 'utf-8') as keywordsfile:
  keywordsfile.write('#This file is generated automatically from the OpenMolcas documentation\n')
  keywordsfile.write('#\n')
  # Sort by filename
  for fname in sorted(xmldocs):
    keywordsfile.write('[%s]\n' % os.path.splitext(os.path.basename(fname))[0].upper())
    for doc in xmldocs[fname]:
      # Reformat the text for the Help file
      text = xmldoc.reformat_Help(doc)
      if (text):
        keywordsfile.write('\n'.join(text)+'\n')

# Create XML file
with codecs.open(os.path.join(datadir, 'keyword.xml'), 'w', 'utf-8') as keywordsfile:
  keywordsfile.write('<!-- This file is generated automatically from the OpenMolcas documentation -->\n')
  keywordsfile.write('<ROOT>\n')
  # Sort by filename
  for fname in sorted(xmldocs):
    for doc in xmldocs[fname]:
      # Reformat the text for XML
      text = xmldoc.reformat_XML(doc)
      if (text):
        keywordsfile.write('\n'.join(text)+'\n')
  keywordsfile.write('</ROOT>\n')

