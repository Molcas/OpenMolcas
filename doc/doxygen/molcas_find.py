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
#***********************************************************************

from __future__ import print_function
import os.path
import sys
import re
import xml.etree.ElementTree as ET

basedir = 'xml'

# Read argument
try:
  argument = sys.argv[1]
except IndexError:
  argument = None

database = ET.parse(os.path.join(basedir,'index.xml'))
root = database.getroot()

# First create a list of files, containing the directory where each find is located.
file_list={}
# Iterate directory items
for diritem in root.iterfind('compound[@kind="dir"]'):
  dirname = diritem.find('name').text
  # Open the XML file for each directory
  dirxml = diritem.get('refid') + '.xml'
  dirtree = ET.parse(os.path.join(basedir,dirxml))
  # For each file in this directory, assign the director name
  for fileitem in dirtree.getroot().iter('innerfile'):
    file_list[fileitem.get('refid')] = dirname

# Create a dictionary linking each element with its parent
parent_map = {c:p for p in root.iter() for c in p}

# Build the dictionary of functions, containg the file where each function is defined
name_list={}
# Iterate all "function" members belonging to "file" compounds
for member in root.iterfind('compound[@kind="file"]/member[@kind="function"]'):
  name = member.find('name').text
  # Find the directory and filename
  parent = parent_map[member]
  filename = parent.find('name').text
  # Ignore .h files
  if (not re.search(r"\.h$",filename)):
    dirname = file_list[parent.get('refid')]
    filename = os.path.join(dirname,filename)
    # Add an entry to the dictionary
    name_list.setdefault(name,[]).append(filename)

#===============================================================================

# No argument, print the full list
if (argument == None):
  for (name,files) in sorted(name_list.items()):
    for filename in sorted(files):
      print("{}:{}".format(name,filename))

# -Conf, find conflicts
elif (re.match('-Conf',argument)):
  check_list = {}
  expr = re.compile(r"(clones|^ga|^ma)")
  for (name,files) in sorted(name_list.items()):
    # Ignore files matching the above pattern
    for filename in sorted(files,reverse=True):
      if expr.search(filename):
        files.remove(filename)
    if (len(files) > 1):
      print("Conflict: {} found in {}".format(name,", ".join(files)))

# Search name, print the names matching the argument
else:
  expr = re.compile(argument)
  for name in sorted(filter(expr.search,name_list)):
    for filename in sorted(name_list[name]):
      print("{}:{}".format(name,filename))
