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
# Lots of problems with this script!!!

from __future__ import print_function
from __future__ import unicode_literals
import os.path
import sys
import re
import xml.etree.ElementTree as ET
import htmlentitydefs

basedir = 'xml'
extension = '.xml'
htmldir = 'html'
htmlextension = '.xhtml'

# Read argument
try:
  argument = sys.argv[1]
except IndexError:
  argument = None

def in_red(text):
  return '\033[31m{}\033[0m'.format(text)

def in_brown(text):
  return '\033[33m{}\033[0m'.format(text)

def in_blue(text):
  return '\033[34m{}\033[0m'.format(text)

database = ET.parse(os.path.join(basedir,'index'+extension))
root = database.getroot()

# First create a list of files, containing the directory where each find is located.
file_list={}
# Iterate directory items
for diritem in root.iterfind('compound[@kind="dir"]'):
  dirname = diritem.findtext('name')
  # Open the XML file for each directory
  dirxml = diritem.get('refid') + extension
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
  name = member.findtext('name')
  refid = member.get('refid')
  # Find the directory and filename
  parent = parent_map[member]
  filename = parent.findtext('name')
  fileref = parent.get('refid')
  # Ignore .h files
  if (not re.search(r"\.h$",filename)):
    dirname = file_list[parent.get('refid')]
    filename = os.path.join(dirname,filename)
    # Add an entry to the dictionary
    name_list.setdefault(name,[]).append([filename,fileref,refid])

#===============================================================================

def get_info(fileref,refid):
  import HTMLParser
  h = HTMLParser.HTMLParser()
  filedata = ET.parse(os.path.join(basedir,fileref+extension))
  root = filedata.getroot()
  memberdef = root.find('.//memberdef[@id="{}"]'.format(refid))
  # Insert Unicode characters
  for item in memberdef.findall('.//grave'):
    item.text = h.unescape('&{}grave;'.format(item.get('char')))
  for item in memberdef.findall('.//acute'):
    item.text = h.unescape('&{}acute;'.format(item.get('char')))
  for item in memberdef.findall('.//circ'):
    item.text = h.unescape('&{}circ;'.format(item.get('char')))
  for item in memberdef.findall('.//tilde'):
    item.text = h.unescape('&{}tilde;'.format(item.get('char')))
  for item in memberdef.findall('.//umlaut'):
    item.text = h.unescape('&{}uml;'.format(item.get('char')))
  for item in memberdef.findall('.//ring'):
    item.text = h.unescape('&{}ring;'.format(item.get('char')))
  for item in memberdef.findall('.//cedil'):
    item.text = h.unescape('&{}cedil;'.format(item.get('char')))
  for item in memberdef.findall('.//slash'):
    item.text = h.unescape('&{}slash;'.format(item.get('char')))
  for item in memberdef.findall('.//larr'):
    item.text = h.unescape('&larr;')
  for item in memberdef.findall('.//rarr'):
    item.text = h.unescape('&rarr;')
  for item in memberdef.findall('.//times'):
    item.text = h.unescape('&times;')
  for item in memberdef.findall('.//le'):
    item.text = h.unescape('&le;')
  for item in memberdef.findall('.//ge'):
    item.text = h.unescape('&ge;')
  for item in memberdef.findall('.//ne'):
    item.text = h.unescape('&ne;')
  for item in memberdef.findall('.//sp'):
    item.text = ' '
  # Find detailed section
  detailednode = memberdef.find('detaileddescription')
  if (len(''.join(detailednode.itertext()).strip()) == 0): detailednode = memberdef.find('inbodydescription')
  # Process definition
  definition = memberdef.findtext('definition')
  argsstring = memberdef.findtext('argsstring')
  print(' {}: {}{}'.format(in_blue('Definition'),definition,argsstring))
  # Process arguments
  parameterlist = memberdef.find('.//parameterlist')
  params = []
  size = 0
  # Get parameter definitions from XML
  if (parameterlist is not None):
    for item in parameterlist.iterfind('parameteritem'):
      names = []
      for param in item.iterfind('parameternamelist/parametername'):
        if (param.get('direction') is None):
          names.append(param.text)
        else:
          names.append('{} [{}]'.format(param.text,param.get('direction')))
      names = ', '.join(names)
      description = ''.join(item.find('parameterdescription').itertext())
      description = description.strip()
      description = re.sub(r"\n",' ',description)
      params.append([names,description])
      size = max(size,len(names))
  # Get only parameter names from XML
  else:
    for item in memberdef.iterfind('param'):
      name = item.findtext('defname')
      if (name is None): name = item.findtext('declname')
      params.append([name,''])
      size = max(size,len(name))
  # Get parameter types from HTML (due to bug in doxygen's XML output)
  newrefid = re.match(r"(.*)_",refid).group(1)
  parser = ET.XMLParser()
  parser.parser.UseForeignDTD(True)
  parser.entity.update(htmlentitydefs.entitydefs)
  filedata = ET.parse(os.path.join(htmldir,newrefid+htmlextension), parser)
  root = filedata.getroot()
  # Since some parameters may be listed together, we need some kind of index
  index = []
  for i in range(len(params)):
    for j in range(params[i][0].count(', ')+1):
      index.append(i)
  # The types are read from the XML and added to the description
  i = 0
  for item in root.findall('.//*[@class="memproto"]//*[@class="paramtype"]'):
    typ = re.sub(r"(, $|\xa0)",'',''.join(item.itertext()))
    # Define implicit types (assuming default implicit)
    if (len(typ) == 0):
      if (re.match(r"[ijklmn]",params[i][0][0], flags=re.IGNORECASE)):
        typ = '(integer)?'
      else:
        typ = '(real*8)?'
    params[index[i]][1] += '  {}'.format(in_brown(typ))
    i += 1
  # Print arguments
  print(' {}:'.format(in_blue('Arguments')))
  formatstring = '  {{:{}s}}: {{}}'.format(size)
  for (names,description) in params:
    print(formatstring.format(names,description))
  print()
  # Process brief description
  briefdescription = ''.join(memberdef.find('briefdescription').itertext()).strip()
  if (len(briefdescription) > 0):
    print(' {}:'.format(in_blue('Brief description')))
    print('  {}'.format(briefdescription))
  # Process authors
  authors = []
  for item in memberdef.findall('.//simplesect[@kind="author"]/para'):
    authors.append(''.join(item.itertext()).strip())
  for i in range(len(authors)):
    if (i == 0):
      print(' {}: {}'.format(in_blue('Author'),authors[i]))
    else:
      print('         {}'.format(authors[i]))
  if (len(authors) > 0): print()
  # Remove unwanted items from the detailed description
  for item in detailednode.iter():
    if (item.tag in ('parameterlist','simplesect')):
      tail = item.tail
      item.clear()
      item.tail = tail
  detaileddescription = ''.join(detailednode.itertext()).strip()
  if (len(detaileddescription) > 0):
    print(' {}:'.format(in_blue('Detailed description')))
    print('  {}'.format(detaileddescription))
  print()

#===============================================================================

# No argument, print the full list
if (argument == None):
  pass

# Search name, print the names matching the argument
else:
  name = argument
  if (name in name_list.keys()):
    for (filename,fileref,refid) in (sorted(name_list[name], key=lambda x: x[0])):
      print('  {} in {}'.format(in_red(name),in_red(filename)))
      get_info(fileref,refid)
  else:
    print('Name {} not found'.format(in_red(name)))
