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
# Copyright (C) 2020, Ignacio Fdez. Galván                             *
#***********************************************************************

import re
import textwrap
from xml.etree import ElementTree as ET

# Return the number of leading spaces and the "item label" format, if any
# The "item label" is the first word if it is an asterisk or a single
# letter or number followed by . or )
#
def _indentlabel(line):
  match = re.match(r'(\s*)(\S*)?', line)
  indent = len(match.group(1))
  label = match.group(2)
  match = re.match(r'[a-zA-Z0-9][.)]$', label)
  if label == '*':
    pass
  elif match:
    if re.match(r'[a-z]', label):
      label = 'a' + label[1]
    elif re.match(r'[A-Z]', label):
      label = 'A' + label[1]
    elif re.match(r'[0-9]', label):
      label = '1' + label[1]
  else:
    match = re.match(r'(\s*)(\S+\s+-+)', line)
    if match:
      label = len(match.group(2))*'_'
    else:
      label = ''
  return indent, label

# Return the number of spaces that should be used for "hanging indent",
# i.e., the number of leading spaces (plus the length of the item label plus one)
#
def _hangindent(label):
  if label[2]:
    return label[1] + len(label[2]) + 1
  else:
    return label[1]

# Finish a paragraph-level block, by assigning the current buffer text
#
def _finish_paragraph(element, paragraph):
  if element is not None:
    element.text = '\n'.join(paragraph).rstrip()
  elif paragraph:
    raise Exception
  return None, []

# Append an element updating the parent map
#
def _append_element(element, parent, parent_map):
  parent.append(element)
  parent_map[element] = parent

# Parse a help text, returning an HTML-like tree
#
def parse_help_text(text):
  el = None
  par = []
  level = [(None, -1, '')]
  newpar = False
  doc = ET.Element('root')
  parent_map = {doc: None}

  for line in text.split('\n'):
    ind, lab = _indentlabel(line)
    blank = not line.strip()
    # A pre-formatted block continues until a line with less indent is found
    if el is not None and el.tag == 'pre' and (ind >= level[-1][1] or blank):
      par.append(line)
      continue
    # A paragraph-block is finished by a blank line
    # or a line that does not match the current hanging indent level
    if blank or ind != _hangindent(level[-1]):
      el, par = _finish_paragraph(el, par)
      newpar = True
    # If a new paragraph starts, find out which type and level
    if not blank and newpar:
      newpar = False
      # Start from innermost level and go outwards, until a match is found
      for i in range(len(level)):
        parent = level[-1][0]
        if parent is None:
          parent = doc
        # Normal paragraph at the current level
        if not lab and ind == _hangindent(level[-1]):
          break
        # List item at the current level
        elif (ind, lab) == level[-1][1:3]:
          # Create a new <li> element in the parent list, and update the current parent
          li = ET.Element('li')
          _append_element(li, parent_map[parent], parent_map)
          parent = li
          level[-1] = (parent,) + level[-1][1:]
          break
        # New list, or more indent
        elif (lab and ind == level[-1][1]) or ind > level[-1][1]:
          # New list: create new <ul> and <li> elements
          if lab:
            ul = ET.Element('ul')
            ul.attrib['label'] = lab
            _append_element(ul, parent, parent_map)
            parent = ul
            li = ET.Element('li')
            _append_element(li, parent, parent_map)
            parent = li
          # Indented paragraph: child of the last added element
          elif len(parent) > 0:
            parent = parent[-1]
          level.append((parent, ind, lab))
          break
        del level[-1]
      # A child of <p> is actually a sibling <pre>
      if parent.tag == 'p':
        parent = parent_map[parent]
        el = ET.Element('pre')
      # Everything else is a <p> (possibly inside <li>)
      else:
        el = ET.Element('p')
      _append_element(el, parent, parent_map)
    # Whatever the case, add the current line to the buffer
    if not blank:
      par.append(line)
  # Finish the last element
  el, par = _finish_paragraph(el, par)

  return doc

# Format an HTML element as plain text
#
def _format_element(element, parent_map, textwidth, indent):
  out = ''
  # Increase indent in lists by multiples of 2
  if element.tag in ['ol', 'ul']:
    indent += 2
  # Increase indent of preformatted texts by 4
  elif element.tag == 'pre':
    indent += 4
  # If there are child elements, call recursively
  if len(element) > 0:
    for el in element:
      out += _format_element(el, parent_map, textwidth, indent)
  # If there are no child elements, format the inner text
  else:
    # Pre-formatted text: Verbatim lines (replacing common indentation)
    if element.tag == 'pre':
      for line in textwrap.dedent(element.text).split('\n'):
        if line:
          out += indent*' ' + line + '\n'
        else:
          out += '\n'
    # Other elements, use textwrap
    elif element.text:
      tw = textwrap.TextWrapper(break_long_words=False, break_on_hyphens=False, width=textwidth)
      tw.initial_indent = indent*' '
      tw.subsequent_indent = tw.initial_indent
      # List items: Hanging indent in first child, extra indent in the rest
      parent = parent_map[element]
      text = element.text.lstrip()
      lab = ''
      if parent.tag == 'li':
        hang = len(parent_map[parent].attrib['label'])
        tw.subsequent_indent += (hang+1)*' '
        if element is parent[0]:
          # Split the label, so the spaces here are not altered
          lab = text[0:hang] + ' '
          text = text[hang:]
        else:
          tw.initial_indent += hang*' '
      text = lab + ' '.join(text.split())
      out += tw.fill(text) + '\n'
    out += '\n'

  return out

# Return a plain-text-formatted text from an HTML-like tree
#
def format_help_text(doc, textwidth=80):
  parent_map = {e:p for p in doc.iter() for e in p}
  return _format_element(doc, parent_map, textwidth, 0).strip()
