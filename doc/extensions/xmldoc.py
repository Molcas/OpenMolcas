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
# Copyright (C) 2015,2018, Ignacio Fdez. Galván                        *
#***********************************************************************

try:
  from docutils.parsers.rst import Directive
  do_directive = True
except ImportError:
  do_directive = False
import re
import codecs
import os.path

X_clean = re.compile(r'^%+\+*')
X_line = re.compile(r'(<\/?(help|key|group|select|module|emil|command).*?>)', flags=re.IGNORECASE)
X_inhelp = re.compile(r'<help>', flags=re.IGNORECASE)
X_inhelp2 = re.compile(r'<\/help>', flags=re.IGNORECASE)
X_comment = re.compile(r'\s*<!--.*?-->')
X_tag = re.compile(r'<')
X_help = re.compile(r'(<\/?help)', flags=re.IGNORECASE)
X_key = re.compile(r'(<\/?(key|command))', flags=re.IGNORECASE)
X_group = re.compile(r'(<\/?(group|select))', flags=re.IGNORECASE)

H_head = re.compile(r'%%(keyword|description):', flags=re.IGNORECASE)
H_remove = re.compile(r'(%|<!--$|-->$)')
H_keyw = re.compile(r'%%keyword:\s*(.*?)\s*<(.*?)>\s*(.*)', flags=re.IGNORECASE)
H_desc = re.compile(r'%%description', flags=re.IGNORECASE)

nowrap = re.compile(br'^\|+(.*)'.decode('raw_unicode_escape'))
nowrap2 = re.compile(br'(\u00a6.*?\u00a6)'.decode('raw_unicode_escape'))
nowrap3 = re.compile(br'\u00a6(.*?)\u00a6'.decode('raw_unicode_escape'))

# Create the xmldoc directive
#
if (do_directive):
  class XMLDocDirective(Directive):
    has_content = True

    def run(self):
      self.assert_has_content()
      # Add a global "XMLDocs" attribute to the environment
      env = self.state.document.settings.env
      if not hasattr(env, 'XMLDocs'):
        env.XMLDocs = []
      # The "XMLDocs" attribute contains all the xmldoc pieces
      env.XMLDocs.append({
        'docname': env.docname,
        'lineno': self.lineno,
        'content': [line for line in self.content],
      })
      return []

# Write all the XML pieces at the end
#
def write_XMLDocs(app, exception):
  env = app.builder.env
  if hasattr(env, 'XMLDocs'):
    data_dir = app.config.data_dir
    with codecs.open(os.path.join(data_dir, 'keyword.xml'), 'w', 'utf-8') as keywordsfile:
      keywordsfile.write('<ROOT>\n')
      # Sort by docname, then by lineno
      docs = list(set([piece['docname'] for piece in env.XMLDocs]))
      for doc in sorted(docs):
        for piece in sorted([pc for pc in env.XMLDocs if pc['docname'] == doc], key=lambda x: x['lineno']):
          # Reformat the text for XML
          text = reformat_XML(piece['content'])
          if (text):
            keywordsfile.write('\n'.join(text)+'\n')
      keywordsfile.write('</ROOT>\n')

# Reformat before writing to the XML
#
def reformat_XML(piece):
  import textwrap
  # Remove lines we do not want in the XML, and do some initial cleanup
  text = ' '.join([nowrap.sub(br'\u00a6\1\u00a6'.decode('raw_unicode_escape'),X_clean.sub('',line)) for line in piece if not H_head.match(line)])
  # XML tags in separate lines
  text = X_line.sub(r'\n\1\n', text)
  # no-wrap lines in separate lines
  text = nowrap2.sub(r'\n\1\n', text)
  text = text.split('\n')
  inhelp = False
  for i in range(len(text)):
    line = text[i].strip()
    if (X_inhelp.match(line)): inhelp = True
    if (X_inhelp2.match(line)): inhelp = False
    # Wrap text lines inside <HELP>, discard them outside
    if (inhelp):
      if (line and not nowrap2.match(line)):
        line = textwrap.fill(X_comment.sub('', line), 59, initial_indent=' ', subsequent_indent=' ', break_on_hyphens=False)
    elif (not X_tag.match(line)):
      line = ''
    # Indent XML tags
    line = X_help.sub(r'         \1', line)
    line = X_key.sub(r'      \1', line)
    line = X_group.sub(r'   \1', line)
    text[i] = line
  # Remove empty lines
  text = [nowrap3.sub(r' \1',line) for line in filter(bool, text)]
  return text

# Write all the xmldoc pieces at the end
#
def write_Help(app, exception):
  env = app.builder.env
  if hasattr(env, 'XMLDocs'):
    data_dir = app.config.data_dir
    with codecs.open(os.path.join(data_dir, 'keyword.db'), 'w', 'utf-8') as keywordsfile:
      keywordsfile.write('#This file generated automatically from MOLCAS documentation\n')
      keywordsfile.write('#\n')
      # Sort by docname, then by lineno
      docs = list(set([piece['docname'] for piece in env.XMLDocs]))
      for doc in sorted(docs):
        keywordsfile.write('[%s]\n' % os.path.basename(doc).upper())
        for piece in sorted([pc for pc in env.XMLDocs if pc['docname'] == doc], key=lambda x: x['lineno']):
          # Reformat the text for the Help file
          text = reformat_Help(piece['content'])
          if (text):
            keywordsfile.write('\n'.join(text)+'\n')

def reformat_Help(piece):
  import textwrap
  # Find the start of the help section
  index = None
  header = None
  for i in range(len(piece)):
    if (H_head.match(piece[i])):
      index = i
      header = help_header(piece[index])
      break 
  if (index is None):
    text = ''
  else:
    # Remove lines we do not want in the help file, and do some initial cleanup
    text = ' '.join([nowrap.sub(br'\u00a6\1\u00a6'.decode('raw_unicode_escape'),X_clean.sub('',line)) for line in piece[index:] if not (H_remove.match(line) or X_line.match(line))])
  # no-wrap lines in separate lines
  text = nowrap2.sub(r'\n\1\n', text)
  text = text.split('\n')
  for i in range(len(text)):
    line = text[i].strip()
    # Wrap text lines
    if (line and not nowrap2.match(line)):
      line = textwrap.fill(line, 59, initial_indent=' ', subsequent_indent=' ', break_on_hyphens=False)
    text[i] = line
  # Remove empty lines
  text = [nowrap3.sub(r' \1',line) for line in filter(bool, text)]
  if (header): text.insert(0, header)
  return text

def help_header(line):
  if (H_desc.match(line)):
    return ':Description:D'
  elif (H_keyw.match(line)):
    tmp = H_keyw.match(line)
    typ = tmp.group(2)[0].upper()
    if (typ == 'C'): typ = 'B'
    gui = tmp.group(3)
    if (not gui): gui = ' '
    return ':%s:%s:%s:' % (tmp.group(1), typ, gui)
  else:
    return ''

# Purge the xmldoc pieces from a file when it is changed
#
def purge_XMLDocs(app, env, docname):
  if hasattr(env, 'XMLDocs'):
    env.XMLDocs = [piece for piece in env.XMLDocs
                   if piece['docname'] != docname]

# Setup
#
if (do_directive):
  def setup(app):
    app.add_directive('xmldoc', XMLDocDirective)
    app.connect('env-purge-doc', purge_XMLDocs)
    app.connect('build-finished', write_XMLDocs)
    app.connect('build-finished', write_Help)
    app.add_config_value('data_dir', 'xmldoc_files', '')
