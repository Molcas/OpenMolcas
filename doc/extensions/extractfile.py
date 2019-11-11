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
# Copyright (C) 2015, Ignacio Fdez. Galván                             *
#***********************************************************************

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from sphinx.util.nodes import set_source_info
from sphinx.directives.code import container_wrapper, CodeBlock
import codecs
import os, os.path

# Create the extractfile directive
#
class ExtractFileDirective(Directive):
  has_content = True
  required_arguments = 1
  optional_arguments = 0
  final_argument_whitespace = False
  option_spec = {
  }

  def run(self):
    code = u'\n'.join(self.content)
    literal = nodes.literal_block(code, code)
    filename = self.arguments[0]
    addCode(self.state.document.settings.env, filename, code)
    return [literal]

# This is just a copy of the CodeBlock directive,
# with the filename argument and the additional extracting code
#
class ExtractCodeBlock(CodeBlock):
  has_content = True
  required_arguments = 0
  optional_arguments = 1
  final_argument_whitespace = False
  option_spec = {
    'force': directives.flag,
    'linenos': directives.flag,
    'dedent': int,
    'lineno-start': int,
    'emphasize-lines': directives.unchanged_required,
    'caption': directives.unchanged_required,
    'class': directives.class_option,
    'name': directives.unchanged,
    ###########################################
    'filename': directives.unchanged_required,
    ###########################################
  }
  def run(self):
    document = self.state.document
    code = '\n'.join(self.content)
    location = self.state_machine.get_source_and_line(self.lineno)
    linespec = self.options.get('emphasize-lines')
    if linespec:
      try:
        nlines = len(self.content)
        hl_lines = parselinenos(linespec, nlines)
        if any(i >= nlines for i in hl_lines):
          logger.warning(__('line number spec is out of range(1-%d): %r') %
                         (nlines, self.options['emphasize-lines']),
                         location=location)
        hl_lines = [x + 1 for x in hl_lines if x < nlines]
      except ValueError as err:
        return [document.reporter.warning(err, line=self.lineno)]
    else:
      hl_lines = None
    if 'dedent' in self.options:
      location = self.state_machine.get_source_and_line(self.lineno)
      lines = code.split('\n')
      lines = dedent_lines(lines, self.options['dedent'], location=location)
      code = '\n'.join(lines)
    literal = nodes.literal_block(code, code)  # type: Element
    if 'linenos' in self.options or 'lineno-start' in self.options:
      literal['linenos'] = True
    literal['classes'] += self.options.get('class', [])
    literal['force'] = 'force' in self.options
    if self.arguments:
      literal['language'] = self.arguments[0]
    else:
      literal['language'] = self.env.temp_data.get('highlight_language',
                                                   self.config.highlight_language)
    extra_args = literal['highlight_args'] = {}
    if hl_lines is not None:
      extra_args['hl_lines'] = hl_lines
    if 'lineno-start' in self.options:
      extra_args['linenostart'] = self.options['lineno-start']
    self.set_source_info(literal)
    caption = self.options.get('caption')
    if caption:
      try:
        literal = container_wrapper(self, literal, caption)
      except ValueError as exc:
        return [document.reporter.warning(exc, line=self.lineno)]
    self.add_name(literal)
    ###########################################
    filename = self.options.get('filename')
    addCode(self.state.document.settings.env, filename, code)
    ###########################################
    return [literal]
  def set_source_info(self, literal):
    try:
      super(ExtractCodeBlock, self).set_source_info(literal)
    except AttributeError:
      set_source_info(self, literal)

# Add code to the environment
#
def addCode (env, filename, code):
  # Add a global "FilesToExtract" attribute to the environment
  if not hasattr(env, 'FilesToExtract'):
    env.FilesToExtract = []
  filenames = [item['filename'] for item in env.FilesToExtract]
  i = 0
  # Prevent multiple files with the same name
  while (filename in filenames):
    i += 1
    filename = '%s.%i' % (filename, i)
  env.FilesToExtract.append({
    'docname': env.docname,
    'filename': filename,
    'content': code,
  })

# Extract the files to extract
#
def extractFiles (app, exception):
  env = app.builder.env
  if hasattr(env, 'FilesToExtract'):
    ext_dir = app.config.extract_dir
    for item in env.FilesToExtract:
      filename = os.path.join(ext_dir, item['filename'])
      directory = os.path.dirname(filename)
      if (not os.path.exists(directory)):
          os.makedirs(directory)
      try:
        with codecs.open(filename, 'w', 'utf-8') as extractedfile:
          extractedfile.write(item['content'])
      except:
        print ('Error creating file %s' % filename)

# Purge the extractfile items from a file when it is changed
#
def purge_FilesToExtract(app, env, docname):
  if hasattr(env, 'FilesToExtract'):
    env.FilesToExtract = [item for item in env.FilesToExtract
                          if item['docname'] != docname]

# Setup
#
def setup(app):
  app.add_directive('extractfile', ExtractFileDirective)
  app.add_directive('extractcode-block', ExtractCodeBlock)
  app.connect('env-purge-doc', purge_FilesToExtract)
  app.connect('build-finished', extractFiles)
  app.add_config_value('extract_dir', 'extracted_files', '')
