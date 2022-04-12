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
# Copyright (C) 2015,2016,2022, Ignacio Fdez. Galván                   *
#***********************************************************************

import sphinx
from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.statemachine import ViewList
from sphinx.util.nodes import set_source_info
from sphinx.builders.latex.nodes import captioned_literal_block

class container(nodes.General, nodes.Element): pass
class float_container(container): pass
class table_container(container): pass
class figure_container(container): pass
class code_container(container): pass

# Create the float directive
#
class FloatDirective(Directive):
  has_content = True
  required_arguments = 0
  optional_arguments = 0
  final_argument_whitespace = False
  option_spec = {
    'type': directives.unchanged_required,
    'caption': directives.unchanged_required,
    'caption-top': directives.flag,
    'name': directives.unchanged,
  }

  def run(self):
    if (self.options.get('type') == 'table'):
      content_node = table_container('')
    elif (self.options.get('type') == 'figure'):
      content_node = figure_container('')
    elif (self.options.get('type') == 'code'):
      content_node = code_container('')
    else:
      content_node = float_container('')
    self.state.nested_parse(self.content, self.content_offset, content_node)
    set_source_info(self, content_node)
    self.add_name(content_node)
    caption = self.options.get('caption')
    if caption:
      parsed = nodes.Element()
      self.state.nested_parse(ViewList([caption], source=''), self.content_offset, parsed)
      caption_node = nodes.caption(parsed[0].rawsource, '', *parsed[0].children)
      caption_node.source = parsed[0].source
      caption_node.line = parsed[0].line
      if ('caption-top' in self.options):
        content_node.insert(0, caption_node)
      else:
        content_node += caption_node
    return [content_node]

def visit_float_html(self, node):
  self.body.append(self.starttag(node, 'div', CLASS='float-wrapper'))

def depart_float_html(self, node):
  self.body.append('</div>\n')

def visit_float_latex(self, node):
  ids = ''
  if node['ids']:
    ids += self.hypertarget(node['ids'][0], anchor=False)
  if (isinstance(node, table_container)):
    floattype = 'table'
    self.body.append('\\sphinxcapstartof{%s}\n' % floattype)
    self.context.append(ids + '\n')
  else:
    if (isinstance(node, figure_container)):
      floattype = 'figure'
    elif (isinstance(node, code_container)):
      floattype = 'code'
    else:
      floattype = 'float'
    self.body.append('\\begin{%s}\n' % floattype)
    if any(isinstance(child, nodes.caption) for child in node):
      self.body.append('\\capstart\n')
    self.context.append(ids + '\\end{%s}\n' % floattype)

def depart_float_latex(self, node):
  self.body.append(self.context.pop())

# Patch caption formatting in containers
#
visit_caption_orig = sphinx.writers.latex.LaTeXTranslator.visit_caption
depart_caption_orig = sphinx.writers.latex.LaTeXTranslator.depart_caption
def visit_caption_patched(self, node: nodes.Element) -> None:
    if isinstance(node.parent, container):
      if node.parent[0] is node:
        self.body.append('\\sphinxthecaptionisattop\n')
      self.in_caption += 1
      if isinstance(node.parent, table_container):
        self.body.append('\\sphinxcaption{')
      else:
        self.body.append('\\caption{')
    else:
      visit_caption_orig(self, node)
def depart_caption_patched(self, node: nodes.Element) -> None:
  if isinstance(node.parent, container):
    self.body.append('}')
    if node.parent[0] is node:
      self.body.append('\n\\sphinxaftertopcaption')
    self.in_caption -= 1
  else:
    depart_caption_orig(self, node)
sphinx.writers.latex.LaTeXTranslator.visit_caption = visit_caption_patched
sphinx.writers.latex.LaTeXTranslator.depart_caption = depart_caption_patched

# Setup
#
def setup(app):
  app.add_directive('float', FloatDirective)
  app.add_node(float_container,
               html=(visit_float_html, depart_float_html),
               latex=(visit_float_latex, depart_float_latex))
  app.add_enumerable_node(table_container, 'table',
               html=(visit_float_html, depart_float_html),
               latex=(visit_float_latex, depart_float_latex))
  app.add_enumerable_node(figure_container, 'figure',
               html=(visit_float_html, depart_float_html),
               latex=(visit_float_latex, depart_float_latex))
  app.add_enumerable_node(code_container, 'code-block',
               html=(visit_float_html, depart_float_html),
               latex=(visit_float_latex, depart_float_latex))
