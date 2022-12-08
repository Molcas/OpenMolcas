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
# Copyright (C) 2022, Ignacio Fdez. Galván                             *
#***********************************************************************

import sphinx

# Patch to work around Sphinx bug #9529
if (3, 5, 0, '', 0) <= sphinx.version_info < (4, 5, 0, '', 0):
  from sphinx.writers.latex import CR
  from docutils import nodes
  from docutils.nodes import Element
  from typing import cast
  visit_start_of_file_orig = sphinx.writers.latex.LaTeXTranslator.visit_start_of_file
  visit_footnote_orig = sphinx.writers.latex.LaTeXTranslator.visit_footnote
  def visit_start_of_file_patched(self, node: Element) -> None:
    visit_start_of_file_orig(self, node)
    self.body.append(CR + r'\sphinxstepscope' + CR)
  def visit_footnote_patched(self, node: Element) -> None:
    visit_footnote_orig(self, node)
    if 'auto' in node:
      label = cast(nodes.label, node[0])
      self.body.append(r'\phantomsection\label{\thesphinxscope.%s}%%' % label.astext() + CR)
  sphinx.writers.latex.LaTeXTranslator.visit_start_of_file = visit_start_of_file_patched
  sphinx.writers.latex.LaTeXTranslator.visit_footnote = visit_footnote_patched

def setup(app):
  pass
