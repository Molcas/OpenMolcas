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

from docutils.transforms import Transform
from docutils import nodes

# Change the "program" role to uppercase
#
class UppercaseProgram(Transform):
  default_priority = 10

  def apply(self):
    for node in self.document.traverse(nodes.inline):
      if ('program' in node['classes']):
        for child in node.traverse(nodes.Text):
          node.replace(child, nodes.Text(child.astext().upper()))

# Setup
#
def setup(app):
  app.add_transform(UppercaseProgram)
