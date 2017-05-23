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
# Copyright (C) 2015,2016, Ignacio Fdez. Galván                        *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)

from abstract_flow import *

def Python_Parse(filename):
  '''Simple method to return an "abstract flow" object with python code'''
  return Python(filename), {}
