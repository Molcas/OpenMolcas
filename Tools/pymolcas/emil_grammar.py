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
# Copyright (C) 2015-2017, Ignacio Fdez. Galván                        *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)

from re import sub
from pyparsing import *

# helpful methods (python 2/3 compatibility)
def chomp(s):
  return s[:-1] if s.endswith('\n') else s

def chompAction(s, l, t):
  try:
    return list(map(lambda s: chomp(unicode(s)), t))
  except NameError:
    return list(map(chomp, t))

def removeEMILEnd(s):
  return sub(r'\s*<*\s*$','',s)

def removeEMILEndAction(s, l, t):
  try:
    return list(map(lambda s: removeEMILEnd(unicode(s)), t))
  except NameError:
    return list(map(removeEMILEnd, t))

################################################################################
# Attempt to specify a grammar for EMIL (extended Molcas input language)
################################################################################

# do not skip newlines when parsing
ParserElement.setDefaultWhitespaceChars(' \t')

# define a class that grabs the rest of the line, dropping the newline character
NL = Suppress( LineEnd() )
restOfLineNL = restOfLine + NL

# empty lines (between programs/commands)
EmptyLines = OneOrMore( NL )

# the same, but removing the optional "<" and spaces at the end of the line
restOfEMILLine = restOfLine.setParseAction(removeEMILEndAction)
restOfEMILLineNL = restOfEMILLine + NL

# and another that stops at a semicolon
RestComment = Suppress( Literal('//') + restOfLineNL )
restOfLineSemiColon = SkipTo( ( NL | RestComment ) | Suppress(Literal(';')), include=True )

# this seems to catch lines starting with spaces more reliably
#LineStartWS = Suppress( LineStart().leaveWhitespace() )
LineStartWS = Suppress( ZeroOrMore(White(' ')) )

# comments (line and block, and trailing newlines) are removed
LineComment = ( Literal('*') | Literal('//') ) + restOfLineNL
BlockComment = Literal('/*') + SkipTo( StringEnd() | Literal('*/'), include=True ) + Optional(NL)
TrailingNL = OneOrMore( NL ) + StringEnd()
Comment = Suppress( TrailingNL | BlockComment | LineComment )

# convenience alias
AnyLine = restOfLineSemiColon

# one-line EMIL commands start with one or more '>' and extend to the end of the line
EMILMark = (OneOrMore('>') + Optional(White())).setParseAction(replaceWith('>'))
EMILCommand = LineStartWS + Group(EMILMark + restOfEMILLineNL)

# the FILE command requires a closing line
EMILEOF = EMILMark + CaselessLiteral('EOF') + restOfEMILLineNL
EMILFileMark = ( originalTextFor( CaselessLiteral('FILE') + restOfLine ) + NL ).setParseAction(removeEMILEndAction)
EMILFile = LineStartWS + Group(EMILMark + EMILFileMark + SkipTo(EMILEOF).leaveWhitespace().setParseAction(chompAction)) + Suppress(EMILEOF)

# at the moment only FILE is a block command
EMILBlock = EMILFile

EMIL = EMILBlock | EMILCommand

# modules start with their name (preceded with '&') and contain every line
# until another module starts or an EMIL command is found
ModuleName = LineStartWS + Literal('&') + Word( alphanums + '_', min=2 ) + Suppress(restOfLineSemiColon)
ModuleEnd = StringEnd() | EMIL | ModuleName
Module = Group( ModuleName + ZeroOrMore( ~ModuleEnd + ( Comment | AnyLine ) ) )

# the input file should contain only modules, commands and comments
EMIL_Grammar = ZeroOrMore( EMIL | Module | Comment | EmptyLines )

################################################################################
