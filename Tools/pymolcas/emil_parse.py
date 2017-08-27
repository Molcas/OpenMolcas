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

from pyparsing import ParseException, __version__ as pyparsing_version
from emil_grammar import EMIL_Grammar
from abstract_flow import *
from re import match, search, IGNORECASE, compile as re_compile

# set of precompiled regular expressions
re_export = re_compile(r'export\s+(\S*)\s*=\s*(.*?)\s*$', flags=IGNORECASE)
re_eval = re_compile(r'eval\s+(\S*)\s*=\s*(.*?)\s*$', flags=IGNORECASE)
re_exit = re_compile(r'exit(\s+(\d+|\$\w*)?)?\s*$', flags=IGNORECASE)
re_shell = re_compile(r'(?:shell|unix)\s+(\S.*?)\s*$', flags=IGNORECASE)
re_exec = re_compile(r'exec\s+(\S.*?)\s*$', flags=IGNORECASE)
re_include = re_compile(r'include\s+(\S.*?)\s*$', flags=IGNORECASE)
re_dowhile = re_compile(r'do\s*while\s*$', flags=IGNORECASE)
re_enddo = re_compile(r'end\s*do\s*$', flags=IGNORECASE)
re_foreach = re_compile(r'foreach\s+(\w+)\s+in\s+\(\s*(.*?)\s*\)\s*$', flags=IGNORECASE)
re_endforeach = re_compile(r'end\s*foreach\s*$', flags=IGNORECASE)
re_dogeo = re_compile(r'do\s*geo\s*$', flags=IGNORECASE)
re_if = re_compile(r'if\s*\(\s*(.*?)\s*(=+|!=|ne|-file)\s*(\S.*?)\s*\)\s*$', flags=IGNORECASE)
re_ifgoto = re_compile(r'if\s*\(\s*(.*?)\s*(=+|!=|ne|-file)\s*(\S.*?)\s*\)\s*go\s*to\s+(\S+)\s*$', flags=IGNORECASE)
re_endif = re_compile(r'end\s*if\s*$', flags=IGNORECASE)
re_echo = re_compile(r'echo\s+(on|off)\s*$', flags=IGNORECASE)
re_file = re_compile(r'file\s+(\S.*?)\s*$', flags=IGNORECASE)
re_rm = re_compile(r'rm\s+(-?force\s+)?(\S.*?)\s*$', flags=IGNORECASE)
re_copy = re_compile(r'(copy|save|clone|collect)\s+(-?force\s+)?(\S*)\s+(\S*)\s*$', flags=IGNORECASE)
re_label = re_compile(r'label\s+(\S+)\s*$', flags=IGNORECASE)
re_goto = re_compile(r'go\s*to\s+(\S+)\s*$', flags=IGNORECASE)

class EMILException(Exception):
  pass

# TODO: check for "ln" once (if) LINK is supported
def check_commandline(item, commandline, forbidden=None):
  name = match(r'\S*', item).group(0).upper()
  if (forbidden is not None):
    if (search(forbidden, commandline)):
      raise EMILException('Forbidden characters in "{0}" command'.format(name))
  if (match(r'cp\s', commandline)):
    raise EMILException('Use COPY instead of "{0} cp"'.format(name))
  if (match(r'rm\s', commandline)):
    raise EMILException('Use RM instead of "{0} rm"'.format(name))
  if (match(r'export\s', commandline)):
    raise EMILException('Use EXPORT instead of "{0} export"'.format(name))

################################################################################
# Parse the result of an initial EMIL_Grammar run into abstract flow tokens
################################################################################

def EMIL_Parse(input_file):
  '''Parse an input file in EMIL language, returns a tuple where the first
     element is Group of abstract flow tokens, and the second element is
     a possibly empty dictionary of files embedded files (filename: content)
  '''

  # Open the file first, to ensure it is utf8-encoded
  with utf8_open(input_file, 'r') as f_in:
    try:
      items = EMIL_Grammar.parseFile(f_in, parseAll=True).asList()
    except ParseException as e:
      if (pyparsing_version >= '2.0.2'):
        raise EMILException('Parsing error: {0}\n{1}'.format(e, e.markInputline()))
      else:
        raise EMILException('Parsing error: {0}'.format(e))

  blocks = [[]]
  files = {}
  level = 0
  leveltype = [['main']]

  # Process each parsed item, adding the corresponding "Statement" subclass
  # to the token list
  for item in items:

    if (item[0] == '&'):
      blocks[level].append(Program(item))

    elif (item[0] == '>'):

      # TODO: variant/option for binary files?
      # >>> FILE
      # this is not added to the token list, but returned as a separate dict
      # so the wrapper can create the files before running the tokens
      if (re_file.match(item[1])):
        re_match = re_file.match(item[1])
        filename = re_match.group(1)
        files[filename]=item[2]

      # >>> EXIT
      elif (re_exit.match(item[1])):
        re_match = re_exit.match(item[1])
        val = re_match.group(2)
        blocks[level].append(Break(val))

      # >>> EXPORT
      elif (re_export.match(item[1])):
        re_match = re_export.match(item[1])
        blocks[level].append(Assignment(re_match.group(1), re_match.group(2)))

      # >>> EVAL
      elif (re_eval.match(item[1])):
        re_match = re_eval.match(item[1])
        blocks[level].append(Assignment(re_match.group(1), re_match.group(2), literal=False))

      # >>> UNIX
      # >>> SHELL
      elif (re_shell.match(item[1])):
        re_match = re_shell.match(item[1])
        check_commandline(item[1], re_match.group(1))
        blocks[level].append(System(re_match.group(1)))

      # >>> EXEC
      elif (re_exec.match(item[1])):
        re_match = re_exec.match(item[1])
        check_commandline(item[1], re_match.group(1), r'[&`;>\\]')
        blocks[level].append(System(re_match.group(1), parallel=True))

      # TODO: variant for including code inside a module
      # >>> INCLUDE
      elif (re_include.match(item[1])):
        re_match = re_include.match(item[1])
        blocks[level].append(Include(re_match.group(1)))

      # >>> DO WHILE
      elif (re_dowhile.match(item[1])):
        blocks.append([])
        level += 1
        leveltype.append(['dowhile'])
      # >>> END DO
      elif (re_enddo.match(item[1])):
        if ((leveltype[level][0] != 'dowhile') and (leveltype[level][0] != 'foreach')):
          raise EMILException('Grouping error')
        # Do not allow empty DO loops
        if (len(blocks[level]) > 0):
          if (leveltype[level][0] == 'foreach'):
            var = leveltype[level][1]
            val = leveltype[level][2]
          else:
            var = None
            val = None
          for item in blocks[level]:
            item.level = level
          blocks[level-1].append(Group(blocks[level], rerun=True, var=var, val=val))
        level -= 1
        blocks.pop()
        leveltype.pop()

      # >>> IF
      # >>> IF () GOTO
      elif (re_if.match(item[1]) or re_ifgoto.match(item[1])):
        re_match = re_if.match(item[1]) or re_ifgoto.match(item[1])
        var = re_match.group(1)
        op = re_match.group(2).lower()
        val = re_match.group(3)
        name = None
        if (re_match.lastindex > 3):
          name = re_match.group(4)
        if (var.lower() == 'iter'):
          var = '$MOLCAS_ITER'
        if (op == '-file'):
          if (var != ''):
            raise EMILException('Error parsing IF command')
          test = Expression(var, 'file', val)
        elif ((op == 'ne') or (op == '!=')):
          test = Expression(var, 'ne', val)
        else:
          test = Expression(var, 'eq', val)
        if (name):
          item = Jump(name)
          item.level = level+1
          blocks[level].append(Group([item], condition=test))
        else:
          blocks.append([])
          level += 1
          leveltype.append(['if', test])
      # >>> END IF
      elif (re_endif.match(item[1])):
        if (leveltype[level][0] != 'if'):
          raise EMILException('Grouping error')
        if (len(blocks[level]) > 0):
          test = leveltype[level][1]
          for item in blocks[level]:
            item.level = level
          blocks[level-1].append(Group(blocks[level], condition=test))
        level -= 1
        blocks.pop()
        leveltype.pop()

      # >>> ECHO (undocumented)
      elif (re_echo.match(item[1])):
        re_match = re_echo.match(item[1])
        blocks[level].append(Setting('echo', re_match.group(1).lower()))

      # >>> RM
      elif (re_rm.match(item[1])):
        re_match = re_rm.match(item[1])
        blocks[level].append(ParTask('rm', (re_match.group(1) is not None), [re_match.group(2)]))
        
      # >>> COPY, SAVE, CLONE, COLLECT
      elif (re_copy.match(item[1])):
        re_match = re_copy.match(item[1])
        blocks[level].append(ParTask(re_match.group(1).lower(), (re_match.group(2) is not None), [re_match.group(3), re_match.group(4)]))
        
      # >>> FOREACH
      elif (re_foreach.match(item[1])):
        re_match = re_foreach.match(item[1])
        var, val = re_match.groups()
        blocks.append([])
        level += 1
        leveltype.append(['foreach', var, val])
      # >>> END FOREACH (for consistency)
      elif (re_endforeach.match(item[1])):
        if (leveltype[level][0] != 'foreach'):
          raise EMILException('Grouping error')
        # Allow empty foreach loops
        #if (len(blocks[level]) > 0):
        if (True):
          var = leveltype[level][1]
          val = leveltype[level][2]
          for item in blocks[level]:
            item.level = level
          blocks[level-1].append(Group(blocks[level], rerun=True, var=var, val=val))
        level -= 1
        blocks.pop()
        leveltype.pop()

      # >>> LABEL
      elif (re_label.match(item[1])):
        re_match = re_label.match(item[1])
        blocks[level].append(Label(re_match.group(1)))

      # >>> GOTO (extension)
      elif (re_goto.match(item[1])):
        re_match = re_goto.match(item[1])
        blocks[level].append(Jump(re_match.group(1)))

      # TODO:
      # >>> LINK
      # >>> DO GEO
      # >>> VERBATIM
      # >>> ELSE (extension)
      # >>> IF ... > / IF ... < (extension)

      # >>> DO GEO unsupported for now
      elif (re_dogeo.match(item[1])):
        raise EMILException('"DO GEO" is unsupported')

      else:
        raise EMILException('Unknown or unsupported statement:\n  > {0}'.format(item[1]))

  if (level != 0):
    raise EMILException('Grouping error')

  return Group(blocks[0]), files
