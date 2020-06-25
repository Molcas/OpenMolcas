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
# Copyright (C) 2018, Ignacio Fdez. Galván                             *
#***********************************************************************

import os
import errno

def read_answer(answers, retry):
  while True:
    try:
      ans = raw_input()
    except NameError:
      ans = input()
    ans = ans.lower()
    if (ans == ''):
      return ans
    elif (ans in answers):
      return answers[ans]
    else:
      print(retry)

def read_number():
  while True:
    try:
      ans = raw_input()
    except NameError:
      ans = input()
    if (ans == ''):
      return ans
    try:
      ans = int(ans)
    except:
      ans = -1
    if (ans > 0):
      return ans
    else:
      print('Please enter a positive number')

def read_string():
  while True:
    try:
      ans = raw_input()
    except NameError:
      ans = input()
    if (' ' not in ans):
      return ans
    else:
      print('It is a bad idea to use spaces in the path')

def write_molcasrc(molcasrc, program='molcas'):
  config = []

  print("\n"
  "This function will help you create a custom configuration file (molcasrc) for\n"
  "OpenMolcas. The file contains settings for some useful environment variables,\n"
  "you can always modify the file, or override any environment variable for\n"
  "specific calculations.")

  print("\n"
  "MOLCAS_PROJECT\n"
  "--------------\n"
  "Specifies how the project name is set if the \"Project\" variable is not defined.\n"
  "\n"
  "Enter a number, or press enter to keep the default (1):\n"
  "1) NAME:    Use the input filename (minus extension) as project name.\n"
  "2) NAMEPID: Same as 1, but also use a unique name for the scratch directory.\n"
  "3) TIME:    Use the current time as a project name.")
  ans = read_answer({'1': 'NAME', '2': 'NAMEPID', '3': 'TIME'}, 'Please enter 1, 2 or 3')
  if (ans != ''):
    config.append('MOLCAS_PROJECT={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_MEM\n"
  "----------\n"
  "Defines the approximate maximum memory used by OpenMolcas during a calculation.\n"
  "\n"
  "Enter a number (in MiB), or press enter to keep the installation default:")
  ans = read_number()
  if (ans != ''):
    config.append('MOLCAS_MEM={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_MAXITER\n"
  "--------------\n"
  "Defines the maximum number of iterations for a \"Do While\" loop in the input.\n"
  "\n"
  "Enter a number, or press enter to keep the default (50):")
  ans = read_number()
  if (ans != ''):
    config.append('MOLCAS_MAXITER={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_NEW_DEFAULTS\n"
  "-------------------\n"
  "Specify whether or not to use the new OpenMolcas defaults, in particular,\n"
  "with new defaults enabled:\n"
  "\n"
  "* RICD will be enabled (disable with NOCD)\n"
  "* IPEA shift will be disabled (set previous default with IPEA=0.25)\n"
  "\n"
  "Enter y or n, or press enter to keep the default (n):\n"
  "y) YES: The new defaults will be used.\n"
  "n) NO:  The previous defaults will be kept.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_NEW_DEFAULTS={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_TRAP\n"
  "-----------\n"
  "Specifies whether a calculation should stop when a module reports a failure.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) YES: The calculation stops after a failure.\n"
  "n) NO:  The calculation continues after a failure.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_TRAP={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_WORKDIR\n"
  "--------------\n"
  "Defines the parent directory inside which a scratch directory for each\n"
  "calculation will be created. It can be an absolute path, a relative (with\n"
  "respect to the submit directory) path, or the special value \"PWD\", which is\n"
  "equivalent to using \".\" as relative path.\n"
  "\n"
  "Enter a string, or press enter to keep the system's default ($TMPDIR):")
  ans = read_string()
  if (ans != ''):
    config.append('MOLCAS_WORKDIR={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_NEW_WORKDIR\n"
  "------------------\n"
  "Specifies whether the scratch directory will be cleaned *before* a\n"
  "calculation. This can also be specified with the -new and -old flags\n"
  "for {}.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (n):\n"
  "y) YES: The scratch directory is cleaned before a calculation.\n"
  "n) NO:  The scratch directory is not cleaned before a calculation.".format(program))
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_NEW_WORKDIR={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_KEEP_WORKDIR\n"
  "-------------------\n"
  "Specifies whether the scratch directory will be cleaned *after* a calculation.\n"
  "This can also be overridden with the -clean flag for {}.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) YES: The scratch directory is kept after a calculation.\n"
  "n) NO:  The scratch directory is cleaned after a calculation.".format(program))
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_KEEP_WORKDIR={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_OUTPUT\n"
  "-------------\n"
  "Specifies where to save generated files, like orbitals and molden files.\n"
  "It can be an absolute or relative (with respect to the submit directory) path,\n"
  "or the special values \"NAME\" (a subdirectory will be created with the name of\n"
  "the project) or \"PWD\" (files will be saved in the submit directory).\n"
  "\n"
  "Enter a string, or press enter to keep the default (PWD):")
  ans = read_string()
  if (ans != ''):
    config.append('MOLCAS_OUTPUT={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_SAVE\n"
  "-----------\n"
  "Defines how to alter filenames to avoid overwriting existing files.\n"
  "\n"
  "Enter a number or, press enter to keep the dafault (1):\n"
  "1) REPL: Overwrite existing files.\n"
  "2) ORIG: Rename existing files to add the extension \".orig\".\n"
  "3) INCR: New files will have a extension with incremental numbers.")
  ans = read_answer({'1': 'REPL', '2': 'ORIG', '3': 'INCR'}, 'Please enter 1, 2 or 3')
  if (ans != ''):
    config.append('MOLCAS_SAVE={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_MOLDEN\n"
  "-------------\n"
  "Specifies whether Molden format files with molecular orbitals will be created.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) ON:  Molden orbital files are created.\n"
  "n) OFF: Molden orbital files are not created.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_MOLDEN={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_PRINT\n"
  "------------\n"
  "Defines the default print level for OpenMolcas.\n"
  "\n"
  "Enter a number, or press enter to keep the default (2):\n"
  "0) SILENT\n"
  "1) TERSE\n"
  "2) NORMAL\n"
  "3) VERBOSE\n"
  "4) DEBUG\n"
  "5) INSANE")
  ans = read_answer({'0': 'SILENT', '1': 'TERSE', '2': 'NORMAL', '3': 'VERBOSE', '4': 'DEBUG', '5': 'INSANE'}, 'Please enter 0, 1, 2, 3, 4 or 5')
  if (ans != ''):
    config.append('MOLCAS_PRINT={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_ECHO_INPUT\n"
  "-----------------\n"
  "Determines if the pre-processed input will be included in the output.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) YES: Include the input at the beginning of the output.\n"
  "n) NO:  Do not include the input in the output.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_ECHO_INPUT={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_COLOR\n"
  "------------\n"
  "Specifies whether to use a simple markup for important information in the\n"
  "output. This markup can be displayed by e.g. vim (see the Tools/syntax\n"
  "directory).\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) YES: Simple markup will be added.\n"
  "n) NO:  No markup will be added.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_COLOR={}'.format(ans))
    print('')

  print("\n"
  "MOLCAS_REDUCE_PRT\n"
  "-----------------\n"
  "Specifies whether the print level will be reduced inside a \"Do While\" loop.\n"
  "\n"
  "Enter y or n, or press enter to keep the default (y):\n"
  "y) YES: Print level is reduced inside a \"Do While\" loop.\n"
  "n) NO:  Print level is the same as outside.")
  ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
  if (ans != ''):
    config.append('MOLCAS_REDUCE_PRT={}'.format(ans))
    print('')

  if (config):
    print("\n"
    "Based on your answers, we recommend the following molcasrc file:\n"
    "-------------------------")
    print("\n".join(config))
    print('-------------------------')

    print("\n"
    "The file {} will be created or overwritten.\n"
    "Is this OK? (y/N)".format(molcasrc))
    ans = read_answer({'y': 'YES', 'n': 'NO'}, 'Please answer y or n')
    if (ans == 'YES'):
      try:
        dirrc = os.path.dirname(molcasrc)
        if (not os.path.exists(dirrc)):
          try:
            os.makedirs(dirrc)
          except OSError as e:
            if (e.errno != errno.EEXIST):
              raise
        with open(molcasrc, 'w') as f:
          f.write("# molcasrc configuration file written by {} -setup\n\n".format(program))
          f.write("\n".join(config))
          f.write("\n")
      except:
        raise
        return False
  else:
    print("\n"
    "All defaults were selected, no molcasrc file will be written.")

  return True
