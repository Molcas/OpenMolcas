#!/usr/bin/env python3
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
import sys
sys.dont_write_bytecode = True

warning = ''
stamp = 'd41d8cd98f00b204e9800998ecf8427e'

def main(my_name):
  if sys.hexversion < 0x03040000:
    sys.exit("Python 3.4 or newer is required to run this program.")

  import os
  import os.path
  import codecs
  import argparse
  import hashlib

  # Default signal handlers
  from signal import signal, SIGPIPE, SIG_DFL
  signal(SIGPIPE, SIG_DFL)

  # Unbuffered output, and utf8 encoding
  # TODO: use binary streams throughout to get rid of encoding troubles?
  sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
  sys.stderr = os.fdopen(sys.stderr.fileno(), 'wb', 0)
  try:
    sys.stdout = codecs.getwriter('utf8')(sys.stdout.buffer)
    sys.stderr = codecs.getwriter('utf8')(sys.stderr.buffer)
  except AttributeError:
    sys.stdout = codecs.getwriter('utf8')(sys.stdout)
    sys.stderr = codecs.getwriter('utf8')(sys.stderr)

  with open(my_name, 'rb') as thisfile:
    stamp = hashlib.md5(thisfile.read()).hexdigest()

  if ('PYMOLCAS_STARTED' not in os.environ):
    if (warning):
      print(warning)
    os.environ['PYMOLCAS_STARTED'] = ''

  # Define this as the driver, so scripts can use the same driver
  os.environ['MOLCAS_DRIVER'] = sys.argv[0]

  # Command-line arguments
  parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=42,width=120))
  parser.add_argument('-env', '--environment', help='display information about environment', action='store_true')
  parser.add_argument('-clean', '--clean_scratch', help='clean scratch area after calculation', action='store_true')
  parser.add_argument('-new', '--new_scratch', help='clean scratch area before calculation', action='store_true')
  parser.add_argument('-old', '--old_scratch', help='reuse scratch area (default)', action='store_true')
  parser.add_argument('-ign', '--ignore_environment', help='run ignoring resource files', action='store_true')
  parser.add_argument('-np', '--nprocs', help='number of parallel (MPI) processes', type=int)
  parser.add_argument('-v', '--version', help='print version of the driver', action='store_true')
  parser.add_argument('-o', '--output', help='redirect output stream to FILE', metavar='FILE')
  parser.add_argument('-e', '--error', help='redirect error stream to FILE', metavar='FILE')
  parser.add_argument('-oe', '-eo', '--outputerror', help='redirect output and error streams to FILE', metavar='FILE')
  parser.add_argument('-f', '--files', help='redirect output/error streams to standard files (*.log and *.err)', action='store_true')
  parser.add_argument('-b', '--buffer', help='buffer size for input/error files (-1: default, 0: none, 1: line, n: n bytes)', type=int, default=-1)
  parser.add_argument('--banner', help='print banner (includes version and compilation info)', action='store_true')
  parser.add_argument('--not-here', help='do not try to find MOLCAS in the current directory (and parents)', action='store_true')
  parser.add_argument('-l', '--license', help=argparse.SUPPRESS or 'check if there is a valid license', action='store_true')
  parser.add_argument('filename', help='input file or script', nargs='?', metavar='input_file | script')
  parser.add_argument('extra', help='additional arguments passed to the script', nargs=argparse.REMAINDER, metavar='...')
  parser.usage = '{0} [options] [input_file | script ...]'.format(parser.prog)
  args = vars(parser.parse_args())

  from molcas_aux import find_molcas, find_sources, attach_streams
  from molcas_wrapper import Molcas_wrapper, MolcasException

  # Checking for version right at the beginning, in case MOLCAS cannot be found
  if (args['version']):
    print('python driver version = {0}'.format(Molcas_wrapper.version))
    print('(after the original perl EMIL interpreter of Valera Veryazov)')
    sys.exit(0)

  xbin_list={}
  find_molcas(xbin_list, here=(not args['not_here']))
  find_sources()

  # If this a program defined in xbin, call it
  import subprocess
  if (args['filename'] in xbin_list):
    target = xbin_list[args['filename']]
    print('"{0}" aliased to {1}'.format(args['filename'], target))
    command = '{0} {1}'.format(args['filename'], args['extra'])
    command = [target] + args['extra']
    sys.exit(subprocess.call(command))

  # If this is not calling a program in sbin, pass the extra options to the main program
  # Unfortunately we cannot use Molcas.in_sbin yet, because Molcas is not initialized
  if (args['filename']):
    # Specifically, verify should work with the default environment
    if (args['filename'] == 'verify'):
      args['ignore_environment'] = True
    if (args['extra']):
      in_sbin = False
      for path in ['MOLCAS', 'OPENMOLCAS_SOURCE', 'MOLCAS_SOURCE']:
        filetest = os.path.join(os.environ[path], 'sbin', args['filename'])
        if (os.path.isfile(filetest) and os.access(filetest, os.X_OK)):
          in_sbin = True
      if (not in_sbin):
        for k,v in vars(parser.parse_args(args['extra'])).items():
          if (v):
            args[k] = v

  if args['outputerror']:
    args['output'] = args['outputerror']
    args['error'] = args['outputerror']

  # Define environment variables according to input arguments
  if (args['nprocs']):
    os.environ['MOLCAS_NPROCS'] = str(args['nprocs'])

  if (args['clean_scratch']):
    os.environ['MOLCAS_KEEP_WORKDIR'] = 'NO'

  if (args['new_scratch']):
    os.environ['MOLCAS_NEW_WORKDIR'] = 'YES'

  if (args['old_scratch']):
    os.environ['MOLCAS_NEW_WORKDIR'] = 'NO'

  # Initialize the system
  if (args['files'] and args['filename']):
    fn = os.path.splitext(args['filename'])[0]
    if (not args['output']):
      args['output'] = '{0}.log'.format(fn)
    if (not args['error']):
      args['error'] = '{0}.err'.format(fn)

  try:
    Molcas = Molcas_wrapper(warning=warning, stamp=stamp)
  except MolcasException as message:
    print(message, file=sys.stderr)
    sys.exit(1)

  if (not args['ignore_environment']):
    Molcas.read_environment()

  if (args['environment']):
    Molcas.show_environment()
    sys.exit(0)

  if (args['banner']):
    Molcas.print_banner()
    sys.exit(0)

  if (args['license']):
    rc = Molcas.check_license()
    if (rc is None):
      print('molcas.exe not found')
    elif (rc):
      print('The license is expired, missing or not valid')
    else:
      print('The license is valid')
      if (Molcas.licensee):
        print('Licensed to {0}'.format(Molcas.licensee))
    sys.exit(0)

  if (not args['filename']):
    parser.description = 'MOLCAS has been found at {0}'.format(Molcas.molcas)
    parser.print_help()
    sys.exit(0)
  else:
    # could this be a program from sbin?
    if (Molcas.in_sbin(args['filename'])):
      # make sure calls to "molcas" use this file
      try:
        rc = Molcas.run_sbin([args['filename']] + args['extra'])
        if (rc is not None):
          sys.exit(rc)
      except KeyboardInterrupt:
        sys.exit(1)
    # if not, it must be an input file
    try:
      attach_streams(output=args['output'], error=args['error'], buffer_size=args['buffer'])
    except IOError as e:
      print(e, file=sys.stderr)
      sys.exit(1)
    try:
      Molcas.read_input(args['filename'])
      Molcas.auto()
    except MolcasException as message:
      print(message, file=sys.stderr)
      # Skip verification with unsupported features
      # TODO: remove when ready
      if ('is unsupported' in str(message) or
          'Unknown module' in str(message)):
        Molcas.rc = '_RC_NOT_AVAILABLE_'
      else:
        Molcas.rc = '_RC_JOB_KILLED_'
    except KeyboardInterrupt:
      for f in [sys.stdout, sys.stderr]:
        print('\nInterrupted by user\n', file=f)
      Molcas.rc = '_RC_JOB_KILLED_'

  return(Molcas.rc_num)

# Run pymolcas if not called via import
if (__name__ == '__main__'):
  import shutil
  sys.exit(main(shutil.which(sys.argv[0])))
