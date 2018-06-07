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
# Copyright (C) 2015-2018, Ignacio Fdez. Galván                        *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)
from six import text_type

from os import environ, getcwd, symlink, access, F_OK, X_OK, pathsep, defpath, curdir
from os.path import join, split, isfile, isdir, exists, expanduser, realpath, dirname, normcase
from sys import platform, executable
from tempfile import mkdtemp
from shutil import rmtree
from re import sub, match
from subprocess import Popen, PIPE
from json import loads

def set_utf8(var, val, dummy_val='UNKNOWN_VARIABLE'):
  '''Auxiliary function to set environment variables
     making sure they are utf-8-encoded (works with python 2 and 3).
  '''
  try:
    environ[var] = val
  except TypeError:
    environ[var] = text_type(val)
  except UnicodeEncodeError:
    environ[var] = val.encode('utf-8')
  if (val == dummy_val):
    del environ[var]

def get_utf8(var, default=None, dummy_val='UNKNOWN_VARIABLE'):
  '''Auxiliary function to get environment variables
     making sure they are utf-8-encoded (works with python 2 and 3).
  '''
  try:
    val = environ.get(var, default).decode('utf-8')
  except AttributeError:
    val = environ.get(var, default)
  if (val == dummy_val):
    return default
  else:
    return val

def expandvars(string, default=None, skip_escaped=False):
  '''Expand environment variables of form $var and ${var}.
     If parameter 'skip_escaped' is True, all escaped variable references
     (i.e. preceded by backslashes) are skipped.
     Unknown variables are set to 'default'.
     If 'default' is None, they are not expanded.
  '''
  def replace_var(m):
    return environ.get(m.group(2) or m.group(1), m.group(0) if default is None else default)
  reVar = (r'(?<!\\)' if skip_escaped else '') + r'\$(\w+|\{([^}]*)\})'
  return sub(reVar, replace_var, string)

def utf8_open(*args, **kwargs):
  '''Open files using io.open, with utf-8 encoding by default (not in binary mode)'''
  import io
  if ((len(args) > 1) and ('b' in args[1])):
    return io.open(*args, **kwargs)
  else:
    return io.open(*args, encoding='utf-8', **kwargs)

def dotmolcas(filename):
  return join(expanduser('~/.Molcas'), filename)

def find_molcas(xbin_list=None, here=True):
  '''Find a molcas installation and define MOLCAS
     based on a tag file or the path in some configuration file.
     Also add the molcas libraries to LD_LIBRARY_PATH
  '''
  # walk up the tree to find .molcashome
  # (this overwrites MOLCAS if defined)
  if (here):
    path = getcwd()
    (head, tail) = split(path)
    if (tail == ''):
      (head, tail) = split(head)
    while (tail != ''):
      if (isfile(join(path, '.molcashome'))):
        set_utf8('MOLCAS', path)
        break
      path = head
      (head, tail) = split(path)

  # look up pre-defined locations
  pathfiles = [dotmolcas('molcas_default'), dotmolcas('molcas'), join('/etc', 'alternatives', 'molcas')]
  if (get_utf8('MOLCAS', default='') == ''):
    for fn in pathfiles:
      if (isfile(fn)):
        with utf8_open(fn, 'r') as f:
          set_utf8('MOLCAS', f.read().strip())
        break

  path = get_utf8('MOLCAS', default='')

  # set up LD_LIBRARY_PATH
  if (path != ''):
    llp = get_utf8('LD_LIBRARY_PATH', default='')
    lp = join(path, 'lib')
    if (llp == ''):
      set_utf8('LD_LIBRARY_PATH', lp)
    else:
      set_utf8('LD_LIBRARY_PATH', '{0}:{1}'.format(lp, llp))

  # process custom environment commands
  filename = dotmolcas('molcas.shell')
  if isfile(filename):
    source = '. {0}'.format(filename)
    dump = '{0} -c "import os, json;print(json.dumps(dict(os.environ)))"'.format(executable)
    pipe = Popen(['/bin/sh', '-c', '{0} && {1}'.format(source, dump)], stdout=PIPE)
    try:
      env = loads(pipe.stdout.read().decode('utf-8'))
      environ.update(env)
    except:
      print('Warning: {0} could not be processed\n'.format(filename))

  # read xbin.cfg for custom executables
  # (no expansion here)
  if (xbin_list is not None):
    filename = dotmolcas('xbin.cfg')
    if isfile(filename):
      with utf8_open(filename, 'r') as xf:
        for line in xf:
          a, b = line.split('=')
          xbin_list[a.strip()] = b.strip()

def find_sources():
  '''Find the source directories from a build directory,
     if using CMake, and sets the *_SOURCE variables.
     Otherwise they default to $MOLCAS
  '''
  MOLCAS = get_utf8('MOLCAS')
  cmake_src = ''
  OPENMOLCAS_SOURCE = ''
  MOLCAS_SOURCE = ''
  try:
    with utf8_open(join(MOLCAS, 'CMakeCache.txt'), 'r') as f:
      for line in f:
        m = match(r'Molcas_SOURCE_DIR:.*?=(.*)', line)
        if (m):
          cmake_src = m.group(1)
        m = match(r'OPENMOLCAS_DIR:.*?=(.*)', line)
        if (m):
          OPENMOLCAS_SOURCE = m.group(1)
        m = match(r'EXTRA:.*?=(.*)', line)
        if (m):
          MOLCAS_SOURCE = m.group(1)
  except:
    pass
  if (cmake_src == ''):
    cmake_src = MOLCAS
  if (MOLCAS_SOURCE == ''):
    MOLCAS_SOURCE = cmake_src
  if (OPENMOLCAS_SOURCE == ''):
    OPENMOLCAS_SOURCE = cmake_src
  if (get_utf8('MOLCAS_SOURCE', default='') == ''):
    set_utf8('MOLCAS_SOURCE', MOLCAS_SOURCE)
  if (get_utf8('OPENMOLCAS_SOURCE', default='') == ''):
    set_utf8('OPENMOLCAS_SOURCE', OPENMOLCAS_SOURCE)

def attach_streams(output, error, buffer_size=-1):
  '''Attach output and error streams to files, with optional buffer size'''
  import sys
  if (output is not None):
    sys.stdout = utf8_open(output, 'w', buffer_size)
  if (error is not None):
    if (error == output):
      sys.stderr = sys.stdout
    else:
      sys.stderr = utf8_open(error, 'w', buffer_size)

def which(cmd, mode=F_OK | X_OK, path=None):
  '''Copy of shutil.which, since it is not available for python < 3.3
     Modified to return realpath'''
  def _access_check(fn, mode):
    return (exists(fn) and access(fn, mode) and not isdir(fn))
  if dirname(cmd):
    if _access_check(cmd, mode):
      return realpath(cmd)
    return None
  if path is None:
    path = environ.get('PATH', defpath)
  if not path:
    return None
  path = path.split(pathsep)
  if platform == 'win32':
    if curdir not in path:
      path.insert(0, curdir)
    pathext = environ.get('PATHEXT', '').split(pathsep)
    if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
      files = [cmd]
    else:
      files = [cmd + ext for ext in pathext]
  else:
    files = [cmd]
  seen = set()
  for dir in path:
    normdir = normcase(dir)
    if normdir not in seen:
      seen.add(normdir)
      for thefile in files:
        name = join(dir, thefile)
        if _access_check(name, mode):
          return realpath(name)
  return None
