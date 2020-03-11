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
# Copyright (C) 2015-2020, Ignacio Fdez. Galván                        *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)
try:
  from builtins import bytes
except ImportError:
  from future.builtins import bytes
try:
  from six import text_type
except ImportError:
  text_type = str

from os import environ, access, W_OK, X_OK, listdir, remove, getpid, getcwd, makedirs, symlink, devnull
from os.path import isfile, isdir, isabs, join, basename, splitext, getmtime, abspath, exists, relpath, realpath
from datetime import datetime
from shutil import copy2, move, rmtree, Error
from subprocess import check_output, STDOUT, CalledProcessError
from re import compile as re_compile, search, sub, MULTILINE, IGNORECASE
from io import BytesIO
from resource import getrusage, RUSAGE_CHILDREN
from glob import glob
from errno import EEXIST, ENOENT
from shlex import split as shsplit
from contextlib import contextmanager
from textwrap import fill

from tee import teed_call
from emil_parse import EMIL_Parse, EMILException
from python_parse import Python_Parse
from molcas_aux import *
from check_test import *
from validate import *

# python2 has no FileNotFoundError, so we will have to check
# the errno attribute of the raised exception
try:
  #python3
  FileNotFoundError
except:
  #python2 (IOError when opening files, OSError when removing files)
  FileNotFoundError = (IOError, OSError)

hcbanner = '''#
      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      &&&                                                      &&&
      &&&                                                      &&&
      &&&    OpenMolcas source code or broken installation     &&&
      &&&                                                      &&&
      &&&                                                      &&&
      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 \n
            version $V$
        patch level $P$
                    $Px$
 \n
                 Copyright (C) The OpenMolcas Authors
           For the author list and the recommended citation,
                    consult the file CONTRIBUTORS.
'''
#===============================================================================

# TODO:
# MOLCAS_Count
# MOLCAS_EMIL_DEBUG
# MOLCAS_IN_GEO
# MOLCAS_ISDEV
# MOLCAS_LASTMOD
# MOLCAS_LINK
# MOLCAS_LOG
# MOLCAS_MOLDEN
# MOLCAS_OLDPWD
# MOLCAS_PARNELL
# MOLCAS_SERIAL
# MOLCAS_SUBMIT_DIR
# MOLCAS_SUBMIT_PWD
# MOLCAS_TASK_INPUT
# MOLCAS_TIME
# MOLCAS_UNIX_SECURE
# MOLCAS_ZOMBIE

class MolcasException(Exception):
  pass

class Molcas_wrapper(object):

  version = 'py2.08'
  rc = 0

  def __init__(self, **kwargs):
    # Get Molcas istallation
    self.molcas = abspath(get_utf8('MOLCAS', default=''))
    if (self.molcas == ''):
      raise MolcasException('"MOLCAS" is not defined')
    if (not isfile(join(self.molcas, '.molcashome'))):
      raise MolcasException('"{0}" is not a valid MOLCAS installation'.format(self.molcas))
    # Get location of sources (other than MOLCAS)
    self.sources = []
    self.sources.append(abspath(get_utf8('OPENMOLCAS_SOURCE', default='')))
    self.sources.append(abspath(get_utf8('MOLCAS_SOURCE', default='')))
    self.sources = [i for i in self.sources if (i != '' and i != self.molcas)]
    # Some settings
    self.allow_shell = True
    self.echo = True
    self.input_filename = ''
    self.only_validate = False
    if ('allow_shell' in kwargs):
      self.allow_shell = bool(kwargs['allow_shell'])
    if ('echo' in kwargs):
      self.echo = bool(kwargs['echo'])
    if ('input_file' in kwargs):
      self.input_filename = kwargs['input_file']
    if ('warning' in kwargs):
      self.warning = kwargs['warning']
    if ('stamp' in kwargs):
      self.stamp = kwargs['stamp']
    if ('validate' in kwargs):
      self.only_validate = kwargs['validate']
    self.licensee = None
    self.keywords = None
    self.rc = None
    self._ready = False
    self._goto = False
    self._parse_rte()
    self._parse_codes()
    self._parse_alias()
    self._parse_keywords()
    self._resources = (0,0,0)
    self._check_count = 0
    self._loop_level = 0
    self._nest_level = 0
    self._last_module = True
    self.parnell = join(self.molcas, 'bin', 'parnell.exe')
    if (self.input_filename):
      self.read_input(self.input_filename)

  @contextmanager
  def Loop(self):
    self.enter_loop()
    try:
      yield
    finally:
      self.exit_loop()

  def _parse_rte(self):
    '''Parse the molcas.rte file'''
    rte = {}
    filename = join(self.molcas, 'molcas.rte')
    try:
      r_f = utf8_open(filename, 'r')
    except IOError:
      pass
    else:
      with r_f:
        comment_match = re_compile(r'#.*$')
        rte_match = re_compile(r'(\w+)\s*=\s*([\'\"])(.*)\2')
        for line in r_f:
          line = comment_match.sub('', line)
          match = rte_match.search(line)
          if (match):
            rte[match.group(1)] = match.group(3)
    self._rte = rte

  def _parse_codes(self):
    '''Parse the list of return codes (in plain or perl syntax)'''
    code_list = {}
    filename = join(self.molcas, 'data', 'rcodes.txt')
    rc_match = re_compile(r'(\w+)\s*=\s*(\d+)')
    try:
      c_f = utf8_open(filename, 'r')
    except IOError:
      filename = join(self.molcas, 'data', 'warnings.plx')
      rc_match = re_compile(r'\$(\w+)\s*=\s*(\d+)\s*;')
      try:
        c_f = utf8_open(filename, 'r')
      except IOError:
        c_f = utf8_open(devnull, 'r')
    with c_f:
      for line in c_f:
        match = rc_match.search(line)
        if (match):
          code_list[int(match.group(2))] = match.group(1)
    self._rc_codes = code_list

  def _parse_alias(self):
    '''Parse the list of aliases (in plain or perl syntax)'''
    self.alias = {}
    filename = join(self.molcas, 'data', 'modalias.txt')
    alias_match = re_compile(r'([^#].*=.*)')
    pair_match = re_compile(r'\s*(\w*)\s*=\s*(\w*)\s*')
    try:
      a_f = utf8_open(filename, 'r')
    except IOError:
      filename = join(self.molcas, 'data', 'alias.plx')
      alias_match = re_compile(r'\s*%NickName\s*=\s*\((.*)\);')
      pair_match = re_compile(r'q\/(.*)\/=>q\/(.*)\/')
      try:
        a_f = utf8_open(filename, 'r')
      except IOError:
        return
    with a_f:
      for line in a_f:
        match = alias_match.match(line)
        if (match):
          pairs = match.group(1).split(',')
          for i in pairs:
            match = pair_match.match(i)
            if (match):
              self.alias[match.group(1)] = match.group(2)

  def _parse_keywords(self):
    try:
      self.keywords = read_db(join(self.molcas, 'data', 'keyword.xml'))
    except NameError:
      self.keywords = 'no_lxml'
    except:
      pass

  def read_environment(self):
    '''Read environment variables from rc files (in priority order)'''
    rcfiles = [dotmolcas('molcasrc'), join(self.molcas, 'molcasrc')]
    var_match = re_compile(r'\s*([^#].*?)\s*=\s*(.*)')
    for rcfile in rcfiles:
      try:
        r_f = utf8_open(rcfile, 'r')
        for line in r_f:
          match = var_match.match(line)
          if (match):
            if (match.group(1) not in environ):
              set_utf8(match.group(1), match.group(2))
        r_f.close()
      except (OSError, IOError):
        pass

  def setup(self):
    if (not self._rte):
      raise MolcasException('"{0}" is not a valid MOLCAS installation'.format(self.molcas))
    if (self._ready):
      print('*******************************')
      print('* ERROR -- setup called twice *')
      print('*******************************')
      self.rc = '_RC_GENERAL_ERROR_'
      return
    # Define some hard-coded/initial environment variables
    set_utf8('LANG', 'C')
    set_utf8('LC_ALL', '')
    set_utf8('EMIL_RC2', '0')
    set_utf8('EMIL_InLoop', '0')
    set_utf8('MOLCAS_STRUCTURE', '0')
    if (self.check_license() == 0):
      set_utf8('MOLCAS_STAMP', '5')
    else:
      set_utf8('MOLCAS_STAMP', 'UNKNOWN_VARIABLE')
    # Get current directory
    self.currdir = getcwd()
    set_utf8('CurrDir', self.currdir)
    set_utf8('MOLCAS_SUBMIT_DIR', self.currdir)
    if ((self.currdir is None) or (self.currdir == '')):
      raise MolcasException('"CurrDir" is not defined')
    # Get project name
    usepid = False
    top_project = get_utf8('MOLCAS_PROJECT', default='NAME')
    if (top_project.upper() == 'NAMEPID'):
      top_project = splitext(self.input_filename)[0]
      usepid = True
    elif (top_project.upper() == 'NAME'):
      top_project = splitext(self.input_filename)[0]
    elif (top_project.upper() == 'TIME'):
      top_project = '{:%Y_%m_%d_%H_%M}'.format(datetime.now())
    self.project = get_utf8('Project', default=top_project)
    self.project = basename(self.project)
    if ((self.project is None) or (self.project == '')):
      self.project = 'Noname'
    set_utf8('Project', self.project)
    # Get scratch directory
    self._tmpdir = get_utf8('TMPDIR', default='/tmp')
    self._pid = get_utf8('MOLCAS_PID', default=getpid())
    self.scratch = get_utf8('WorkDir', default='')
    if (self.scratch == ''):
      top_workdir = get_utf8('MOLCAS_WORKDIR', default=self._tmpdir)
      if (isabs(top_workdir)):
        self.scratch = join(top_workdir, self.project)
        if (usepid):
          self.scratch = '{0}.{1}'.format(self.scratch, self._pid)
      else:
        if (top_workdir.upper() == 'PWD'):
          self.scratch = self.currdir
        else:
          self.scratch = join(self.currdir, top_workdir)
    self.scratch = abspath(self.scratch)
    set_utf8('WorkDir',self.scratch)
    if ((self.scratch is None) or (self.scratch == '')):
      raise MolcasException('"WorkDir" is not defined')
    if (not self.only_validate):
      if (self.parallel_task(['base', self.scratch]) != 0):
        raise MolcasException('parnell failed to create a WorkDir at {0}'.format(self.scratch))
    # Get output directory
    self.output = get_utf8('MOLCAS_OUTPUT', default=get_utf8('PBS_O_WORKDIR', default=self.currdir))
    if (not isabs(self.output)):
      if (self.output.upper() == 'PWD'):
        self.output = self.currdir
      elif (self.output.upper() == 'WORKDIR'):
        self.output = self.scratch
      elif (self.output.upper() == 'NAME'):
        self.output = join(self.currdir, self.project)
      else:
        self.output = join(self.currdir, self.output)
    set_utf8('MOLCAS_OUTPUT', self.output)
    if ((self.output is None) or (self.output == '')):
      raise MolcasException('Output directory is not defined')
    if (not self.only_validate):
      try:
        makedirs(self.output)
      except OSError as e:
        if ((e.errno == EEXIST) and (isdir(self.output))):
          pass
        else:
          raise
      if (not isdir(self.output) or not access(self.output, W_OK)):
        raise MolcasException('"{0}" is not a writable directory'.format(self.output))
    # Get other variables
    if (get_utf8('MOLCAS_MEM', default='') == ''):
      set_utf8('MOLCAS_MEM', self._rte['DEFMOLCASMEM'])
    if (get_utf8('MOLCAS_DISK', default='') == ''):
      set_utf8('MOLCAS_DISK', self._rte['DEFMOLCASDISK'])
    if (get_utf8('MOLCAS_NPROCS', default='') == ''):
      set_utf8('MOLCAS_NPROCS', '1')
    if (get_utf8('SubProject', default='') == ''):
      set_utf8('SubProject', '')
    if (get_utf8('GeoDir', default='') == ''):
      set_utf8('GeoDir', join(self.scratch, self.project+'.GEO'))
    self.save_mode = get_utf8('MOLCAS_SAVE', default='repl').lower()
    if (get_utf8('MOLCAS_NEW_WORKDIR', default='NO').upper() == 'YES'):
      self.delete_scratch(True)
    self._is_empty = self._is_scratch_empty()
    self._ready = True

  def rc_to_name(self, rc):
    try:
      rc = int(rc)
      if (rc in self._rc_codes):
        return self._rc_codes[rc]
      else:
        return rc
    except:
      return rc

  @property
  def rc_num(self):
    try:
      rc = int(self.rc)
    except:
      num = [k for k,v in self._rc_codes.items() if v==self.rc]
      if (num):
        rc = int(num[0])
      else:
        rc = -1
    return rc

  def print_banner(self):
    try:
      with utf8_open(join(self.molcas, 'data', 'banner.txt'), 'r') as banner_file:
        banner = banner_file.read().rstrip('\n')
    except FileNotFoundError as e:
      if (e.errno == ENOENT):
        banner = hcbanner.rstrip('\n')
      else:
        raise
    tag = '(unknown)'
    tag_x = ''
    try:
      with utf8_open(join(self.molcas, '.molcasversion'), 'r') as version_file:
        for line in version_file:
          if (search('\.x\d', line)):
            tag_x = line.rstrip()
          else:
            tag = line.rstrip()
    except FileNotFoundError as e:
      if (e.errno == ENOENT):
        try:
          command = ["git", "describe", "--always", "--match", "v*", "--dirty"]
          line = check_output(command, stderr=STDOUT).decode('utf-8')
          if (search('\.x\d', line)):
            tag_x = line.rstrip()
          else:
            tag = line.rstrip()
        except:
          tag = '(unknown)'
      else:
        raise
    if ((tag_x != '') and (tag == '(unknown)')):
      tag = tag_x
      tag_x = ''
    v_match = re_compile('v(\d+\.\d+)[\.-](.*)')
    match = v_match.match(tag)
    if (match):
      version = match.group(1)
      patch = match.group(2)
    else:
      version = tag
      patch = ''
    match = v_match.match(tag_x)
    if (match):
      patch_x = match.group(2)
    else:
      patch_x = ''
      with utf8_open(join(self.molcas, '.molcashome')) as homefile:
        if (homefile.read().rstrip() not in ['openmolcas', '']):
          patch_x = 'x'
    date = datetime.now()
    # Remove comments and empty lines (not if there are spaces)
    banner = sub(r'#.*$', '', banner, flags=MULTILINE)
    banner = sub(r'^\n', '', banner, flags=MULTILINE)
    # Replace version and patch numbers
    patch_re = re_compile(r'^(.*)\$P\$(.*)$', flags=MULTILINE)
    patch_match = patch_re.search(banner)
    if (patch_match):
      fmt = '{{0}}\n{{1:>{0}}} {{2}}'.format(len(patch_match.group(1))-1)
      s = fmt.format(patch_match.group(0), '', '$Px$')
      if (banner.find('$Px$') < 0):
        banner = patch_re.sub(s, banner)
    banner = sub(r'\$V\$', version, banner)
    banner = sub(r'\$P\$', patch, banner)
    if (patch_x):
      banner = sub(r'\$Px\$', patch_x, banner)
    else:
      banner = sub(r'.*\$Px\$.*\n', patch_x, banner)
    # Christmas tree between Dec. 18th and Jan. 8th
    # (disabled)
    #if ((date.month == 12 and date.day > 17) or
    #    (date.month == 1 and date.day < 9)):
    #  banner = banner.replace('@', '^')
    #  banner = banner.replace('!', 'o')
    #else:
    #  banner = banner.replace('@', ' ')
    #  banner = banner.replace('!', ' ')
    # Do the main character substitution
    rep = get_utf8('MOLCAS_BANNER')
    if (patch_x == ''):
      ch0 = '&'
      ch1 = '+'
      ch2 = '*'
      if (not rep):
        rep = 'OPENMOLCAS'
    else:
      ch0 = '&'
      ch1 = '*'
      ch2 = '+'
      if (not rep):
        rep = 'MOLCAS'
    maxlen = 0
    new_banner = []
    for line in banner.split('\n'):
      maxlen = max(maxlen, len(line))
      lline = list(line.rstrip())
      j = 0
      for i in range(len(lline)):
        if ((lline[i] == ch0) or (lline[i] == ch1)):
          lline[i] = rep[j]
          j = (j+1) % len(rep)
        elif (lline[i] == ch2):
          lline[i] = ' '
      line = ''.join(lline)
      new_banner.append(line.rstrip())
    # Reverse on Apr. 1st
    if (date.month == 4 and date.day == 1):
      fmt = '{{:<{0}}}'.format(maxlen)
      new_banner = [fmt.format(line)[::-1].rstrip() for line in new_banner]
    banner = '\n'.join(new_banner)
    print(banner)
    if hasattr(self, 'stamp'):
      fmt = '{:^72}'
      print(fmt.format('').rstrip())
      print(fmt.format('*************************************************').rstrip())
      print(fmt.format('* pymolcas version {0:<28} *'.format(self.version)).rstrip())
      print(fmt.format('*   build {0}      *'.format(self.stamp)).rstrip())
      print(fmt.format('*   (after the EMIL interpreter by V. Veryazov) *').rstrip())
      print(fmt.format('*************************************************').rstrip())
    try:
      with utf8_open(join(self.molcas, 'data', 'info.txt'), 'r') as info_file:
        print(info_file.read().rstrip('\n'))
    except:
      pass

  def file_exists(self, filename):
    filetest = expandvars(filename, default='UNKNOWN_VARIABLE')
    if (not isabs(filetest)):
      filetest = join(self.scratch, filetest)
    return isfile(filetest)

  def add_files(self, files):
    for fname, content in files.items():
      filename = expandvars(fname, default='UNKNOWN_VARIABLE')
      with utf8_open(join(self.scratch, filename), 'w') as f:
        print(content, file=f)

  def read_input(self, input_file):
    self.input_filename = input_file
    self.setup()
    if (not isfile(input_file)):
      raise MolcasException('File "{0}" does not exist'.format(input_file))
    try:
      self.flow, files = Python_Parse(input_file)
    except SyntaxError:
      try:
        self.flow, files = EMIL_Parse(input_file)
      except EMILException as message:
        raise MolcasException('EMIL: {0}'.format(message))
    if (not self.only_validate):
      self.add_files(files)

  def print_input(self):
    if (get_utf8('MOLCAS_ECHO_INPUT', default='YES').upper() != 'NO'):
      foo = text_type(self.flow)
      print('++ ---------   Input file   ---------\n')
      print(foo)
      print('\n-- ----------------------------------')

  def print_environment(self):
    # These variables must always be defined
    lines = {}
    lines['Project'] = self.project
    lines['Submitted from'] = self.currdir
    lines['Scratch area'] = self.scratch
    if (self.output == self.scratch):
      lines['Save outputs to'] = 'WORKDIR'
    else:
      lines['Save outputs to'] = self.output
    lines['Molcas'] = self.molcas
    keys = ['Project',
            'Submitted from',
            'Scratch area',
            'Save outputs to',
            'Molcas']
    # Get all the other MOLCAS* environment variables
    molcas_env = {k:get_utf8(k) for k in environ.keys() if k.startswith('MOLCAS')}
    keys.extend(sorted(molcas_env.keys()))
    for key in molcas_env.keys():
      lines[key] = molcas_env[key]
    # Do not print the variables in the list
    do_not_print = ['MOLCAS_PROJECT',
                    'MOLCAS_WORKDIR',
                    'MOLCAS',
                    'MOLCAS_STAMP',
                    'MOLCAS_SUBMIT_DIR',
                    'MOLCAS_INFO',
                    'MOLCAS_TIME',
                    'MOLCAS_DISK',
                    'MOLCAS_MEM',
                    'MOLCAS_ITER',
                    'MOLCAS_OLDPWD',
                    'MOLCAS_OUTPUT',
                    'MOLCAS_TEST',
                    'MOLCAS_LOG',
                    'MOLCAS_ECHO_INPUT',
                    'MOLCAS_ISDEV',
                    'MOLCAS_SILENT_LIC']
    for key in do_not_print:
      if key in lines:
        keys.remove(key)
        del lines[key]
    # Format string and length
    maxlen = max(len(x) for x in keys)
    ini = '  |'
    fmt = ini + ' {{0:>{0}}} = {{1}}'.format(maxlen)
    fmt2 = ini + '  {{0:>{0}}}: {{1}}'.format(maxlen)
    maxlen = maxlen + max(len(x) for x in lines.values())
    # Print the environment summary
    print('   '+'-'*(maxlen+4))
    print(ini)
    for key in keys[0:5]:
      print(fmt2.format(key, lines[key]))
    print(ini)
    if (self._is_empty):
      print(ini+' Scratch area is empty')
    else:
      print(ini+' Scratch area is NOT empty')
    print(ini)
    for key in keys[5:]:
      print(fmt.format(key, lines[key]))
    print(ini)
    print('   '+'-'*(maxlen+4))

  def write_environment(self):
    molcas_env = {k:get_utf8(k) for k in environ.keys() if k.startswith('MOLCAS')}
    # Write the molcas.env file
    fmt = '{0}={1}'
    env_file = join(self.scratch, 'molcas.env')
    with utf8_open(env_file, 'w') as env:
      print('# this is a generated file! do not edit.', file=env)
      print(fmt.format('WorkDir', self.scratch), file=env)
      print(fmt.format('Project', self.project), file=env)
      print(fmt.format('SubProject', get_utf8('SubProject', '')), file=env)
      print(fmt.format('GeoDir', get_utf8('GeoDir', '')), file=env)
      for key in sorted(molcas_env.keys()):
        print(fmt.format(key, molcas_env[key]), file=env)
    self.parallel_task(['c', '1', env_file, self.scratch])

  def show_environment(self):
    molcas_env = {k:get_utf8(k) for k in environ.keys() if k.startswith('MOLCAS')}
    keys = sorted(molcas_env.keys())
    text = ''
    line = ''
    if 'MOLCAS_STAMP' in keys:
      keys.remove('MOLCAS_STAMP')
    for key in keys:
      prev_line = line
      line = '{0}={1}'.format(key, molcas_env[key])
      if (len(line) >= 60):
        if (len(prev_line) > 0):
          text += '{0}\n'.format(prev_line)
        text += '{0}\n'.format(line)
        line = ''
      elif (len(prev_line) > 0):
        text += '{0:60} {1}\n'.format(prev_line, line)
        line = ''
    text += line
    print(text)

  def _is_scratch_empty(self):
    # Check if there are any files in $WorkDir, other than those listed
    ignore = re_compile(r'(molcas|paraops|\.input$|\.inp$|\.log$|\.err$)')
    try:
      files = [f for f in listdir(self.scratch) if (not ignore.search(f))]
    except FileNotFoundError as e:
      if (e.errno == ENOENT):
        files = []
    if (len(files) == 0):
      return True
    else:
      return False

  def add_resources(self, res):
    self._resources = tuple(map(sum, zip(self._resources, res)))

  def reset_resources(self):
    self._resources = (0,0,0)

  def enter_loop(self, loop_type=None):
    self._loop_level += 1
    set_utf8('EMIL_InLoop', self._loop_level)
    if (loop_type == 'do'):
      # This simply marks if we ever enter a do loop
      set_utf8('MOLCAS_STRUCTURE', '1')

  def exit_loop(self):
    self._loop_level -= 1
    if (self._loop_level < 0):
      raise MolcasException('Bad nesting level')
    set_utf8('EMIL_InLoop', self._loop_level)

  def startup(self):
    if (self._ready):
      if (not self.only_validate):
        self.write_environment()
      print('   This run of MOLCAS is using the pymolcas driver')
      if (self.licensee):
        print('   Licensed to: {0}'.format(self.licensee))
      self.print_banner()
      print('\n')
      if (hasattr(self, 'warning') and self.warning):
        print(self.warning)
      self.print_environment()
      if (not self.only_validate):
        print('')
        self.print_input()

  def end(self):
    if (hasattr(self, 'rc')):
      # In case of error in parallel, collect the "stdout" files from the slaves,
      # since they may have important information (create an empty "stdout" file
      # to avoid further errors)
      if ((self.rc_to_name(self.rc) == '_RC_INTERNAL_ERROR_') and not self.is_serial):
        f = 'stdout'
        open(join(self.scratch, f), 'a').close()
        self.parallel_task(['c', '2', f, self.currdir], force=True)
      for l in self._final_rc(self.rc_num):
        print(l)
      if (self.rc_num == -1):
        print('\nAborting...')
    if (self._resources != (0,0,0)):
      print('    Timing: Wall={0:.2f} User={1:.2f} System={2:.2f}'.format(*self._resources))
    if (get_utf8('MOLCAS_KEEP_WORKDIR', default='YES').upper() == 'NO'):
      self.delete_scratch()

  def _final_rc(self, rc):
    rc_form = re_compile('rc={0}\s(.*)'.format(rc))
    text = []
    try:
      with utf8_open(join(self.molcas, 'data', 'landing.txt'), 'r') as l:
        for line in l:
          match = rc_form.match(line)
          if match:
            text.append(match.group(1))
    except:
      pass
    if not text:
      text.append('Non-zero return code')
    if text:
      maxlen = max(len(l) for l in text)
      line = '.{}.'.format('#'*(maxlen+4))
      text = ['', line] + ['.# {} #.'.format(l) for l in text] + [line, '']
    return text

  def validate(self):
    import abstract_flow
    if (self.only_validate):
      self.startup()
      print('\nThe input file will be validated against the documented syntax, no')
      print('files will be created or modified. Note that there is no guarantee that')
      print('a validated input will actually run as intended, but a validation error')
      print('will most likely result in a run error.')
    self.rc = 0
    # flatten input to keep only programs (recursively inside groups)
    flat = self.flow.contents[:]
    redo = True
    while (redo):
      redo = False
      new = []
      for item in flat:
        if (isinstance(item, abstract_flow.Program)):
          new.append(item)
        else:
          redo = True
          if (isinstance(item, abstract_flow.Group)):
            new.extend(item.contents)
      flat = new
    final_rc = 0
    for item in (flat):
      inp = '&' + '\n'.join(item.lines[1:])
      inputlines = sanitize(inp).strip().split('\n')
      inputlines.append('End of input')
      rc, result = validate(inputlines, self.keywords)
      print('\n--')
      for line in result:
        print(line)
      if (rc == 0):
        print('\nLooks good!')
      elif (rc > 0):
        print('\n*** FAILED! ***')
      print('--')
      final_rc = max(final_rc, rc)
    if (final_rc):
      print('\n*************************************')
      print('*** There were validation errors. ***')
      print('*************************************')
    else:
      print('\nValidation passed.')
    self.rc = final_rc

  def auto(self):
    self.startup()
    set_utf8('MOLCAS_ITER', '0')
    if (get_utf8('MOLCAS_VALIDATE') == 'FIRST'):
      self.validate()
      if (self.rc > 0):
        self.rc = '_RC_INPUT_ERROR_'
        self.end()
        return
    try:
      self.run_logue('prologue')
      self.rc = self.flow.run(self)
      if (self.rc is None):
        self.rc = '_RC_ALL_IS_WELL_'
      if (self._goto):
        print('Label "{0}" not found'.format(self._goto))
        self.rc = '_RC_INPUT_ERROR_'
      if ((self.rc_to_name(self.rc) == '_RC_ALL_IS_WELL_') and (not self._last_module)):
        print('*** Extra labels in check file')
        self.rc = '_RC_CHECK_ERROR_'
      self.run_logue('epilogue')
    except KeyboardInterrupt:
      self.rc = '_RC_JOB_KILLED_'
    finally:
      self.end()

  @property
  def is_serial(self):
    serial = get_utf8('MOLCAS_SERIAL')
    nprocs = int(get_utf8('MOLCAS_NPROCS', default='1'))
    return ((serial == 'YES') or (nprocs == 1))

  def wrap_command(self, executable, may_be_serial=True, may_debug=True):
    '''
    Wrap the specified command as appropriate.
    If may_be_serial=False, never uses RUNBINARYSER.
    Returns a list.
    '''
    nprocs = int(get_utf8('MOLCAS_NPROCS', default='1'))
    set_utf8('MOLCAS_NPROCS', nprocs)
    if (may_debug):
      debugger = get_utf8('MOLCAS_DEBUGGER')
    else:
      debugger = ''
    if (debugger):
      command = '{0} $program'.format(debugger)
    elif (may_be_serial and self.is_serial):
      command = self._rte['RUNBINARYSER']
    else:
      command = self._rte['RUNBINARY']
    command = sub(r'\$program', executable, command)
    command = shsplit(expandvars(command, default='UNKNOWN_VARIABLE'))
    return command

  def _set_threads(self):
    threads = get_utf8('MOLCAS_THREADS')
    if (threads is not None):
      try:
        threads = int(threads)
        if (threads == 0):
          if ('OMP_NUM_THREADS' in environ):
            del environ['OMP_NUM_THREADS']
        else:
          set_utf8('OMP_NUM_THREADS', threads)
      except:
        pass

  def run_module(self, name, inp):
    if (self._ready):
      set_utf8('MOLCAS_CURRENT_PROGRAM', name.lower())
      self._set_threads()
      self.write_environment()
      try:
        self._current_module = Molcas_module(self, name, inp)
      except:
        self.rc = '_RC_NOT_AVAILABLE_'
        return self.rc
      rc = self._current_module.run()
      rcname = ''
      if (rc is not None):
        self.rc = rc
        rcname = self.rc_to_name(self.rc)
      if (rcname == '_RC_INVOKED_OTHER_MODULE_'):
        self._nest_level += 1
        set_utf8('EMIL_RC2',self._nest_level)
        stdin_name = 'Stdin.{0}.{1}'.format(self._nest_level, self._loop_level)
        Stdin = join(self.scratch, stdin_name)
        try:
          new_inp, files = EMIL_Parse(Stdin)
        except EMILException as message:
          raise MolcasException('EMIL: {0}'.format(message))
        self.add_files(files)
        rc = new_inp.run(self)
        if (rc is None):
          self.rc = '_RC_ALL_IS_WELL_'
        else:
          self.rc = rc
          rcname = self.rc_to_name(self.rc)
        self._nest_level -= 1
        set_utf8('EMIL_RC2',self._nest_level)
      # I wouldn't put this here, but compatibility obliges
      if ((rcname in ['_RC_ALL_IS_WELL_', '_RC_CONTINUE_LOOP_', '_RC_CONTINUE_UNIX_LOOP_', '_RC_EXIT_EXPECTED_', '_RC_INVOKED_OTHER_MODULE_']) or \
          (get_utf8('MOLCAS_TRAP', default='on').lower() == 'off')):
        self._run_check()
    else:
      print('***************************************************')
      print('*** MOLCAS system not ready -- setup not called ***')
      print('***************************************************')
      self.rc = '_RC_GENERAL_ERROR_'
    return self.rc

  def _run_check(self):
    '''
    Run the check command for verification.
    Give preference to the compiled "CHECK" program, but use a python method if not available.
    '''
    check = get_utf8('MOLCAS_TEST', default='').lower()
    if ((check == 'check') or (check=='gene')):
      set_utf8('MOLCAS_CURRENT_PROGRAM', 'check')
      self._check_count += 1
      self.write_environment()
      check_input = 'check {0}'.format(self._check_count)
      try:
        self._current_module = Molcas_module(self, 'check', check_input)
        rc = self._current_module.run()
        if (rc == '_RC_NOT_AVAILABLE_'):
          raise Exception
      except:
        (rc, last) = check_test(join(self.scratch, 'molcas_info'), join(self.scratch, 'checkfile'), self._check_count)
        self._last_module = last
      if (self.rc_to_name(rc) != '_RC_ALL_IS_WELL_'):
        self.rc = rc
    self.parallel_task(['x', 'molcas_info'], force=True)

  #TODO: buffer parnell calls
  def parallel_task(self, task, force=False):
    task_type = task[0][0]
    # some replacements for parnell calls in the serial case,
    # if no replacement is available, parnell will be used
    if (self.is_serial):
      # copy files: should also work if the destination is a dir
      if (task_type == 'c'):
        orig = task[2]
        if (not isabs(orig)):
          orig = join(self.scratch, orig)
        dest = task[3]
        if (not isabs(dest)):
          dest = join(self.scratch, dest)
        try:
          copy2(orig, dest)
        # would use SameFileError, but that's only available since python 3.4,
        # so use this workaround
        except Error as e:
          if ('same file' not in text_type(e)):
            raise
        except FileNotFoundError as e:
          if (e.errno == ENOENT):
            if (not force):
              print('Error: file "{0}" not found'.format(orig))
            return -1
          else:
            raise
        return 0
      # remove files: the list is colon-separated, and based on the scratch dir
      if (task_type == 'x'):
        rc = 0
        for i in task[1].split(':'):
          try:
            orig = join(self.scratch, i)
            remove(orig)
          except FileNotFoundError as e:
            if (e.errno == ENOENT):
              if (not force):
                print('Error: file "{0}" not found'.format(orig))
              rc -= 1
            else:
              raise
        return rc
    output = BytesIO()
    error = BytesIO()
    if (task_type == 'b'):
      # No working directory for "base"
      cwd = None
    else:
      cwd = self.scratch
    parnell = self.wrap_command(self.parnell, may_be_serial=False, may_debug=False)
    rc = teed_call(parnell + task, cwd=cwd, stdout=output, stderr=error)
    return rc

  def delete_scratch(self, force=False):
    #TODO: use parnell
    if (self._ready or force):
      if (realpath(self.scratch) == realpath(self.currdir)):
        line = '*** WorkDir and CurrDir are the same, not cleaned! ***'
      else:
        # WorkDir may be a link and not a real directory, so remove its contents
        try:
          filelist = listdir(self.scratch)
        except:
          filelist = []
        for i in [join(self.scratch, x) for x in filelist]:
          if (isdir(i)):
            rmtree(i)
          else:
            remove(i)
        line = '*** WorkDir at {0} cleaned ***'.format(self.scratch)
      print('*'*len(line))
      print(line)
      print('*'*len(line))

  def in_sbin(self, prog):
    '''Return the path of a program in sbin if it exists'''
    if (prog == basename(prog)):
      for path in [self.molcas] + self.sources:
        filename = join(path, 'sbin', prog)
        if (isfile(filename) and access(filename, X_OK)):
          return filename
    return False

  def run_sbin(self, command):
    '''Run a program from the sbin directory if it exists'''
    filename = self.in_sbin(command[0])
    if (filename):
      command[0] = filename
      rc = teed_call(command)
      return rc
    else:
      return None

  def run_logue(self, logue_name):
    '''Run a prologue/epilogue, if found'''
    files = []
    logue_dir = get_utf8('MOLCAS_LOGUE_DIR')
    if (logue_dir):
      if (not isabs(logue_dir)):
        logue_dir = join(self.currdir, logue_dir)
      files.append(join(logue_dir, logue_name))
    files.append(dotmolcas(logue_name))
    for filename in files:
      if (isfile(filename) and access(filename, X_OK)):
        # No rc control here, if it fails, so be it
        output = BytesIO()
        error = BytesIO()
        teed_call(['sh', filename], cwd=self.scratch, stdout=output, stderr=error)

  def check_license(self):
    '''Use molcas.exe to check the license'''
    filename = join(self.molcas, 'bin', 'molcas.exe')
    if (isfile(filename) and access(filename, X_OK)):
      try:
        out = check_output(filename, stderr=STDOUT).decode('utf-8')
        rc = 0
      except CalledProcessError as e:
        out = e.output.decode('utf-8')
        rc = e.returncode-1
      # Capture the licensee
      match = search(r'This copy of MOLCAS is licensed to\s*(.*)\n', out)
      if (match):
        self.licensee = match.group(1)
      return rc
    else:
      return None

#===============================================================================

_prgm_match = re_compile(r'\s*\(prgm\)\s+(\S+)\s+(executable|script)')
_file_match = re_compile(r'\s*\(file\)\s+(\S+)\s+(\S+)\s+(\S+)')
_rc_match = re_compile(r'\$(\w+)\s*=\s*(\d+)\s*;')
_emil_newln = re_compile(r'[ \t]*[=][ \t]*')
_emil_leadb = re_compile(r'^[ \t]*', flags=MULTILINE)
_emil_endofinput = re_compile(r'^end\s*of\s*input', flags=MULTILINE|IGNORECASE)

# TODO: custom .prgm, gracefully fail
def parse_prgm(prgm_file):
  '''Parse a .prgm file to find the module executable and files.'''
  prgm = None
  file_list = {}
  with utf8_open(prgm_file, 'r') as p_f:
    for line in p_f:
      match = _prgm_match.match(line)
      if (match):
        prgm = [expandvars(match.group(1), default='UNKNOWN_VARIABLE').replace('"', ''), match.group(2)[0]]
      match = _file_match.match(line)
      if (match):
        file_list[match.group(1)] = [expandvars(match.group(2), default='UNKNOWN_VARIABLE').replace('"', ''), match.group(3)]
  return (prgm, file_list)

def sanitize(text):
  '''Sanitize input, breaking lines at '=', removing leading blanks, etc'''
  out_text = _emil_leadb.sub('', text)
  out_text = _emil_newln.sub('\n', out_text)
  out_text = _emil_endofinput.sub('End of input', out_text)
  return out_text

class Molcas_module(object):
  '''Setup and run a Molcas module'''

  # Module initialization
  # 1st argument: Molcas environment object
  # 2nd argument: module name
  # 3rd argument (optional): input
  def __init__(self, parent, *args):
    if (len(args) < 1):
      raise MolcasException('Bad initialization of module')
    self.parent = parent
    if (not isinstance(self.parent, Molcas_wrapper)):
      raise MolcasException('The module is not contained in a Molcas_wrapper instance')
    name = args[0].lower()
    if name in self.parent.alias:
      self.name = self.parent.alias[name]
    else:
      self.name = name
    # Read the .prgm file to find the executable and the list of module files
    try:
      prgm_file = join(self.parent.molcas, 'data', self.name + '.prgm')
      (self._exec, self._files) = parse_prgm(prgm_file)
    except FileNotFoundError as e:
      if (e.errno == ENOENT):
        raise MolcasException('Unknown module: {0}'.format(self.name))
      else:
        raise
    prgm_file = join(self.parent.molcas, 'data', 'global.prgm')
    self._files.update(parse_prgm(prgm_file)[1])
    self._links = []
    self.rc = 0
    # Pass the input
    if (len(args) >= 2):
      inp = args[1]
      self.module_input = inp

  # Define the module_input property
  @property
  def module_input(self):
    return self._input
  @module_input.setter
  def module_input(self, inp):
    self._input = sanitize(inp)

  # Get the start line with time
  @property
  def start(self):
    return '--- Start Module: {0} at {1} ---'.format(self.name, self._start.ctime())

  # Get the end line with time and return code
  @property
  def stop(self):
    return '--- Stop Module: {0} at {1} /rc={2} ---'.format(self.name, self._stop.ctime(), self.rc_name)

  # Get the line with time spent (empty if less than 1 second)
  @property
  def span(self):
    span = self._stop - self._start
    d = span.days
    h = span.seconds//3600
    m = (span.seconds-h*3600)//60
    s = span.seconds-h*3600-m*60
    d_name = 'day' if (d==1) else 'days'
    h_name = 'hour' if (h==1) else 'hours'
    m_name = 'minute' if (m==1) else 'minutes'
    s_name = 'second' if (s==1) else 'seconds'
    if (d > 0):
      string = '--- Module {{0}} spent {{1}} {0} {{2}} {1} {{3}} {2} {{4}} {3} ---'.format(d_name, h_name, m_name, s_name)
    elif (h > 0):
      string = '--- Module {{0}} spent {{2}} {1} {{3}} {2} {{4}} {3} ---'.format(d_name, h_name, m_name, s_name)
    elif (m > 0):
      string = '--- Module {{0}} spent {{3}} {2} {{4}} {3} ---'.format(d_name, h_name, m_name, s_name)
    elif (s > 0):
      string = '--- Module {{0}} spent {{4}} {3} ---'.format(d_name, h_name, m_name, s_name)
    else:
      string = ''
    return string.format(self.name, d, h, m, s)

  # Get a tuple with wall/user/system timing
  @property
  def timing(self):
    wall_time = (self._stop - self._start).total_seconds()
    user_time = self._rstop.ru_utime - self._rstart.ru_utime
    system_time = self._rstop.ru_stime - self._rstart.ru_stime
    return (wall_time, user_time, system_time)

  # Get the output stream
  @property
  def output(self):
    return (bytes(self.start, 'utf-8') + b'\n' + self._output.getvalue() + b'\n' + bytes(self.stop, 'utf-8'))

  # Get the error stream
  @property
  def error(self):
    return (bytes(self.start, 'utf-8') + b'\n' + self._error.getvalue() + b'\n' + bytes(self.stop, 'utf-8'))

  # Get the named return code (if available)
  @property
  def rc_name(self):
    return self.parent.rc_to_name(self.rc)

  # Private method to write the input in the $WorkDir
  def _write_input(self, inp):
    stdin = join(self.parent.scratch, 'stdin')
    with utf8_open(stdin, 'w') as i_f:
      print('&' + self.name.upper(), file=i_f)
      if (inp):
        print(inp.strip(), file=i_f)
      print('End of input', file=i_f)
    input_file = self.name.upper()[:5] + 'INP'
    if (input_file in self._files):
      self.parent.parallel_task(['c', '1', stdin, self._files[input_file][0]], force=True)
    self.parent.parallel_task(['c', '1', stdin, self.parent.scratch], force=True)

  def _validate_input(self, inp):
    val = get_utf8('MOLCAS_VALIDATE', default='NO').upper()
    if (val not in ['YES', 'CHECK']):
      return
    if (self.name == 'check'):
      return
    rc, result = validate(inp, self.parent.keywords)
    print('\n-- Input validation for {}'.format(self.name))
    for line in result:
      print(line)
    if (rc == 0):
      print('\nLooks good!')
    elif (rc > 0):
      print('\n*** FAILED! ***')
    print('--')
    if (val == 'YES'):
      self.rc = rc

  # Private method to read the return code from the $WorkDir
  def _read_rc(self):
    rc_local = join(self.parent.scratch, 'rc.local')
    rc_global = join(self.parent.scratch, 'rc.global')
    self.parent.parallel_task(['c', '4', rc_local, rc_global], force=True)
    with utf8_open(rc_global, 'r') as rc_f:
      self.rc = int(rc_f.read())

  def _read_extra_prgm(self):
    prgm_file = join(self.parent.scratch, 'extra.prgm')
    if (exists(prgm_file)):
      self._files.update(parse_prgm(prgm_file)[1])

  def _copy_files(self):
    self._read_extra_prgm();
    dest = self.parent.output
    if (dest == self.parent.scratch):
      return
    files_to_copy = sorted([(k,v[0]) for (k,v) in self._files.items() if 's' in v[1]])
    files_to_move = sorted([(k,v[0]) for (k,v) in self._files.items() if 'm' in v[1]])
    files = self._copy_or_move(copy2, dest, files_to_copy)
    files.extend(self._copy_or_move(move, dest, files_to_move))
    if (len(files) > 0):
      listfiles = ' '.join(files)
      print(fill(listfiles, width=180, break_long_words=False, break_on_hyphens=False,
               initial_indent='*** files: ',
            subsequent_indent='           '))
      print('    saved to directory {0}'.format(dest))

  def _copy_or_move(self, action, dest, filelist):
    #TODO: use parnell
    #TODO: support the '.' attribute
    files = []
    for name, path in filelist:
      if '*' in self._files[name][1]:
        path += '*'
      for i in glob(path):
        if (datetime.fromtimestamp(getmtime(i)) > self._start):
          bi = basename(i)
          j = join(dest, bi)
          if (exists(j)):
            if (isfile(j)):
              if (self.parent.save_mode == 'repl'):
                action(i, j)
              elif (self.parent.save_mode == 'orig'):
                orig = j+'.orig'
                if (not exists(orig)):
                  move(j, orig)
                action(i, j)
              elif (self.parent.save_mode == 'incr'):
                fmt = '{0}.#{1}#'
                n = 1
                jj = fmt.format(j, n)
                while (exists(jj)):
                  n += 1
                  jj = fmt.format(j, n)
                action(i, jj)
              files.append(bi)
            else:
              pass
          else:
            action(i, j)
            files.append(bi)
    return files

  def _delete_files(self):
    files_to_delete = [(k,v[0]) for (k,v) in self._files.items() if 'p' in v[1]]
    rmlist = ['extra.prgm']
    for name, path in files_to_delete:
      if '*' in self._files[name][1]:
        path += '*'
      for i in glob(path):
        rmlist.append(relpath(i, self.parent.scratch))
    for i in glob(join(self.parent.scratch, 'purge*')):
      rmlist.append(relpath(i, self.parent.scratch))
    for i in self._links:
      rmlist.append(relpath(i, self.parent.scratch))
    if (len(rmlist) > 0):
      self.parent.parallel_task(['x', ':'.join(rmlist)], force=True)

  def _make_links(self):
    if (self.name == 'check'):
      return
    if (self.name == 'scf'):
      orblist = ['GVORB', 'LOCORB', 'UHFORB', 'SCFORB', 'RASORB', 'PT2ORB', 'UNAORB', 'CIORB01', 'CPFORB', 'SIORB', 'GSSORB']
    elif (self.name == 'grid_it'):
      orblist = ['GVORB', 'LOCORB', 'SIORB', 'PT2ORB', 'CIORB01', 'CPFORB', 'RASORB', 'MP2ORB', 'UHFORB', 'UNAORB', 'SCFORB', 'GSSORB']
    else:
      orblist = ['GVORB', 'LOCORB', 'RASORB', 'PT2ORB', 'UNAORB', 'SCFORB', 'CIORB01', 'CPFORB', 'SIORB', 'GSSORB']
    #TODO: fail gracefully if file does not exist
    prgm_file = join(self.parent.molcas, 'data', 'inporb_files.prgm')
    filelist = parse_prgm(prgm_file)[1]
    orblist = [filelist[x][0] for x in orblist if x in filelist]
    print('')
    inporb = join(self.parent.scratch, 'INPORB')
    if (not exists(inporb)):
      for orbfile in orblist:
        if (exists(orbfile)):
          symlink(orbfile, inporb)
          self.parent.parallel_task(['c', '1', inporb, self.parent.scratch], force=True)
          self._links.append(inporb)
          print('*** symbolic link created: INPORB -> {0}'.format(basename(orbfile)))
          break

  # Run the module
  def run(self):
    if (self._exec is None):
      return '_RC_NOT_AVAILABLE_'

    self._write_input(self.module_input)

    if (self.name == 'gateway'):
      runfiles = ['RUNFILE']
      if ('RUNFILE' in self._files):
        runfiles.append(relpath(self._files['RUNFILE'][0], self.parent.scratch))
      self.parent.parallel_task(['x', ':'.join(runfiles)], force=True)
    elif (self.name == 'loop'):
      set_utf8('MOLCAS_REDUCE_PRT', 'NO')

    command = self.parent.wrap_command(self._exec[0])
    no_tee=(get_utf8('MOLCAS_DEBUGGER') is not None)

    self._output = BytesIO()
    self._error = BytesIO()
    self._start = datetime.now()
    self._make_links()
    print(self.start)
    self._rstart = getrusage(RUSAGE_CHILDREN)
    self._validate_input(join(self.parent.scratch, 'stdin'))
    if (self.rc > 0):
      self.rc = '_RC_INPUT_ERROR_'
    else:
      self.parent.run_logue('module.prologue')
      if (isfile(self._exec[0]) and access(self._exec[0], X_OK)):
        teed_call(command, cwd=self.parent.scratch, stdout=self._output, stderr=self._error, no_tee=no_tee)
        self._read_rc()
      else:
        self.rc = '_RC_NOT_AVAILABLE_'
      self.parent.run_logue('module.epilogue')
    self._rstop = getrusage(RUSAGE_CHILDREN)
    self._stop = datetime.now()
    print(self.stop)
    self._copy_files()
    self._delete_files()
    span = self.span
    if (span):
      print(span)

    self.parent.add_resources(self.timing)
    k = self.output

    return self.rc
