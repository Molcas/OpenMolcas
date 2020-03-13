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
try:
  from builtins import super
except ImportError:
  from future.builtins import super
try:
  from six import text_type, python_2_unicode_compatible
except ImportError:
  text_type = str
  def python_2_unicode_compatible(orig):
    return orig

from os.path import isfile
from re import match
from io import BytesIO

from molcas_aux import *
from tee import teed_call
from simpleeval import simple_eval, SimpleEval

#===============================================================================
# Patch SimpleEval to use the "decimal" module for better precision handling

import ast
import decimal

def _pymolcas_eval(self, expr):
  self.expr = expr
  return text_type(self._eval(ast.parse(expr.strip()).body[0].value))

@staticmethod
def _pymolcas_eval_num(node):
  return decimal.Decimal(text_type(node.n))

SimpleEval.eval = _pymolcas_eval
SimpleEval._eval_num = _pymolcas_eval_num

#===============================================================================

def env_print(env, string):
  if (env.echo):
    print(string)

def check_trap(rc_name=None):
  # this codes should not be affected by MOLCAS_TRAP
  if (rc_name in ['_RC_ALL_IS_WELL_', '_RC_CHECK_ERROR_']):
    return True
  elif (get_utf8('MOLCAS_TRAP', default='on').lower() == 'off'):
    print('*********************************************************')
    print('**                                                     **')
    print('** Non-zero return code, and MOLCAS_TRAP is set to OFF **')
    print('**                                                     **')
    print('*********************************************************')
    return False
  else:
    return True

@python_2_unicode_compatible
class Expression(object):
  def __init__(self, var, op, val):
    self.var = var
    self.op = op
    self.val = val
  def __bool__(self):
    if (self.op == 'eq'):
      return (expandvars(self.var) == expandvars(self.val))
    if (self.op == 'ne'):
      return (expandvars(self.var) != expandvars(self.val))
    if (self.op == 'file'):
      return self.env.file_exists(self.val)
  def __str__(self):
    var = self.var
    if (var == '$MOLCAS_ITER'):
      var = 'ITER'
    if (self.op == 'eq'):
      return '{0} == {1}'.format(var, self.val)
    if (self.op == 'ne'):
      return '{0} != {1}'.format(var, self.val)
    if (self.op == 'file'):
      return '-FILE {0}'.format(self.val)
  #python2
  __nonzero__ = __bool__

class Statement(object):
  def __init__(self, string=None):
    self.rc = None
    self.level = 0
    self.indent = '  '
    if (string is not None):
      self.string = string
  def run(self, env):
    if (env._goto):
      return None
    print('Unsupported statement: {0}'.format(self.string))
    return self.rc

@python_2_unicode_compatible
class Group(Statement):
  def __init__(self, contents, condition=True, rerun=False, var=None, val=None):
    super().__init__()
    self.contents = contents
    self.condition = condition
    self.rerun = rerun
    self.var = var
    self.val = val
    self.maxiter = 1
    self.grouptype = self._find_type()
  def _find_type(self):
    if (isinstance(self.condition, Expression) and (self.rerun is False)):
      return 'if'
    elif ((self.condition is True) and (self.rerun is True)):
      if (self.var is not None):
        return 'foreach'
      else:
        return 'do'
    elif ((self.condition is True) and (self.rerun is False)):
      return 'group'
    else:
      return 'unknown'
  def __str__(self):
    if (self.grouptype == 'group'):
      init = ''
      end = ''
    elif (self.grouptype == 'foreach'):
      init = self.level*self.indent + '>>> FOREACH {0} IN ({1})\n\n'.format(self.var, self.val)
      end = '\n\n' + self.level*self.indent + '>>> END FOREACH'
    elif (self.grouptype == 'do'):
      init = self.level*self.indent + '>>> DO WHILE\n\n'
      end = '\n\n' + self.level*self.indent + '>>> END DO'
    elif (self.grouptype == 'if'):
      if ((len(self.contents) == 1) and isinstance(self.contents[0], Jump)):
        init = self.level*self.indent + '>>> IF ({0}) GOTO {1}'.format(self.condition, self.contents[0].name)
        return init
      init = self.level*self.indent + '>>> IF ({0})\n\n'.format(self.condition)
      end = '\n\n' + self.level*self.indent + '>>> END IF'
    else:
      init = '('
      end = ')'
    return '{0}{1}{2}'.format(init, '\n\n'.join(map(text_type, self.contents)), end)
  def run(self, env):
    if (env._goto):
      return None
    if (self.grouptype == 'foreach'):
      self.values = expandvars(self.val, default='UNKNOWN_VARIABLE')
      re_match = match(r'(-?\d+)\s*\.\.\s*(-?\d+)', self.values)
      if (re_match):
        ini = int(re_match.group(1))
        fin = int(re_match.group(2))
        sign = -1 if (fin < ini) else 1
        self.values = range(ini, fin + sign, sign)
        self.values = list(map(text_type, self.values))
      else:
        self.values = [x.strip() for x in self.values.split(',')]
      self.maxiter = len(self.values)
      env.enter_loop('foreach')
    elif (self.grouptype == 'do'):
      self.previter = get_utf8('MOLCAS_ITER', default='0')
      self.maxiter = int(get_utf8('MOLCAS_MAXITER', default='50'))
      env.enter_loop('do')
    self.thisiter = 1
    if (isinstance(self.condition, Expression)):
      self.condition.env = env
    if (isinstance(self.rerun, Expression)):
      self.rerun.env = env
    while (self.thisiter <= self.maxiter):
      if (self.grouptype == 'if'):
        env_print(env, '\n>>> IF ({0})'.format(self.condition))
      itercontents = self.contents[:]
      if (self.condition):
        i = 0
        if (self.grouptype == 'foreach'):
          Assignment(self.var, self.values[self.thisiter-1]).run(env)
        elif (self.grouptype == 'do'):
          set_utf8('MOLCAS_ITER', self.thisiter)
        no_break = True
        if (self.grouptype == 'do'):
          rc_name = '_RC_CONTINUE_LOOP_'
        else:
          rc_name = '_RC_ALL_IS_WELL_'
        while ((i < len(itercontents)) and (no_break)):
          item = itercontents[i]
          item.group = itercontents
          rc = item.run(env)
          # Transparent return: only update return code if not None
          if (rc is not None):
            self.rc = rc
            rc_name = env.rc_to_name(self.rc)
          no_break = rc_name in ['_RC_ALL_IS_WELL_', '_RC_CONTINUE_LOOP_', '_RC_CONTINUE_UNIX_LOOP_']
          # check MOLCAS_TRAP here only if it's not the last element in the group
          if (i+1 < len(itercontents)):
            if ((not no_break) and (not check_trap(rc_name))):
              no_break = True
              rc_name = '_RC_ALL_IS_WELL_'
          set_utf8('EMIL_RETURNCODE', 0 if self.rc is None else self.rc)
          i += 1
      else:
        env_print(env, '(Skipped)')
      # Termination conditions
      # TODO: should MOLCAS_TRAP be checked here too?
      # TODO: run EMIL commands after a module fails (see FS#37)
      if (env._goto):
        break
      # this seems to be needed (here) for compatibility with old driver
      if (self.grouptype in ['do', 'foreach']):
        env._check_count += 1
      if (not self.rerun):
        break
      # in a DO loop, "all is well" means the loop terminates
      if (self.grouptype == 'do'):
        no_break = rc_name in ['_RC_CONTINUE_LOOP_', '_RC_CONTINUE_UNIX_LOOP_']
      elif (self.grouptype == 'foreach'):
        no_break = rc_name in ['_RC_ALL_IS_WELL_', '_RC_CONTINUE_LOOP_', '_RC_CONTINUE_UNIX_LOOP_']
      else:
        no_break = rc_name in ['_RC_ALL_IS_WELL_']
      # no-trap and loops are tricky
      if ((not no_break) and (not check_trap(rc_name))):
        if ((rc_name not in ['_RC_NOT_CONVERGED_']) or (self.grouptype != 'do')):
          no_break = True
          rc_name = '_RC_ALL_IS_WELL_'
      if (self.grouptype == 'do'):
        env_print(env, '\n>>> END DO')
      elif (self.grouptype == 'foreach'):
        env_print(env, '\n>>> END FOREACH')
      if (not no_break):
        break
      self.thisiter += 1
    if (self.grouptype == 'do'):
      if (not env._goto):
        if (rc_name in ['_RC_CONTINUE_LOOP_', '_RC_CONTINUE_UNIX_LOOP_']):
          self.rc = '_RC_NOT_CONVERGED_'
      set_utf8('MOLCAS_ITER', self.previter)
      env.exit_loop()
    if (self.grouptype == 'foreach'):
      env.exit_loop()
    set_utf8('EMIL_RETURNCODE', self.rc)
    return self.rc

@python_2_unicode_compatible
class Assignment(Statement):
  def __init__(self, var, val, literal=True):
    super().__init__()
    self.literal = literal
    self.var = var
    self.val = val
  def __str__(self):
    if (self.literal):
      return self.level*self.indent + '>>> EXPORT {0} = {1}'.format(self.var, self.val)
    else:
      return self.level*self.indent + '>>> EVAL {0} = {1}'.format(self.var, self.val)
  def run(self, env):
    if (env._goto):
      return None
    if (self.literal):
      val = expandvars(self.val, default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> EXPORT {0} = {1}'.format(self.var, val))
      set_utf8(self.var, val)
    else:
      val = expandvars(self.val, default='0')
      try:
        eval_val = simple_eval(val)
      except:
        eval_val = ''
      env_print(env, '\n>>> EVAL {0} = {1} = {2}'.format(self.var, val, eval_val))
      set_utf8(self.var, eval_val)
    return self.rc

@python_2_unicode_compatible
class Break(Statement):
  def __init__(self, rc):
    super().__init__()
    self.rc = rc
  def __str__(self):
    num = ' {0}'.format(self.rc)
    if (self.rc is None):
      num = ''
    return self.level*self.indent + '>>> EXIT{0}'.format(num)
  def run(self, env):
    if (env._goto):
      return None
    try:
      num = expandvars(self.rc)
      rc = int(num)
    except:
      rc = None
    if (rc is None):
      rc = '_RC_EXIT_'
    if (rc == '0'):
      rc = '_RC_EXIT_EXPECTED_'
    env_print(env, '\n>>> EXIT {0}'.format(env.rc_to_name(rc)))
    return rc

@python_2_unicode_compatible
class Include(Statement):
  def __init__(self, filename):
    super().__init__()
    self.filename = filename
  def __str__(self):
    return self.level*self.indent + '>>> INCLUDE {0}'.format(self.filename)
  def run(self, env):
    if (env._goto):
      return None
    # NOTE: circular dependency here
    import emil_parse
    pos = self.group.index(self)+1
    fname = expandvars(self.filename, default='UNKNOWN_VARIABLE')
    env_print(env, '\n>>> INCLUDE {0}'.format(fname))
    if (isfile(fname)):
      blocks, files = emil_parse.EMIL_Parse(fname)
      env.add_files(files)
      self.group.insert(pos, blocks)
    else:
      self.rc = '_RC_INPUT_EMIL_ERROR_'
    return self.rc

@python_2_unicode_compatible
class Program(Statement):
  def __init__(self, lines):
    super().__init__()
    self.lines = lines
  def __str__(self):
    name = '&' + self.lines[1].upper()
    code = ('\n' + self.level*self.indent).join(self.lines[2:])
    if (code != ''):
      code = '\n' + self.level*self.indent + code
    return self.level*self.indent + name + code
  def run(self, env):
    if (env._goto):
      return None
    name = self.lines[1]
    inp = expandvars('\n'.join(self.lines[2:]), default='UNKNOWN_VARIABLE')
    rc = env.run_module(name, inp)
    if (rc is not None):
      self.rc = rc
    return self.rc

# TODO: convert to ParTask?
@python_2_unicode_compatible
class System(Statement):
  def __init__(self, commandline, parallel=False):
    super().__init__()
    self.commandline = commandline
    self.parallel = parallel
  @property
  def expanded_commandline(self):
    return expandvars(self.commandline, skip_escaped=True)
  @property
  def name(self):
    if (self.parallel):
      return "EXEC"
    else:
      return "SHELL"
  def __str__(self):
    return self.level*self.indent + '>>> {0} {1}'.format(self.name, self.commandline)
  def run(self, env):
    if (env._goto):
      return None
    output = BytesIO()
    error = BytesIO()
    env_print(env, '\n>>> {0} {1}'.format(self.name, self.expanded_commandline))
    if (env.allow_shell):
      if (self.parallel):
        task = ['!'] + self.expanded_commandline.split()
        self.rc = env.parallel_task(task)
      else:
        self.rc = teed_call(self.expanded_commandline, shell=True, cwd=env.scratch, stdout=output, stderr=error)
    else:
      env_print(env, '(Shell commands disabled)')
    # Successful return is transparent: it inherits the previous return code
    if (self.rc == 0):
      self.rc = None
    return self.rc

@python_2_unicode_compatible
class Setting(Statement):
  def __init__(self, var, val):
    super().__init__()
    self.var = var
    self.val = val
  def __str__(self):
    if (self.var == 'echo'):
      val = self.val
      if (val in ['on', 'off']):
        val = val.upper()
      return self.level*self.indent + '>>> ECHO {0}'.format(val)
  def run(self, env):
    if (env._goto):
      return None
    if (self.var == 'echo'):
      if (self.val == 'on'):
        env.echo = True
      elif (self.val == 'off'):
        env.echo = False
      else:
        print(self.val)
    return self.rc

@python_2_unicode_compatible
class Label(Statement):
  def __init__(self, name):
    super().__init__()
    self.name = name
  def __str__(self):
    return self.level*self.indent + '>>> LABEL {0}'.format(self.name)
  def run(self, env):
    if (env._goto):
      if (env._goto == self.name):
        env._goto = False
    return self.rc

@python_2_unicode_compatible
class Jump(Statement):
  def __init__(self, name):
    super().__init__()
    self.name = name
  def __str__(self):
    return self.level*self.indent + '>>> GOTO {0}'.format(self.name)
  def run(self, env):
    if (env._goto):
      return None
    env_print(env, '\n>>> GOTO {0}'.format(self.name))
    env._goto = self.name
    return self.rc

@python_2_unicode_compatible
class ParTask(Statement):
  def __init__(self, action, force, args):
    super().__init__()
    self.action = action
    self.force = force
    self.args = args
  def __str__(self):
    f = '-FORCE ' if self.force else ''
    if (self.action == 'rm'):
      return self.level*self.indent + '>>> RM {0}{1}'.format(f, self.args[0])
    if (self.action == 'copy'):
      return self.level*self.indent + '>>> COPY {0}{1} {2}'.format(f, self.args[0], self.args[1])
    if (self.action == 'save'):
      return self.level*self.indent + '>>> SAVE {0}{1} {2}'.format(f, self.args[0], self.args[1])
    if (self.action == 'clone'):
      return self.level*self.indent + '>>> CLONE {0}{1} {2}'.format(f, self.args[0], self.args[1])
    if (self.action == 'collect'):
      return self.level*self.indent + '>>> COLLECT {0}{1} {2}'.format(f, self.args[0], self.args[1])
  def run(self, env):
    if (env._goto):
      return None
    f = '-FORCE ' if self.force else ''
    if (self.action == 'rm'):
      fname = expandvars(self.args[0], default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> RM {0}{1}'.format(f, fname))
      task = ['x', fname]
    if (self.action == 'copy'):
      f_in = expandvars(self.args[0], default='UNKNOWN_VARIABLE')
      f_out = expandvars(self.args[1], default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> COPY {0}{1} {2}'.format(f, f_in, f_out))
      task = ['c', '1', f_in, f_out]
    if (self.action == 'save'):
      f_in = expandvars(self.args[0], default='UNKNOWN_VARIABLE')
      f_out = expandvars(self.args[1], default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> SAVE {0}{1} {2}'.format(f, f_in, f_out))
      task = ['c', '0', f_in, f_out]
    if (self.action == 'clone'):
      f_in = expandvars(self.args[0], default='UNKNOWN_VARIABLE')
      f_out = expandvars(self.args[1], default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> CLONE {0}{1} {2}'.format(f, f_in, f_out))
      task = ['c', '3', f_in, f_out]
    if (self.action == 'collect'):
      f_in = expandvars(self.args[0], default='UNKNOWN_VARIABLE')
      f_out = expandvars(self.args[1], default='UNKNOWN_VARIABLE')
      env_print(env, '\n>>> COLLECT {0}{1} {2}'.format(f, f_in, f_out))
      task = ['c', '2', f_in, f_out]
    self.rc = env.parallel_task(task, force=self.force)
    if (self.force):
      self.rc = 0
    # Successful return is transparent: it inherits the previous return code
    if (self.rc == 0):
      self.rc = None
    else:
      self.rc = '_RC_INPUT_EMIL_ERROR_'
    return self.rc

@python_2_unicode_compatible
class Python(Statement):
  '''
  This is not really part of the "abstract flow" family, but it's included
  for consistency. Intended to be used with arbitrary python code, instead
  of some specific input language (like EMIL).
  '''
  def __init__(self, filename):
    self.rc = None
    with utf8_open(filename) as f_in:
      self.text = f_in.read()
      # If there's an error in the code SyntaxError will be raised
      # and catched upstream
      self.code = compile(self.text, filename, 'exec')
  def __str__(self):
    return self.text.rstrip()
  def run(self, env):
    exec(self.code, {'MOLCAS': env})
    self.rc = 0
    return self.rc
