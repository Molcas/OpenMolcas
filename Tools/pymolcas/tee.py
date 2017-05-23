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
#***********************************************************************

#***********************************************************************
# This file is based on code posted by J.F. Sebastian on stackoverflow *
# (https://stackoverflow.com/a/4985080)                                *
# and licensed under CC-BY-SA 3.0                                      *
# (https://creativecommons.org/licenses/by-sa/3.0)                     *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)

import sys
from subprocess import Popen, PIPE
from threading import Thread

def tee(infile, *files):
  '''Print "infile" to "files" in a separate thread.'''
  def fanout(infile, *files):
    # Write each input line in all the files,
    # Try as a byte string first, then as utf8
    for line in iter(infile.readline, b''):
      for f in files:
        try:
          f.write(line)
        except TypeError:
          f.write(line.decode('utf8', 'replace'))
    infile.close()
  t = Thread(target=fanout, args=(infile,)+files)
  t.daemon = True
  t.start()
  return t

def teed_call(cmd_args, **kwargs):    
  stdout, stderr = [kwargs.pop(s, None) for s in ['stdout', 'stderr']]
  p = Popen(cmd_args,
            stdout=PIPE if stdout is not None else None,
            stderr=PIPE if stderr is not None else None,
            **kwargs)
  threads = []
  if stdout is not None: threads.append(tee(p.stdout, stdout, sys.stdout))
  if stderr is not None: threads.append(tee(p.stderr, stderr, sys.stderr))
  for t in threads:
    t.join()
  return p.wait()
