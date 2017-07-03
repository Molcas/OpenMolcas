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
# Copyright (C) 2015,2017, Ignacio Fdez. Galván                        *
#***********************************************************************

'''
This exports pymolcas to a single-file python script that can be copied around.
It still uses python and it's trivial to recover the original files (although
without comments). Obfuscation is not the goal here, just getting something
that's easy to run, move and distribute.
'''

import sys, zlib, base64, os, stat
sys.dont_write_bytecode = True

files = ['tee', 'molcas_aux', 'emil_grammar', 'abstract_flow', 'emil_parse', 'python_parse', 'check_test', 'molcas_wrapper', 'pymolcas']
try:
  exe_name = sys.argv[1]
except:
  exe_name = 'pymolcas'
obfuscate = False # only obfuscate the "wrapper", not the modules' code (yet)
try:
  import pyminifier.minification as pm
  compact = True
except ImportError:
  compact = False
compress_and_b64 = True

def wrap(text, width):
  lines = []
  l = 0
  while (l < len(text)):
    r = min(l+width,len(text)+1)
    lines.append(text[l:r])
    l = r
  return '\n'.join(lines)

def minify(string):
  string = pm.remove_comments_and_docstrings(string)
  string = pm.multiline_indicator.sub('', string)
  string = pm.fix_empty_methods(string)
  string = pm.join_multiline_pairs(string)
  string = pm.join_multiline_pairs(string, '[]')
  string = pm.join_multiline_pairs(string, '{}')
  string = pm.remove_blank_lines(string)
  string = pm.reduce_operators(string)
  string = pm.dedent(string)
  return string

if (compress_and_b64):
  code = 'import zlib,base64;exec(zlib.decompress(base64.b64decode(m[1])),module.__dict__);del zlib,base64'
else:
  code = 'exec(m[1],module.__dict__)'

hexcode = ''.join('{:02x}'.format(c) for c in bytes(code, 'ascii'))

from pymolcas import warning

mods_name = 'M' if obfuscate else 'modules'

failed = False

with open(exe_name, 'w', encoding='utf-8') as f:
  f.write('''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#{0}
{1} = [

'''.format(warning.replace('\n','\n#'), mods_name))

  try:

    for i in files:
      filename = '{0}.py'.format(i)
      with open(filename, 'r', encoding='utf-8') as fin:
        content = fin.read()
      if (compact):
        content = minify(content)
      if (compress_and_b64):
        content = str(base64.b64encode(zlib.compress(bytes(content, 'utf8'))), 'utf8')
        content = wrap(content, 120)
        fmt = "  ['{0}', '''\n{1}\n'''],\n\n"
      else:
        content = bytes(content, 'utf8')
        fmt = "  ['{0}', \n{1}\n],\n\n"
      if (obfuscate):
        name_i = ''.join('{:02x}'.format(c) for c in bytes(i, 'ascii'))
      else:
        name_i = i
      f.write(fmt.format(name_i, content))

    f.write(']\n')

    if (obfuscate):
      f.write('''
checksum = '\''
{0}
'\''.replace('\\n', '')
'''.format(wrap(hexcode,120)))
      mod_name_code = 'binascii.unhexlify(m[0]).decode(\'ascii\')'
      mod_code = 'binascii.unhexlify(checksum)'
    else:
      mod_name_code = 'm[0]'
      mod_code = '\'{0}\''.format(code)

    f.write('''
import sys, shutil, os.path, types, binascii
for m in {0}:
  x = {1}
  module = types.ModuleType(x)
  exec({2})
  sys.modules[x] = module
del {0}, types, binascii

from pymolcas import main

sys.exit(main(os.path.realpath(shutil.which(sys.argv[0]))))
'''.format(mods_name, mod_name_code, mod_code))

  except Exception as e:
    print(str(e))
    failed = True

if (failed):
  os.remove(exe_name)
  sys.exit(1)
else:
  st = os.stat(exe_name)
  os.chmod(exe_name, st.st_mode | stat.S_IEXEC)

sys.exit(0)
