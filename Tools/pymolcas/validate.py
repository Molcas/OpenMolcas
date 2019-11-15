#!/usr/bin/env python
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
# Copyright (C) 2019, Ignacio Fdez. Galván                             *
#***********************************************************************

from __future__ import (unicode_literals, division, absolute_import, print_function)
import re
import shlex
try:
  from lxml import etree as ET
except ImportError:
  pass

# Global variables
class gv:
  current_name = None
  lookup = { 'NSYM': None }
  syms = [ 1, 2, 4, 8 ]

# Basic expression to find keywords:
#  KEYWORD elements
#  GROUP elements
#  SELECT/KEYWORD elements
#  all of the above inside GROUP[KIND=BOX] (just visual grouping)
kw_exp = '(.|GROUP[@KIND="BOX"])/GROUP | (.|GROUP[@KIND="BOX"])/KEYWORD | (.|GROUP[@KIND="BOX"])/SELECT/KEYWORD'

# Compare two strings, case-insensitive, trimmed/padded to specified length
# (by default, 4 characters, typical of Molcas keywords)
def cmp_str(string1, string2, n=4):
  s1 = re.sub(r'\s+', ' ', string1.strip().upper()).ljust(n)[0:n]
  s2 = re.sub(r'\s+', ' ', string2.strip().upper()).ljust(n)[0:n]
  return s1 == s2

# Convert an integer number from Fortran format
def fortran_int(num):
  # pretend environment variables are always good
  if (num[0] == '$'):
    return 1
  return int(num)

# Convert a real number from Fortran format
def fortran_float(num):
  # pretend environment variables are always good
  if (num[0] == '$'):
    return 1.0
  # in case there is no exponent marker
  num = re.sub(r'(\d)([+-]\d)', r'\1e\2', num)
  # convert D to E
  num = num.translate(str.maketrans('dD', 'eE'))
  return float(num)

# Convert a list to integers
def to_int(parts):
  return [fortran_int(i) for i in parts]

# Convert a list to floats
def to_float(parts):
  return [fortran_float(i) for i in parts]

# Split in Fortran style
def fortran_split(string):
  return string.replace(',', ' ').split()

# Split in Fortran style, with quotes
def fortran_split_quotes(string):
  return shlex.split(string.replace(',', ' '))

# Get the first word
def first_word(string):
  return string.split()[0].strip()

# Get the first word as an integer
def first_int(string):
  return fortran_int(first_word(string))

# Is line blank?
def blank(string):
  return string.strip() == ''

# Check if the current line starts a group
#  1: group matches, consume one line
#  0: group doesn't match
# -1: group matches, but unhandled type
def test_group(lines, group):
  name = group.get('NAME').upper().strip()
  kind = group.get('KIND')
  gv.current_name = name
  if (cmp_str(name, lines[0])):
    if (kind in ['BLOCK', 'RADIO']):
      return 1
    else:
      return -1
  return 0

# Check if the current line(s) matches a keyword
#    n: the next n lines match the keyword
#   -n: the next n lines match the keyword and end a group
# None: the keyword doesn't match
def test_keyword(lines, keyword):
  name = keyword.get('NAME').upper().strip()
  kind = keyword.get('KIND')
  size = keyword.get('SIZE')
  choice = keyword.get('LIST')
  gv.current_name = name
  # if the keyword is in a group, these are the ways to end the group
  endlist = ['END']
  parent = keyword.getparent()
  if (parent.tag == 'GROUP'):
    endlist.append('END' + parent.get('NAME')[0])
  l = 0
  if (cmp_str(name, lines[0])):
    ll = test(lines, keyword, kind, size, choice)
    # try possible <ALTERNATE> definitions
    if (ll is None):
      for alt in (keyword.xpath('ALTERNATE')):
        kind_ = alt.get('KIND')
        size_ = alt.get('SIZE')
        choice_ = alt.get('LIST')
        ll = test(lines, keyword, kind_, size_, choice_)
        if (ll is not None):
          l += ll
          break
      else:
        return None
    else:
      l += ll
    # if tag is in endlist, invert
    if (name in endlist):
      l *= -1
  else:
    return None
  return l

# Wrapper function to test a keyword with standard or custom kind
def test(lines, keyword, kind, size, choice=None):
  if (kind == 'CUSTOM'):
    return test_custom(lines, keyword)
  else:
    return test_standard(lines, kind, size, choice=choice)

# Check if the current lines match a keyword with standard types
#    n: the next n lines match the keyword
# None: the keyword doesn't match
def test_standard(lines, kind, size, computed=False, choice=None):

  # first line contains the keyword
  l = 1
  # no additional input
  if (kind == 'SINGLE'):
    pass

  # string types
  elif (kind == 'STRINGS'):
    if (size is None):
      return None
    n = int(size)
    if (computed):
      n *= first_int(lines[l])
      l += 1
    l += n
    if (l >= len(lines)):
      return None
  elif (kind == 'STRING'):
    return test_standard(lines, 'STRINGS', 1)
  elif (kind == 'STRINGS_COMPUTED'):
    return test_standard(lines, 'STRINGS', size, computed=True)
  elif (kind == 'STRINGS_LOOKUP'):
    if (gv.lookup.get(size) is not None):
      return test_standard(lines, 'STRINGS', gv.lookup[size])
    else:
      return None

  # int types
  elif (kind == 'INTS'):
    if (computed):
      off = 1
      n = 0
    else:
      off = 0
      n = int(size)
    ints = 0
    # read the necessary number of integers
    while (ints < n+off):
      if (l >= len(lines)):
        return None
      for part in fortran_split(lines[l]):
        try:
          i = fortran_int(part)
          ints += 1
          if ((gv.lookup.get(gv.current_name) is None) and (ints == 1)):
            gv.lookup[gv.current_name] = i
        except:
          return None
        if (computed and (n == 0)):
          n = i*int(size)
          off = 0
          ints = 0
        if (ints == n+off):
          break
      l += 1
  elif (kind == 'INT'):
    return test_standard(lines, 'INTS', 1)
  elif (kind == 'INTS_COMPUTED'):
    return test_standard(lines, 'INTS', size, computed=True)
  elif (kind == 'INTS_LOOKUP'):
    if (gv.lookup.get(size) is not None):
      return test_standard(lines, 'INTS', gv.lookup[size])
    elif (size == 'NSYM'):
      for n in reversed(gv.syms):
        ll = test_standard(lines, 'INTS', n)
        if (ll):
          l = ll
          gv.lookup[size] = n
          break
    elif (size == 'ANY'):
      n = 1
      ll = 0
      while (ll is not None):
        ll = test_standard(lines, 'INTS', n)
        if (ll):
          l = ll
          n += 1
        else:
          if (n == 1):
            return None

  # real types
  elif (kind == 'REAL'):
    return test_standard(lines, 'REALS', 1)
  elif (kind == 'REALS'):
    if (computed):
      off = 1
      n = 0
    else:
      off = 0
      n = int(size)
    nums = 0
    # read the necessary number of reals
    while (nums < n+off):
      if (l >= len(lines)):
        return None
      for part in fortran_split(lines[l]):
        try:
          if (off == 0):
            i = fortran_float(part)
          else:
            i = fortran_int(part)
          nums += 1
        except:
          return None
        if (computed and (n == 0)):
          n = i*int(size)
          off = 0
          nums = 0
        if (nums == n+off):
          break
      l += 1
  elif (kind == 'REALS_COMPUTED'):
    return test_standard(lines, 'REALS', size, computed=True)
  elif (kind == 'REALS_LOOKUP'):
    if (gv.lookup.get(size) is not None):
      return test_standard(lines, 'REALS', gv.lookup[size])
    elif (size == 'NSYM'):
      for n in reversed(gv.syms):
        ll = test_standard(lines, 'REALS', n)
        if (ll):
          l = ll
          gv.lookup[size] = n
          break
    elif (size == 'DEG_FREEDOM'):
      n = 1
      ll = 0
      while (ll is not None):
        ll = test_standard(lines, 'REALS', 3*n)
        if (ll):
          l = ll
          gv.lookup[size] = n
          n += 1

  # choice type
  elif (kind == 'CHOICE'):
    opts = [i.upper().strip() for i in choice.split(',')]
    opts = [i.split(':')[0].strip() for i in opts if (i != '----')]
    ll = 0
    for opt in opts:
      # an option can contain several lines, separated by ;
      words = [i.strip() for i in opt.split(';')]
      ll = 0
      for i,w in enumerate(words):
        word = first_word(lines[i+1])
        # pretend an environment variable will always match
        if ((word.strip()[0] != '$') and (word.strip().upper() != w)):
          ll = 0
          break
        ll += 1
      if (ll > 0):
        break
    if (ll == 0):
      return None
    l += ll

  # unrecognized type
  else:
    return None
  return l

# Check if the current lines match a keyword with custom type,
# these are keywords with a syntax that doesn't comply with the standard types
#    n: the next n lines match the keyword
# None: the keyword doesn't match
def test_custom(lines, keyword):

  name = keyword.get('NAME')
  module = keyword.get('MODULE')
  node = keyword
  # first line contains the keyword
  l = 1

  if (module == 'CASPT2'):
    if (name == 'AFREEZE'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        nums = to_float(parts[1:3])
        assert (len(nums) == 2)
        l += 1
        parts = fortran_split(lines[l])
        assert (len(parts) == n)
        l += 1
      except:
        return None
    elif (name == 'EFFE'):
      try:
        n = first_int(lines[l])
        ll = test_standard(lines, 'REALS', n*n)
        if (ll):
          l = ll
        else:
          return None
      except:
        return None
    elif (name in ['MULTISTATE', 'XMULTISTATE']):
      try:
        assert (first_word(lines[l]).upper() == 'ALL')
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'CASVB'):
    if (name == 'ORB'):
      try:
        l -= 1
        n = 0
        parts = fortran_split(lines[l])
        parts.pop(0)
        if (len(parts) < 1):
          l += 1
          parts = fortran_split(lines[l])
        n = fortran_int(parts.pop(0))
        if (len(parts) < 1):
          l += 1
          parts = fortran_split(lines[l])
        nums = to_float(parts)
        l += 1
        while (l < len(lines)):
          parts = fortran_split(lines[l])
          try:
            nums = to_float(parts)
            l += 1
          except:
            break
      except:
        return None
    elif (name == 'ORBREL'):
      try:
        parts = fortran_split(lines[l])
        assert (len(parts) >= 3)
        nums = to_int(parts[0:2])
        l += 1
      except:
        return None
    elif (name == 'FIXORB'):
      try:
        assert (first_word(lines[l]).upper() in ['ALL', 'NONE'])
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'CHT3'):
    if (name == 'TITLE'):
      kwlist = [kw.get('NAME') for kw in keyword.getparent().xpath(kw_exp)]
      ll = 1
      # all lines (max 10) until a keyword matches
      while ((l+ll < len(lines)) and (ll < 10)):
        if (any(cmp_str(lines[l+ll], kw) for kw in kwlist)):
          break
        ll += 1
      l += ll
    else:
      return None

  elif (module == 'DMRGSCF'):
    if (name == 'SOCCUPY'):
      ll = 1
      while (l+ll < len(lines)):
        try:
          parts = fortran_split(lines[l+ll])
          assert all([i in '0ud2' for i in parts])
          ll += 1
        except:
          break
      if (ll == 0):
        return None
      l += ll
    else:
      return None

  elif (module == 'ESPF'):
    if (name == 'EXTERNAL'):
      try:
        n = first_int(lines[l])
        l += 1
        if (n == 0):
          try:
            ll = 0
            while (l+ll < len(lines)):
              parts = fortran_split(lines[l+ll])
              i = fortran_int(parts[0])
              nums = to_float(parts[1:11])
              assert (len(nums) == 10)
              l += 1
          except:
            pass
          l += ll
        elif (n > 0):
          for ll in range(n):
            nums = to_float(fortran_split(lines[l])[0:4])
            l += 1
        else:
          return None
      except:
        return None
    else:
      return None

  elif (module == 'EXTF'):
    if (name == 'LINEAR'):
      try:
        n = first_int(lines[l])
        l += 1
        n = first_int(lines[l])
        l += 1
        n = fortran_float(first_word(lines[l]))
        l += 1
        n = first_int(lines[l])
        assert (n in [0, 1])
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'FFPT'):
    if (name in ['DIPO', 'QUAD', 'OCTU', 'EFLD', 'EFGR']):
      ll = 0
      while (l+ll < len(lines)):
        try:
          parts = fortran_split(lines[l+ll])
          i = fortran_float(parts[1])
          if (name in ['DIPO', 'EFLD']):
            assert (parts[0].upper() in ['X', 'Y', 'Z'])
          elif (name in ['QUAD', 'EFGR']):
            assert (parts[0].upper() in ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'])
          elif (name in ['OCTU']):
            assert (parts[0].upper() in ['XXX', 'XXY', 'XXZ', 'XYY', 'XYZ', 'XZZ', 'YYY', 'YYZ', 'YZZ', 'ZZZ'])
          ll += 1
        except:
          break
      if (ll == 0):
        return None
      else:
        l += ll
      if (name != 'DIPO'):
        try:
          parts = fortran_split(lines[l].replace(',', ' '))
          assert (parts[0].upper() == 'ORIG')
          nums = to_float(parts[1:4])
          assert (len(nums) == 3)
          l += 1
        except:
          if (name in ['EFLD', 'EFGR']):
            return None
    elif (name == 'GLBL'):
      while (l < len(lines)):
        try:
          parts = fortran_split_quotes(lines[l])
          assert (len(parts[0]) == 8)
          i = fortran_int(parts[1])
          i = fortran_float(parts[2])
          l += 1
        except:
          break
    elif (name == 'SELE'):
      try:
        n = first_int(lines[l])
        l += 1
        for i in range(n):
          parts = fortran_split(lines[l])
          assert (parts[0].upper() in ['.TRUE.', '.FALSE.', 'T', 'F'])
          n1, n2 = to_int(parts[1:3])
          l += 1
        for i in range(n-1):
          parts = fortran_split(lines[l])
          assert all([ii.upper() in ['.TRUE.', '.FALSE.', 'T', 'F'] for ii in parts[0:i-1]])
          l += 1
        parts = fortran_split(lines[l])
        nums = to_float(parts[0:3])
        assert (len(nums) == 3)
        l += 1
      except:
        return None
    else:
      return None

  elif (module in ['GATEWAY', 'SEWARD']):
    if (name == 'COORD'):
      try:
        assert (first_word(lines[l])[0] != '$')
        n = first_int(lines[l])
      except:
        n = 0
      l += 1
      if (n > 0):
        try:
          l += 1
          for ll in range(n):
            parts = fortran_split(lines[l])
            nums = to_float(parts[1:4])
            assert (len(nums) == 3)
            l += 1
        except:
          return None
    elif (name == 'BASIS (NATIVE)'):
      try:
        while (blank(lines[l])):
          l += 1
        opts = lines[l].split('/')
        l += 1
        end = False
        if ((len(opts) > 1) and (opts[1].strip().upper() == 'INLINE')):
          while (l < len(lines)):
            if (blank(lines[l])):
              l += 1
              continue
            parts = fortran_split(lines[l])
            if (parts[-1] == '/'):
              parts.pop(-1)
            try:
              nums = to_float(parts)
              l += 1
            except:
              break
        if ('FRAGMENT' in opts[0].strip().upper()):
          ncoor = 0
          nener = 0
          while (l < len(lines)):
            parts = fortran_split(lines[l])
            if (blank(lines[l])):
              pass
            elif (parts[0].upper() == 'LBASIS'):
              ll = test_standard(lines[l:], 'STRINGS_COMPUTED', 1)
              if (ll):
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'RELCOORDS'):
              ll = test_standard(lines[l:], 'REALS_COMPUTED', 3)
              if (ll):
                ncoor = fortran_int(lines[l+1].split()[0])
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'ENERGIES'):
              ll = test_standard(lines[l:], 'REALS_COMPUTED', 1)
              if (ll):
                nener = fortran_int(lines[l+1].split()[0])
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'MOCOEFF'):
              ll = test_standard(lines[l:], 'REALS_COMPUTED', nener)
              if (ll):
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'MULLIKEN'):
              ll = test_standard(lines[l:], 'REALS', ncoor)
              if (ll):
                l += ll-1
              else:
                return None
            else:
              break
            l += 1
        elif ('ECP' in opts[0].strip().upper()):
          while (l < len(lines)):
            parts = fortran_split(lines[l])
            if (blank(lines[l])):
              pass
            elif (parts[0].upper() in ['M1', 'M2']):
              ll = test_standard(lines[l:], 'REALS_COMPUTED', 2)
              if (ll):
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'COREREP'):
              ll = test_standard(lines[l:], 'REALS', 1)
              if (ll):
                l += ll-1
              else:
                return None
            elif (parts[0].upper() == 'PROJOP'):
              ll = 1
              n = first_int(lines[l+ll])
              ll += 1
              for i in range(n+1):
                n1, n2 = to_int(fortran_split(lines[l+ll])[0:2])
                ll += test_standard(lines[l+ll:], 'REALS', n2+n1*(1+n2))
              if (ll):
                l += ll-1
            else:
              break
            l += 1
        n = 0
        while (l < len(lines)):
          if (cmp_str(lines[l], 'END')):
            l += 1
            break
          parts = fortran_split(lines[l])
          if (blank(lines[l])):
            pass
          elif any([cmp_str(parts[0], i) for i in ['CHARGE', 'NUCLEAR']]):
            ll = test_standard(lines[l:], 'REALS', 1)
            if (ll):
              l += ll-1
            else:
              return None
          elif any([cmp_str(parts[0], i) for i in ['SPHERICAL', 'CARTESIAN', 'CONTAMINANT']]):
            assert (len(parts) > 1)
          elif any([cmp_str(parts[0], i) for i in ['MUON', 'PSEUDO', 'FRAGMENT']]):
            pass
          elif (cmp_str(parts[0], 'SPEC')):
            l += 1
            while (l < len(lines)):
              if (cmp_str(lines[l], 'END')):
                l += 1
                break
              elif (cmp_str(lines[l], 'VALE')):
                pass
              elif (cmp_str(lines[l], 'EXCH')):
                pass
              elif (cmp_str(lines[l], '1STO')):
                l += 1
              else:
                return None
              l += 1
            else:
              return None
            l -= 1
          else:
            parts = fortran_split(lines[l])
            nums = to_float(parts[1:4])
            assert (len(nums) == 3)
            n += 1
          l += 1
        else:
          return None
        if (n == 0):
          return None
      except:
        return None
    elif (name in ['CONSTRAINTS', 'NGEXCLUDE']):
      try:
        first = True
        while (l < len(lines)):
          lab = first_word(lines[l])
          l += 1
          if ((lab.upper() == 'INVERT') and (name == 'NGEXCLUDE')):
            continue
          if (cmp_str(lab, 'VALU')):
            first = False
            continue
          if (cmp_str(lab, 'END')):
            break
          if (first):
            typ = first_word(lines[l]).upper()
            assert (typ in ['BOND', 'ANGLE', 'LANGLE(1)', 'LANGLE(2)', 'DIHEDRAL', 'OUTOFP', 'DISSOC', 'CARTESIAN', 'EDIFF', 'SPHERE', 'TRANSVERSE', 'FRAGMENT'])
            l += 1
        else:
          return None
      except:
        return None
    elif (name in ['EPOT', 'EFLD', 'FLDG']):
      try:
        n = first_int(lines[l])
        l += 1
        for i in range(n):
          parts = fortran_split(lines[l])
          try:
            m = fortran_float(parts[0])
          except ValueError:
            pass
          else:
            nums = to_float(parts[0:3])
            assert (len(nums) == 3)
          l += 1
      except:
        return None
    elif (name == 'ISOTOPES'):
      try:
        n = first_int(lines[l])
        l += 1
        for i in range(n):
          parts = fortran_split(lines[l])
          if ((len(parts) > 2) and (parts[2].upper() == 'DALTON')):
            m = fortran_int(parts[0])
            m = fortran_float(parts[1])
          else:
            nums = to_int(parts[0:2])
            assert (len(nums) == 2)
          l += 1
      except:
        return None
    elif (name == 'REACTION'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts.pop(2))
        nums = to_float(parts[0:3])
        assert (len(nums) >= 2)
        l += 1
      except:
        return None
    elif (name == 'RELATIVISTIC'):
      try:
        kw = first_word(lines[l]).upper()
        assert (kw[0] == 'R')
        n = int(kw[1:3])
        assert (kw[3] in ['O', 'E', 'S', 'M', 'C'])
        l += 1
      except:
        return None
    elif (name == 'TITLE'):
      kwlist = [kw.get('NAME') for kw in keyword.getparent().xpath(kw_exp)]
      ll = 1
      # all lines (max 10) until a keyword matches
      while ((l+ll < len(lines)) and (ll < 10)):
        if (any(cmp_str(lines[l+ll], kw) for kw in kwlist)):
          break
        ll += 1
      l += ll
    elif (name == 'XBAS'):
      lab = first_word(lines[l])
      l += 1
      if ('.' in lab):
        while (l < len(lines)):
          if (cmp_str(lines[l].split()[0], 'END')):
            break
          l += 1
        else:
          return None
    elif (name == 'XFIELD'):
      nmul = [0, 1, 4, 10]
      npol = [0, 1, 6]
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        if ((len(parts) > 1) and (parts[1].upper() == 'ANGSTROM')):
          parts.pop(1)
        l += 1
        nord = 1
        np = 0
        nfrag = 0
        nread = 0
        if (len(parts) > 1):
          nums = to_int(parts[1:5])
          if (len(nums) > 0):
            nord = nums[0]
          if (len(nums) > 1):
            np = nums[1]
          if (len(nums) > 2):
            nfrag = nums[2]
          if (len(nums) > 3):
            nread = nums[3]
        nint = nfrag+nread
        nfloat = 3+nmul[nord+1]+npol[np]
        for i in range(n):
          parts = fortran_split(lines[l])
          nums = to_int(parts[0:nint])
          assert (len(nums) == nint)
          nums = to_float(parts[nint:nint+nfloat])
          assert (len(nums) == nfloat)
          l += 1
      except:
        return None
    else:
      return None

  elif (module == 'GUESSORB'):
    if (name == 'PRMO'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        assert (n in [1, 2, 3, 4])
        if (len(parts) > 1):
          n = fortran_float(parts[1])
        l += 1
      except:
        return None
    else:
      return None

  elif (module in ['GUGA', 'GUGADRT']):
    if (name == 'REFERENCE'):
      try:
        parts = fortran_split(lines[l])
        n1, n2 = to_int(parts[0:2])
        l += 1
        for i in range(n1):
          s = first_word(lines[l])
          assert ((len(s) >= n2) and all([c in '012' for c in s]))
          l += 1
      except:
        return None
    else:
      return None

  elif (module == 'LOCALISATION'):
    if (name in ['LOCN', 'LOCC']):
      try:
        parts = fortran_split(lines[l])
        n = fortran_float(parts[1])
        n = fortran_int(parts[0])
        l += 1
        parts = fortran_split(lines[l])
        assert (len(parts) >= n)
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'MBPT2'):
    if (name == 'LOVMP2'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_float(parts[1])
        n = fortran_int(parts[0])
        l += 1
        parts = fortran_split(lines[l])
        assert (len(parts) >= n)
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'POLY_ANISO'):
    if (name in ['LIN1', 'PAIR', 'LIN3', 'ALIN', 'LIN9']):
      try:
        n = first_int(lines[l])
        l += 1
        if (name in ['LIN1', 'PAIR']):
          m = 1
        elif (name in ['LIN3', 'ALIN']):
          m = 3
        elif (name == 'LIN9'):
          m = 9
        for i in range(n):
          parts = fortran_split(lines[l])
          nums = to_int(parts[0:2])
          assert (len(nums) == 2)
          nums = to_float(parts[2:2+m])
          assert (len(nums) == m)
          l += 1
      except:
        return None
    elif (name == 'HEXP'):
      try:
        n = first_int(lines[l])
        ll = test_standard(lines[l-1:], 'REALS_COMPUTED', 1)
        if ((ll is None) or (ll < 1)):
          return None
        l += ll-1
        ll = test_standard(lines[l-1:], 'REALS_COMPUTED', n+1)
        if ((ll is None) or (ll < 1)):
          return None
        l += ll-1
      except:
        return None
    elif (name == 'NNEQ'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        gv.lookup['3NNEQ'] = 3*n
        parts = [i.upper().strip() for i in parts[1:3]]
        assert (len(parts) == 2)
        assert all([i in ['T', 'F'] for i in parts])
        more = (parts[0] == 'F')
        l += 1
        parts = fortran_split(lines[l])
        nums = to_int(parts[0:n])
        assert (len(nums) == n)
        l += 1
        parts = fortran_split(lines[l])
        nums = to_int(parts[0:n])
        assert (len(nums) == n)
        l += 1
        if (more):
          parts = fortran_split(lines[l])[0:n]
          assert (len(parts) == n)
          parts = [i.upper().strip() for i in parts]
          assert all([i in ['A', 'B'] for i in parts])
          l += 1
          for i in parts:
            if (i == 'B'):
              n = fortran_float(first_word(lines[l]))
              l += 1
      except:
        return None
    elif (name == 'SYMM'):
      try:
        m = 1
        while (l < len(lines)):
          try:
            n = first_int(lines[l])
          except:
            break
          assert (n == m)
          m += 1
          l += 1
          while (l < len(lines)):
            ll = test_standard(lines[l-1:], 'REALS', 9)
            if ((ll is None) or (ll < 1)):
              break
            l += ll-1
        if (m == 1):
          return None
      except:
        return None
    else:
      return None

  elif (module == 'MCLR'):
    if (name == 'THERMO'):
      try:
        n = first_int(lines[l])
        l += 1
        n = fortran_float(first_word(lines[l]))
        l += 1
        while (l < len(lines)):
          if (cmp_str(lines[l].split()[0], 'END')):
            l += 1
            break
          n = fortran_float(first_word(lines[l]))
          l += 1
        else:
          return None
      except:
        return None
    else:
      return None

  elif (module == 'QMSTAT'):
    if (name == 'DPARAMETERS'):
      while (l < len(lines)):
        try:
          parts = fortran_split(lines[l])
          nums = to_float(parts[0:2])
          assert (len(nums) == 2)
          l += 1
        except:
          break
    elif (name == 'DISPERSION'):
      try:
        for ll in range(2):
          parts = fortran_split(lines[l])
          nums = to_float(parts[0:4])
          assert (len(nums) == 4)
          l += 1
      except:
        return None
      while (l < len(lines)):
        try:
          parts = fortran_split(lines[l])
          nums = to_float(parts[0:4])
          assert (len(nums) == 4)
          l += 1
        except:
          break
    elif (name == 'EXTERNAL'):
      try:
        n = first_int(lines[l])
        l += 1
        for i in range(n):
          parts = fortran_split_quotes(lines[l])
          m = fortran_float(parts[0])
          assert (len(parts[1]) == 8)
          m = fortran_int(parts[2])
          l += 1
      except:
        return None
    else:
      return None

  elif (module == 'QUATER'):
    if (name in ['GEO1', 'GEO2']):
      try:
        n = first_int(lines[l])
        l += 2
        for i in range(n):
          parts = fortran_split(lines[l])
          nums = to_float(parts[1:4])
          assert (len(nums) == 3)
          l += 1
        assert cmp_str(lines[l], 'END')
        l += 1
      except:
        return None
    else:
      return None

  elif (module == 'RASSCF'):
    if (name == 'CIROOT'):
      try:
        parts = fortran_split(lines[l])
        n1, n2 = to_int(parts[0:2])
        if (len(parts) > 2):
          n3 = fortran_int(parts[2])
        else:
          n3 = 0
        ll = 1
        if (n3 != 1):
          ll += test_standard(lines[l+ll-1:], 'INTS', n1) - 1
          if (n1 > 1):
            ll += test_standard(lines[l+ll-1:], 'REALS', n1) - 1
        l += ll
      except:
        return None
    elif (name == 'GASSCF'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        ll = 1
        for i in range(n):
          ll += test_standard(lines[l+ll-1:], 'INTS_LOOKUP', 'NSYM') - 1
          ll += test_standard(lines[l+ll-1:], 'INTS', 2) - 1
        l += ll
      except:
        return None
    elif (name == 'SUPSYM'):
      try:
        nsym = gv.lookup['NSYM']
        if (nsym is None):
          nsym = gv.syms[-1]
          find = True
        else:
          find = False
        for s in range(nsym):
          try:
            n = first_int(lines[l])
            l += 1
            if (n > 0):
              parts = fortran_split(lines[l])
              m = fortran_int(parts[0])
              nums = to_int(parts[1:m+1])
              assert (len(nums) == m)
              l += 1
          except:
            if (find and (s in gv.syms)):
              gv.lookup['NSYM'] = s
              break
            else:
              raise
      except:
        return None
    else:
      return None

  elif (module == 'RASSI'):
    if (name in ['NROFJOBIPHS', 'NR OF JOBIPHS']):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        gv.lookup['NROFJOBIPHS'] = n
        if ((len(parts) > 1) and (parts[1].upper().strip() == 'ALL')):
          ll = 1
        else:
          ll = test_standard(lines, 'INTS_COMPUTED', 1) - l
          parts = fortran_split(' '.join(lines[l:l+ll]))
          nums = to_int(parts[1:])
          for i in range(n):
            ll += test_standard(lines[l+ll-1:], 'INTS', nums[i]) - 1
        l += ll
      except:
        return None
    elif (name == 'PROPERTY'):
      try:
        parts = fortran_split_quotes(lines[l])
        n = fortran_int(parts.pop(0))
        l += 1
        for i in range(n):
          if (len(parts) < 2):
            parts = fortran_split_quotes(lines[l])
            l += 1
          parts.pop(0)
          m = fortran_int(parts.pop(0))
      except:
        return None
    else:
      return None

  elif (module == 'SCF'):
    if (name == 'CONSTRAINTS'):
      try:
        parts = fortran_split(lines[l])
        n = to_int(parts)
        assert (len(n) in gv.syms)
        l += 1
        for m in n:
          parts = fortran_split(lines[l])
          nums = to_int(parts[0:m])
          assert (len(nums) == m)
          assert all([i in [-1, 1] for i in nums])
          l += 1
      except:
        return None
    elif (name == 'OCCUPIED'):
      ll = test_standard(lines[l-1:], 'INTS_LOOKUP', 'NSYM')
      if (ll is not None):
        l += ll-1
        ll = test_standard(lines[l-1:], 'INTS_LOOKUP', 'NSYM')
        if (ll is not None):
          l += ll-1
      else:
        return None
    elif (name in ['OCCNUMBERS', 'MCCNUMBERS']):
      ll = 0
      while (l+ll < len(lines)):
        parts = fortran_split(lines[l+ll])
        try:
          nums = to_float(parts)
        except:
          break
        ll += 1
      l += ll
    else:
      return None

  elif (module == 'SINGLE_ANISO'):
    if (name == 'QUAX'):
      try:
        n = first_int(lines[l])
        l += 1
        if (n == 3):
          for i in range(3):
            parts = fortran_split(lines[l])
            nums = to_float(parts[0:3])
            assert (len(nums) == 3)
            l += 1
      except:
        return None
    else:
      return None

  elif (module == 'SLAPAF'):
    if (name == 'RTRN'):
      try:
        parts = fortran_split(lines[l])
        n = fortran_int(parts[0])
        n = fortran_float(parts[1])
      except:
        return None
      l += 1
    elif (name in ['INTERNAL', 'TSCONSTRAINTS']):
      try:
        first = True
        while (l < len(lines)):
          lab = first_word(lines[l])
          l += 1
          if (cmp_str(lab, 'VARY')):
            first = False
            continue
          if (cmp_str(lab, 'FIX', 3)):
            first = False
            continue
          if (cmp_str(lab, 'VALU')):
            first = False
            continue
          if (cmp_str(lab, 'END')):
            break
          if (first):
            typ = first_word(lines[l]).upper()
            assert (typ in ['BOND', 'ANGLE', 'LANGLE(1)', 'LANGLE(2)', 'DIHEDRAL', 'OUTOFP', 'DISSOC', 'CARTESIAN', 'EDIFF', 'SPHERE', 'TRANSVERSE', 'FRAGMENT'])
            l += 1
        else:
          return None
      except:
        return None
    elif (name == 'THERMO'):
      try:
        n = first_int(lines[l])
        l += 1
        n = fortran_float(first_word(lines[l]))
        l += 1
        while (l < len(lines)):
          if (cmp_str(lines[l].split()[0], 'END')):
            l += 1
            break
          n = fortran_float(first_word(lines[l]))
          l += 1
        else:
          return None
      except:
        return None
    else:
      return None

  elif (module == 'VIBROT'):
    if (name in ['POTE', 'OBSE']):
      if (name == 'OBSE'):
        l += 1
        if (l >= len(lines)):
          return None
      ll = 0
      while (l+ll < len(lines)):
        try:
          parts = fortran_split(lines[l+ll])
          nums = to_float(parts[0:2])
          assert (len(nums) == 2)
          ll += 1
        except:
          if ((name != 'OBSE') and (ll == 0)):
            return None
          l += ll
          try:
            assert cmp_str(lines[l], 'PLOT')
            l += 1
            parts = fortran_split(lines[l])
            nums = to_float(parts[0:3])
            assert (len(nums) == 3)
            l += 1
            break
          except:
            break
    elif (name == 'ATOMS'):
      try:
        parts = fortran_split(lines[l])
        l += 1
        n = fortran_int(parts[0])
        if (n < 0):
          l += 1
        n = fortran_int(parts[2])
        if (n < 0):
          l += 1
      except:
        return None
    else:
      return None

  elif (module == 'WFA'):
    if (name == 'ATLISTS'):
      try:
        n = first_int(lines[l])
        l += 1
        for i in range(n):
          parts = lines[l].split()
          assert (parts.pop() == '*')
          nums = to_int(parts)
          l += 1
      except:
        return None
    elif (name == 'PROPLIST'):
      try:
        parts = lines[l].split()
        assert (parts[-1] == '*')
        l += 1
      except:
        return None
    else:
      return None
  else:
    return None
  return l

# Read and parse an XML keyword database from a file
def read_db(filename):
  root = ET.parse(filename).getroot()

  # Copy all content from INCLUDE references
  for include in root.xpath('//INCLUDE'):
    xmod = root.find('MODULE[@NAME="{0}"]'.format(include.get('MODULE')))
    parent = include.getparent()
    for node in xmod:
      parent.append(ET.fromstring(ET.tostring(node)))

  # Add END keywords in modules and block groups if not already there
  for group in root.xpath('MODULE | //GROUP[@KIND="BLOCK"]'):
    if (group.find('KEYWORD[@NAME="END"]') is None):
      group.append(ET.Element('KEYWORD', NAME='END', KIND='SINGLE'))
    # groups can also end with ENDX (X is first letter in group's name)
    if (group.tag == 'GROUP'):
      end = 'END' + group.get('NAME')[0]
      if (group.find('KEYWORD[@NAME="{}"]'.format(end)) is None):
        group.append(ET.Element('KEYWORD', NAME=end, KIND='SINGLE'))

  return root

# This is a hack for the XYZ input of GATEWAY/SEWARD
# to hide "native input" keywords and enable "XYZ input" ones
def enable_xyz(module):
  for kw in module.xpath('KEYWORD[@EXCLUSIVE="COORD"]'):
    kw.set('NAME', '#{}'.format(kw.get('NAME')))
  for kw in module.xpath('KEYWORD[starts-with(@NAME, "*")]'):
    kw.set('NAME', kw.get('NAME')[1:])

# Validate a list of strings against a keyword database
def validate(inp, db):
  rc = 0
  result = []
  found = []

  # Re-initialize global variables
  gv.current_name = None
  gv.lookup = { 'NSYM': None }
  gv.syms = [ 1, 2, 4, 8 ]

  if (isinstance(inp, list)):
    inputlines = inp
  else:
    try:
      with open(inp, 'r') as f:
        inputlines = f.readlines()
    except:
      rc = -1
      return (rc, ['Input not found, no validation possible'])

  program = inputlines[0].strip().upper()
  if (program.startswith('&')):
    program = program[1:]
  else:
    rc = 2
    return (rc, ['No module name'])
  if (db is None):
    rc = -1
    return (rc, ['Database not found, no validation possible'])
  if (db == 'no_lxml'):
    rc = -2
    return (rc, ['lxml python module not found, no validation possible'])
  module = db.find('MODULE[@NAME="{0}"]'.format(program))
  if (module is None):
    rc = 3
    return (rc, ['No module named "{0}"'.format(program)])

  # This is a hack for the XYZ input of GATEWAY/SEWARD
  # to hide "XYZ input" keywords (they will be enabled if COORD, GROMACS or TINKER is found)
  if (module.get('NAME') in ['GATEWAY', 'SEWARD']):
    for kw in module.xpath('KEYWORD[@NAME="BASIS (XYZ)"] | KEYWORD[@NAME="GROUP"]'):
      kw.set('NAME', '*{}'.format(kw.get('NAME')))

  result.append('Input for: {0}'.format(program))

  stack = [None]
  group = module

  nline = 1
  bad = False
  # check line by line until the end of the file or the module input
  while (nline < len(inputlines) and group is not None):
    # skip blank lines (between keywords)
    if (blank(inputlines[nline])):
      nline += 1
      continue
    # check all allowed keywords in the current module/group
    anykw = False
    for kw in group.xpath(kw_exp):
      name = kw.get('NAME')
      # skip GUI-only keywords
      if (kw.get('LEVEL') == 'GUI'):
        continue
      if (kw.tag == 'KEYWORD'):
        # n should return the number of lines (if any) matching the keyword
        n = test_keyword(inputlines[nline:], kw)
        # for radio groups, the first matching keyword ends the group
        if ((group.get('KIND') == 'RADIO') and (n is not None)):
          n *= -1
        anykw = anykw or cmp_str(name, inputlines[nline])
        # if there's a match, print the line(s)
        if (n):
          if ((n < 0) and (len(stack) > 1)):
            result.append('  < leaving group {}'.format(group.get('NAME')))
          else:
            result.append('  keyword {}'.format(name))
          for i in range(abs(n)):
            result.append('    |{0}'.format(inputlines[nline+i]).strip('\n'))
          # a negative n ends the group, so pop the stack
          if (n < 0):
            group = stack.pop(-1)
          bad = False
          found.append(name)
          if ((program in ['GATEWAY', 'SEWARD']) and (name in ['COORD', 'XBAS', 'TINKER'])):
            enable_xyz(kw.getparent())
          break
      else:
        n = test_group(inputlines[nline:], kw)
        # n is 1 if a group starts
        if (n > 0):
          stack.append(group)
          group = kw
          result.append('  > entering group {}'.format(name))
          for i in range(abs(n)):
            result.append('    |{0}'.format(inputlines[nline+i].strip('\n')))
          bad = False
          found.append(name)
          if ((program in ['GATEWAY', 'SEWARD']) and (name == 'GROMACS')):
            enable_xyz(kw.getparent())
          break
        # otherwise, no line should be consumed
        elif (n < 1):
          n = 0
    # if nothing has matched, it's either an unknown keyword or syntax error
    else:
      # consume one line
      n = 1
      # print message only for first one of consecutive bad lines
      if (not bad):
        if (anykw):
          result.append('  *** syntax error at "{}" ***'.format(inputlines[nline].strip()))
        else:
          result.append('  *** unknown keyword at "{}" ***'.format(inputlines[nline].strip()))
      result.append('    *** {0}'.format(inputlines[nline].strip('\n')))
      bad = True
      rc = 1
    nline += abs(n)

  # Check for required keywords
  for kw in module.xpath('.//KEYWORD[@INPUT="REQUIRED"]'):
    # skip keywords imported from another module
    if (kw.get('MODULE') != module.get('NAME')):
      continue
    name = kw.get('NAME')
    # skip disabled keywords
    if (name[0] in ['*', '#']):
      continue
    # special case
    if ((name == 'COORD') and any([i in found for i in ['BASIS (NATIVE)', 'XBAS', 'GROMACS', 'TINKER']])):
      continue
    if (name not in found):
      result.append('*** Keyword {} is required, but was not found'.format(name))
      rc = 5

  # Check select groups
  for group in module.xpath('.//SELECT'):
    names = []
    for kw in group.xpath('.//KEYWORD'):
      if (kw.get('NAME') in found):
        names.append(kw.get('NAME'))
    if (len(names) > 1):
      result.append('*** Keywords {} are not compatible'.format(', '.join(names)))
      rc = 6

  # Check keyword-specific required/exclusive keywords
  for name in found:
    kw = module.find('.//KEYWORD[@NAME="{}"]'.format(name))
    try:
      req = kw.get('REQUIRE').split(',')
      for r in req:
        alt = r.split('.OR.')
        if (not any([i in found for i in alt])):
          result.append('*** Keyword {} is required by {}, but was not found'.format('|'.join(alt), name))
          rc = 7
    except AttributeError:
      pass
    try:
      exc = kw.get('EXCLUSIVE').split(',')
      for i in exc:
        if (i in found):
          result.append('*** Keyword {} is not compatible with {}'.format(i, name))
          rc = 8
    except AttributeError:
      pass

  # Undo possible changes to keyword names
  for kw in module.xpath('KEYWORD[starts-with(@NAME, "*")] | KEYWORD[starts-with(@NAME, "#")]'):
    kw.set('NAME', '{}'.format(kw.get('NAME')[1:]))

  return (rc, result)
