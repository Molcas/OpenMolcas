#!/usr/bin/env python3

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
# Copyright (C) 2021, Ignacio Fdez. Galván                             *
#***********************************************************************

# Help script for reformatting a basis set file in Molcas format.
# It writes the numbers in a consistent format and fails if precision
# is lost.
# Use with compact=True and remove_comments=True to facilitate comparison
# between different versions.

import sys
import re
import numpy as np
from decimal import Decimal

# Output formats: (total_places, decimal_places)
#           if compact, just use minimal space
compact  =  False
#           format for basis set exponents
exp_form =  [20, 8]
#           format for basis set coefficients, orbital energies and Fock operator
coe_form =  [15, 12]
#           format for pseudopotential parameters
pp_form =   [14, 8]
#           format for "M1" and "M2" parameters
ecpm_form = [15, 12]
#           format for "PROJOP" exponents
ecpe_form = [20, 8]
#           format for "PROJOP" coefficients and constants
ecpc_form = [15, 12]
#           maximum line length
max_line =  180

uncontracted = False
remove_comments = False
insert_comments = False

l_lab = 'spdfghi'

if compact:
  exp_form[0] = None
  coe_form[0] = None
  pp_form[0] = None
  ecpm_form[0] = None
  ecpe_form[0] = None
  ecpc_form[0] = None

# Read next line, ignores (and prints unchanged) comments and empty lines
def nextline(f):
  global uncontracted
  line = f.readline()
  while line and (line.startswith('*') or line.startswith('#') or line.startswith('!') or line.startswith('\n')):
    if re.match(r'#Contraction *UNC', line, re.I):
      uncontracted = True
    if not remove_comments:
      print(line.rstrip())
    line = f.readline()
  if line:
    return line.rstrip()
  else:
    return None

def count_items(parts):
  n = 0
  for i in parts:
    if i.startswith('*'):
      break
    n += 1
  return n

def read_nums(f, n):
  nums = []
  while len(nums) < n:
    line = nextline(f)
    nums.extend(line.split())
    if nums[-1] == '/':
      l = n-len(nums)+1
      nums = nums[:-1]+l*['0']
  assert(len(nums) == n)
  return nums

def print_line(nums, form):
  parts = [f' {format_num(i, *form)}' for i in nums]
  if max_line is None:
    print(''.join(parts))
  else:
    line = ''
    for i in parts:
      if len(line) + len(i) > max_line:
        print(line)
        line = ''
      line += i
    print(line)

################################################################################
# Set of functions to "pretty print" a number in some reserved space,
# ensuring that no precision is lost.
# Quite hackish and inefficient, but as long as it works...

# Return the "scale" of a number, i.e, the exponent of its leading digit
def scale(num):
  if num == 0:
    return 0
  return num.adjusted()

# Significant places of a number
def sig_places(num):
  digits = list(num.as_tuple().digits)
  while len(digits) > 1:
    if digits[-1] != 0:
      break
    digits.pop()
  return len(digits)

# Try to fit a number in n total places with m decimals.
# Returns a string of asterisks if it fails
def fit_num(f, n, m):
  s = sig_places(f)
  e = scale(f)
  sg = 1 if f < 0 else 0
  if (s-e-1 <= m) and (e+1 <= n-m-sg-1):
    out = f'{f:{n}.{m}f}'
  else:
    out = '*'*n
  return out

# Try to fit a number in n total places with m decimals, with a exponent (shift)
# Returns a string of asterisks if it fails
def print_shifted(f, n, m, shift=0):
  fs = f.scaleb(shift)
  if shift == 0:
    es = 0
  else:
    es = np.int(np.floor(np.log10(abs(shift))))+2
    if shift > 0:
      es += 1
  s = sig_places(fs)
  if shift == 0:
    out = fit_num(fs, n, m)
  else:
    out = fit_num(fs, n-es, m-es)
    out += f'e{-shift}'
  if len(out) > n or out[0] == '*':
    out = '*'*n
  return out

# Simplify a number that may come from an "overprecise" float,
# i.e. round to fewer decimals as long as it represents the same float number
def simplify(f):
  out = f
  ref = float(f)
  s, n, e = f.as_tuple()
  m = len(n)
  for i in range(m):
    a = f.quantize(f.scaleb(i+1))
    if float(a) != ref:
      break
    out = a
  if out != f:
    print(f'simplify {f} to {out}', file=sys.stderr)
  return out

# Regexp to fix exponents from Fortran prints, which may have D instead of E
fortfixexp = re.compile(r'([\d.])[dD]?(((?<=[dD])[+-]?|[+-])\d)')

# Format a number using the minumum necessary space
# (but still using at least one explicit zero before/after decimal point)
def format_minimal(num):
  t = num.as_tuple()
  d = [i for i in t.digits]
  e = t.exponent
  # remove trailing zeros
  while (len(d) > 1) and d[-1] == 0:
    d.pop()
    e += 1
  # "text" contains all significant digits
  text = ''.join([f'{i}' for i in d])
  # special case for zero
  if text == '0':
    return text
  # large number, use exponent if large enough
  if ((len(d) > 1) and (e > 3)) or (e > 4):
    if len(d) == 1:
      text += f'.0e{e}'
    else:
      s = len(d)-1
      # prefer fewer digits in exponent
      if (e+s > 9) and (e+1 < 10):
        s = 9-e 
      elif (e+s > 99) and (e+1 < 100):
        s = 99-e
      e_ = e+s
      text = text[0:-s] + '.' + text[-s:] + f'e{e_}'
  # not too large, not too small, print without exponent
  elif e >= 0:
    text += e*'0'
  elif len(d)+e > 0:
    text = text[:e] + '.' + text[e:]
  elif ((len(d) > 1) and (-e-len(d) < 3)) or (e > -5):
    text = '0.' + (-len(d)-e)*'0' + text
  # small number, use exponent
  else:
    if len(d) == 1:
      s = -2
      e_ = e
      z = '0'
    else:
      s = len(d)-1
      e_ = e+s
      z = ''
    text = text[0:-s] + '.' + text[-s:] + z + f'e{e_}'
  out = '-' if t.sign else ''
  out += text
  return out

# Format a decimal number in n total places with m decimals, adjusting
# the exponent if necessary.
# If n = None, formats the number with the minimal space needed.
# Raises exception if failed.
def format_num(num, n=None, m=None):
  num = fortfixexp.sub(r'\1e\2', num)
  f = Decimal(num)
  f = simplify(f)
  if n is None:
    out = format_minimal(f)
  else:
    out = print_shifted(f, n, m)
    # if it doesn't fit without exponent, try shifting until it fits
    if out[0] == '*':
      e = scale(f)
      if e > 0:
        for i in range(-e, 0):
          out = print_shifted(f, n, m, shift=i)
          if out[0] != '*':
            break
      else:
        for i in range(-e, -e+(n-m-1)):
          out = print_shifted(f, n, m, shift=i)
          if out[0] != '*':
            break
  if out[0] == '*':
    raise Exception(f'Cannot fit number in ({n}.{m}): {num}')
  # ensure that the final number represents the same float
  if float(num) != float(out):
    raise Exception(f'Changed value in ({n}.{m}): {num} -> {out}')
  return out
################################################################################

with open(sys.argv[1]) as f:
  line = nextline(f)
  while line is not None:
    if line.startswith('/'):
      # basis label
      print(line)
      parts = line.split('.')
      prim = parts[3]
      prim = re.sub(r'[spdfghi]', ' ', prim).split()
      n_p = len(prim)-1
      n_prim = [int(i) for i in prim]
      func = parts[4]
      func = re.sub(r'[spdfghi]', ' ', func).split()
      n_f = [int(i) for i in func]
      if uncontracted:
        assert(np.all(prim == func))
      # two reference lines
      line = nextline(f)
      print(line)
      line = nextline(f)
      print(line)
      options = ''
      line = nextline(f)
    elif re.match(r'Options', line, re.I):
      # all lines until "EndOptions"
      print(line)
      line = nextline(f)
      while line is not None and not re.match(r'End([ O]|$)', line, re.I):
        print(line)
        if re.match(r'FockOperator', line, re.I):
          options += 'f'
        if re.match(r'OrbitalEnergies', line, re.I):
          options += 'e'
        line = nextline(f)
      print(line)
      line = nextline(f)
    elif re.match(r'Spectral', line, re.I):
      # all lines until "End of Spectral"
      print(line)
      line = nextline(f)
      while line is not None and not re.match(r'End([ o]|$)', line, re.I):
        print(line)
        if re.match(r'Exte', line, re.I):
          line = nextline(f)
          maxl = int(line)
          print(f' {maxl:4d}')
          for i in range(maxl+1):
            line = nextline(f)
            mp, = [int(j) for j in line.split()]
            print(f' {mp:4d}')
            exps = read_nums(f, mp)
            for e in exps:
              print(format_num(e, *exp_form))
        elif line == 'SOC':
          line = nextline(f)
          maxl = int(line)
          print(f' {maxl:4d}')
          for i in range(maxl+1):
            line = nextline(f)
            mp, mo, md = [int(j) for j in line.split()]
            if insert_comments:
              print(f'* {l_lab[i]}-type functions')
            print(f' {mp:4d} {mo:4d} {md:4d}')
            exps = read_nums(f, mp)
            for e in exps:
              print(format_num(e, *exp_form))
            for c in range(mp):
              coefs = read_nums(f, mo)
              print_line(coefs, coe_form)
        line = nextline(f)
      print(line)
      line = nextline(f)
    elif re.match(r'PP', line, re.I):
      parts = [i.strip() for i in line.strip(';').split(',')]
      parts[0] = 'PP'
      parts[1] = parts[1][0].upper() + parts[1][1:].lower()
      maxl = int(parts[3])
      print(', '.join(parts)+' ;')
      for nl in range(maxl+1):
        line = nextline(f).split(';')[0]
        n = int(line)
        print(f'{n:3d} ;')
        for m in range(n):
          line = nextline(f).strip(';')
          parts = line.split(',')
          i = int(parts[0])
          f1 = parts[1]
          f2 = parts[2]
          print(f'{i:2d}, {format_num(f1, *pp_form)}, {format_num(f2, *pp_form)} ;')
      line = nextline(f)
    elif re.match(r'M[12]', line, re.I):
      print(line.strip().upper())
      line = nextline(f)
      n = int(line)
      print(f'{n:3d}')
      if n > 0:
        for i in [1, 2]:
          nums = read_nums(f, n)
          print_line(nums, ecpm_form)
      line = nextline(f)
    elif re.match(r'COREREP', line, re.I):
      print('COREREP')
      line = nextline(f)
      print(line)
      line = nextline(f)
    elif re.match(r'PROJOP', line, re.I):
      print('PROJOP')
      line = nextline(f)
      maxl = int(line)
      print(f' {maxl:4d}')
      for i in range(maxl+1):
        line = nextline(f)
        parts = [int(j) for j in line.split()]
        mp = parts[0]
        mo = parts[1]
        occ = parts[2:]
        assert(len(occ) <= mo)
        if len(occ) == mo:
          if all([j == 4*i+2 for j in occ]):
            occ = []
        print(f' {mp:4d} {mo:4d}' + ''.join([f' {j:4d}' for j in occ]))
        coefs = read_nums(f, mo)
        print_line(coefs, ecpc_form)
        exps = read_nums(f, mp)
        for e in exps:
          print(format_num(e, *ecpe_form))
        for c in range(mp):
          coefs = read_nums(f, mo)
          print_line(coefs, ecpc_form)
      line = nextline(f)
    else:
      # here go the exponents and coefficients for the basis set
      parts = line.split()
      Z = float(parts[0])
      maxl = n_p
      print(f'   {Z:5.1f} {maxl:3d}')
      line = nextline(f)
      for nl in range(maxl+1):
        if insert_comments:
          print(f'* {l_lab[nl]}-type functions')
        parts = line.split()
        n_p = int(parts[0])
        assert(n_prim[nl] == n_p)
        assert(count_items(parts) <= 2)
        if count_items(parts) > 1:
          assert(n_f[nl] == int(parts[1]))
        print(f' {n_p:4d} {n_f[nl]:4d}')
        exps = read_nums(f, n_p)
        for e in exps:
          print(format_num(e, *exp_form))
        if not uncontracted:
          for c in range(n_p):
            coefs = read_nums(f, n_f[nl])
            print_line(coefs, coe_form)
        line = nextline(f)
        if 'f' in options:
          n = int(line)
          print(f' {n:4d}')
          for c in range(n):
            coefs = read_nums(f, n)
            print_line(coefs, coe_form)
          line = nextline(f)
        if 'e' in options:
          n = int(line)
          print(f' {n:4d}')
          coefs = read_nums(f, n)
          if n:
            print_line(coefs, coe_form)
          line = nextline(f)
