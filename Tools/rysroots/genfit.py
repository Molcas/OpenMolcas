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
# Copyright (C) 2017, Ignacio Fdez. Galván                             *
#***********************************************************************

# Compute fitting polynomials for Rys roots and weights.
# The TMax values come from the asymptotic.py script.
# For each order, the range [0,TMax] is split in intervals and for each
# interval a set of fitting polynomials is generated that approximate the
# roots and weights within the interval with the desired accuracy.
# Each interval is expanded as much as possible to keep the fitting accuracy.
# The final part of the output is compatible with the "rysrw" file for Molcas.

from mpmath import mp
from rys_aux import *
import pickle

# Set floating point precision
# Single precision (32 bits): 24-bit significand, 6 decimal places
# Double precision (64 bits): 53-bit significand, 15 decimal places
# Quadruple precision (128 bits): 113-bit significand, 33 decimal places
mp.dps = 33

mRys = 9
fit_degree = 6

order = list(range(1, mRys+1))
only_n = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# Requested accuracy (input)
accuracy = list(map(mp.mpf, ['1e-16', '1e-16', '1e-15', '1e-15', '1e-14', '1e-14', '1e-13', '1e-13', '1e-12']))
# x step (input)
dx = list(map(mp.mpf, ['0.008', '0.007', '0.009', '0.009', '0.012', '0.012', '0.017', '0.016', '0.023']))
fine_step = mp.mpf('0.001')
# x value from which the asymptotic formula gives a relative error in roots and weights lower than the requested accuracy (computed with 0.01 precision)
TMax = list(map(mp.mpf, ['38.80', '46.21', '50.71', '57.03', '60.32', '66.01', '68.69', '74.03', '76.26']))

# This would be with the same accuracy at all orders
#accuracy = [mp.mpf('1e-16') for o in order]
#TMax = list(map(mp.mpf, ['38.80', '46.21', '53.24', '59.64', '65.69', '71.51', '77.16', '82.68', '88.09']))

# Number of data points before reaching "TMax"
nMap = [int(mp.floor(a/b)+1) for a,b in zip(TMax, dx)]
# x values for the data points
xlists = [[] for x in range(len(order))]
for i in range(len(order)):
  for j in range(nMap[i]):
    prev = j*dx[i]
    xlists[i].append(list(frange(prev, (j+1)*dx[i], fine_step)))

# reload the data from previous runs
try:
  with open('tmpfile', 'rb') as f:
    data = pickle.load(f)
except:
  data = {}

for n,acc in zip(order, accuracy):
  if (n not in only_n):
    continue
  xlist = []
  for l in xlists[n-1]:
    xlist.extend(l)

  # Compute roots and weights for different values of x
  # Note that there may be some "duplication" if the end of one interval
  # does not exactly match the beginning of the next (which may occur with frange)
  rtswts = {}
  for x in xlist:
    if (x in rtswts):
      continue
    roots, weights = rysroots(n, x)
    rtswts[x] = [roots, weights]
    print('x = {0}'.format(figs(x)))
    for i in range(len(roots)):
      print('{0}: {1} {2}'.format(i+1, figs(mp.sqrt(roots[i])), figs(weights[i], exp=True)))
  xlist = xlists[n-1]

  # Now proceed with the fitting
  l_min = int(mp.ceil((fit_degree+1)*fine_step/dx[n-1]))+1
  nx = len(xlist)
  num = []
  coeff = []
  center = []
  ini = 0
  print()
  print('{0} points'.format(nx))
  l_prev = 0
  while (ini < nx-1):
    # TODO: bad luck if there are not enough points left to fit
    if (ini+l_min > nx):
      raise
    # l is the number of data points to fit (minimum is l_min)
    # ini is the initial point to include in the fit, so the range is ini to ini+l-1 (ini:ini+l for python)
    l_max = nx-ini
    # start with the same l that worked last time
    l = min(max(l_min, l_prev), l_max)
    while (True):
      print('{0} to {1}'.format(ini, ini+l-1))
      flat_x = [xlist[ini][0]]
      flat_x.extend([x for l in xlist[ini:ini+l] for x in l[1:]])
      x0 = (flat_x[0]+flat_x[-1])/two
      xoff = [x-x0 for x in flat_x]
      mat = []
      for x in xoff:
        mat.append([x**i for i in range(fit_degree, -1, -1)])
      mat = mp.matrix(mat)
      max_error = zero
      cff = []
      # fit for each root (nn) and for roots (i=0) and weights (i=1)
      for nn in range(n):
        for i in range(2):
          val = [rtswts[x][i][nn] for x in flat_x]
          try:
            fit, _ = mp.qr_solve(mat, val)
          except:
            fit = mp.lu_solve(mat, val)
          cff.append(fit.T.tolist()[0])
          error = max([abs((mp.polyval(fit, x)-v)/v) for x,v in zip(xoff, val)])
          max_error = max(max_error, error)
      # if the error is too large, the number of poinst in the fit may have to be reduced
      if (max_error > acc):
        # error if it's already the minimum
        if (l == l_min):
          raise
        # if the last good fit was just one point less, that's the final one
        if (l-l_save == 1):
          l = l_save
          break
        # otherwise reduce to midpoint between current and last good fit
        l_max = l-1
        l = min(int(mp.ceil((l+l_save)/two)), l-1)
      # we have to make sure the saved data correspond to the last good fit
      else:
        cff_save = cff[:]
        x0_save = x0
        l_save = l
        # if l is already at maximum, we're done, otherwise be conservative and step up by 1
        if (l == nx-ini):
          break
        l +=1
    print('fit from {0} to {1}'.format(ini, ini+l-1))
    num.append(l_save)
    ini += l_save
    l_prev = l_save
    coeff.append(cff_save)
    center.append(x0_save)

  # the nesting order for printout is:
  #   roots, weights
  #     0th...6th order coefficients
  #       root 1...n
  #         polynomial 1...npol
  # while in coeff we have:
  #   polynomial 1...npol
  #     root + weight 1...n
  #       6th...0th order coefficients
  coeflist = []
  for i in range(2):
    for j in range(fit_degree, -1, -1):
      for k in range(n):
          for l in range(len(num)):
            coeflist.append(coeff[l][2*k+i][j])
  # map specifying which approximating polynomial to apply for each x
  Map = []
  for i in range(len(num)):
    Map.extend([i+1]*num[i])

  data[n] = {'nx': nx, 'npol': len(num), 'Map': Map, 'center': center, 'coeflist': coeflist}
  # dump the data in a temporary file, so we can run this in several steps
  with open('tmpfile', 'wb') as f:
    pickle.dump(data, f)

# Finally print the data we have available so far
# This is the part that should go in the rysrw file
print()
print(mRys, fit_degree)
printlist(accuracy, 3, f=lambda x: figs(x, n=6, exp=True))
printlist(TMax, 3, f=lambda x: figs(x, n=6))
printlist(dx, 3, f=lambda x: figs(x, n=6))
printlist([i[1]['nx'] for i in sorted(data.items())], 3)
printlist([i[1]['npol'] for i in sorted(data.items())], 3)
for i in sorted(data.items()):
  printlist(i[1]['Map'], 10)
  printlist(i[1]['center'], 6, f=lambda x: figs(x, n=10, exp=True))
for i in sorted(data.items()):
  printlist(i[1]['coeflist'], 4, f=lambda x: figs(x, exp=2, sign=True))
