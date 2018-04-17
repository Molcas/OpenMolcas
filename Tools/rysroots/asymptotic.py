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

# Compute the minimum parameter x for which the roots and weights of Rys polynomials
# can be obtained with the desired accuracy from the asymptotic formula.
# The list of orders for Rys polynomials is specified in "order", the desired
# accuracy in "accuracy", and the precision for the result is the last element in "step".

from mpmath import mp
from rys_aux import *

# Set floating point precision
# Single precision (32 bits): 24-bit significand, 6 decimal places
# Double precision (64 bits): 53-bit significand, 15 decimal places
# Quadruple precision (128 bits): 113-bit significand, 33 decimal places
mp.dps = 33

# Order of the Rys polynomials (1 to 9)
order = list(range(1, 10))
# Requested accuracy (input)
# Higher orders will typically have a smaller contribution, so the accuracy can be lower
accuracy = list(map(mp.mpf, ['1e-16', '1e-16', '1e-15', '1e-15', '1e-14', '1e-14', '1e-13', '1e-13', '1e-12']))
# Precision for the limits (several passes)
step = list(map(mp.mpf, ['10', '1', '0.01']))

# Uncomment to get same accuracy all along
#accuracy = [mp.mpf('1e-16') for o in order]

# Do several passes, increasing x by the specified step and stopping when
# the error is below the requested accuracy
# Start from zero and then use the result from the previous pass
init = [mp.mpf(0) for i in order]
for s in step:
  asymp = []
  for n,x0,acc in zip(order, init, accuracy):

    print('Starting n={0}'.format(n))

    # Generate roots and weights for the asymptotic formula (Hermite)
    #   alpha(i) = 0
    #   beta(0)  = sqrt(pi)
    #   beta(i)  = i/2
    alpha_herm = [zero for i in range(2*n)]
    beta_herm = [mp.sqrt(pi)]
    for i in range(1, 2*n):
      beta_herm.append(i/two)
    roots_herm, weights_herm = gauss(alpha_herm, beta_herm, mp.mpf('1e-20'))

    i = 0
    while True:
      x = x0+i*s
      i += 1

      # Generate roots and weights of Rys polynomial of order n for parameter x
      roots, weights = rysroots(n, x)

      # Compute asymptotic (x -> inf) Rys roots and weights
      # from the second half of the Hermite solution for degree 2n
      #   r_rys = r_h^2 / x
      #   w_rys = w_h / sqrt(x)
      y = one
      if (x != 0):
        y = one/x
      roots_asymp = [r*r*y for r in roots_herm[n:]]
      weights_asymp = [w*mp.sqrt(y) for w in weights_herm[n:]]

      # Compute maximum error in roots and weights, and stop when below accuracy
      error_roots = max([abs(r1-r2)/r2 for r1,r2 in zip(roots_asymp, roots)])
      error_weights = max([abs(w1-w2)/w2 for w1,w2 in zip(weights_asymp, weights)])
      error_tot = max(error_roots, error_weights)
      print('x = {0}: error = {1}'.format(figs(x, 2), figs(error_tot, 4, exp=True)))
      if (error_tot < acc):
        asymp.append(x)
        break

  # Take a step back for the next pass
  init = [a-s for a in asymp]

  print()
  print('Step: {0}'.format(mp.nstr(s)))
  for n,a in zip(order, asymp):
    print('n = {0}: x = {1}'.format(n, figs(a, 2)))
