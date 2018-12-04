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

# Generate alpha and beta recursion coefficients for the Rys polynomials
# Output compatible with the "abdata" file for Molcas

from mpmath import mp
from rys_aux import *

# Set floating point precision
# Single precision (32 bits): 24-bit significand, 6 decimal places
# Double precision (64 bits): 53-bit significand, 15 decimal places
# Quadruple precision (128 bits): 113-bit significand, 33 decimal places
mp.dps = 33

tabini = -2
tabend = 678
preci = mp.mpf('1e-20')
terms = 20

# Print header
print('''Data tables from genab.py
=========================

Each table gives the recursion parameters ALPHA(K), BETA(K) and the
zeroth order polynomial P0 needed to generate Rys polynomials up to
order MAXDEG for different values of the continuous parameter T.
NTAB1 and NTAB2 are the first and last table index. PRECIS is the
accuracy to which the values were computed. For each table with index
ITAB, the corresponding T value is:

  T = (sqrt(72900-(260-ITAB)*ITAB)+ITAB-270)/10

and the inverse formula for later interpolation:

  ITAB = 5*T+200*T/(14+T)

NTAB1, NTAB2, MAXDEG:
{0} {1} {2}
PRECIS
{3}'''.format(tabini, tabend, terms-1, preci))

# Compute and print alpha and beta values for different x values
for i in range(tabini, tabend+1):
  # Get the x value from the index
  x = mp.mpf(mp.sqrt(72900-(260-i)*i)+i-270)/10
  alpha, beta = rysab(terms, x, preci=preci)
  # Transform to the Molcas convention
  beta = [mp.sqrt(b) for b in beta]
  p0 = one/beta[0]
  beta[0] = zero

  # Print the results
  print('TAB POINT NR., T VALUE, P0 VALUE:')
  print('{0} {1} {2}'.format(i, figs(x), figs(p0)))
  print('ALPHA ARRAY:')
  printlist(alpha, 4, f=figs)
  print('BETA ARRAY:')
  printlist(beta, 4, f=figs)
