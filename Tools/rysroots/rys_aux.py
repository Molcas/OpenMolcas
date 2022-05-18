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
#===============================================================================
# Auxiliary functions for the rysroots tools
#===============================================================================

from mpmath import mp
from decimal import Decimal

# Useful constants
zero  = mp.mpf(0)
one   = mp.mpf(1)
two   = mp.mpf(2)
four  = mp.mpf(4)
quart = one/four
pi    = mp.pi

#===============================================================================
# Wrapper function to return an mp float converted to string, with the specified
# number of decimal places. "exp" can be a boolean to force exponential notation
# or an integer to specify the minimum number of digits in the exponent.
def figs(x, n=20, exp=0, sign=False):
  if (sign and x >= 0.0):
    s = ' '
  else:
    s = ''
  # mp always returns '0.0' for zero. Bypass it using standard formats
  if (x == 0.0):
    fmt = '{{0:.{0}f}}'.format(n)
    # python uses two-digit exponent, but mp uses the least amount. Hard code it for zero
    if (exp):
      fmt += 'e+'+'0'*exp
    return s+fmt.format(0.0)
  # Convert significant figures to decimal places:
  # With exponent, integer part should be 1 figure, so it's n+1
  # Without exponent, floor(log10(abs(x)))+1 gives the number of integer figures (except for zero)
  if (exp):
    tmp = s+mp.nstr(mp.mpf(x), n+1, strip_zeros=False, min_fixed=mp.inf, max_fixed=-mp.inf, show_zero_exponent=True)
    # Add zeros to the exponent if necessary
    exppos = tmp.find('e')+2
    diff = exp-(len(tmp)-exppos)
    if (diff > 0):
      tmp = tmp[:exppos] + '0'*diff + tmp[exppos:]
    return tmp
  else:
    return s+mp.nstr(mp.mpf(x), int(mp.floor(mp.log10(abs(x))))+n+1, strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)
#===============================================================================

#===============================================================================
# Function to print a list in chunks of n elements
# The argument f is a function that converts each element to string
def printlist(l, n, f=str):
  for i in range(0, len(l), n):
    num = min(n, len(l)-i)
    print(' '.join(map(f, l[i:i+num])))
#===============================================================================

#===============================================================================
# range analogue for real numbers (with exact short decimal representation)
def frange(start, end, step):
  # Use a parallel Decimal representation to keep operations exact,
  # but yield the original representation.
  # This could mean that the final value is slightly larger than "end"
  a = Decimal(str(start))
  b = Decimal(str(end))
  c = Decimal(str(step))
  if (step < 0.0):
    f = -1
  else:
    f = 1
  x,d = (start,a)
  i = 0
  while (f*d <= f*b):
    yield x
    i += 1
    x,d = (start+i*step, a+i*c)
#===============================================================================

#===============================================================================
# Given alpha and beta recursion coefficients for a set of orthogonal
# polynomials, compute the corresponding roots and weights for a Gaussian
# quadrature. The roots and weights are obtained from the eigenvalues and
# eigenvectors of a tridiagonal matrix containing the alpha and beta
# coefficients (Jacobi matrix).
#
# Code translated directly from:
#   Walter Gautschi. ACM Trans. Math. Softw. 20 (1994) 21-62
#   doi:10.1145/174603.174605
# (routine supplied by G. H. Golub according to a footnote)
def gauss(alpha, beta, eps, nloop=30):
  n = len(alpha)
  if (len(beta) != n):
    raise

  # initialize
  root = alpha[:]
  weight = [one]+[zero for i in range(n-1)]
  aux = list(map(mp.sqrt, beta[1:]))
  aux.append(zero)

  for l in range(n):
    for j in range(nloop):
      for m in range(l, n):
        if (m == n-1):
          break
        if (abs(aux[m]) <= eps*(abs(root[m])+abs(root[m+1]))):
          break
      dp = root[l]
      if (m == l):
        break
      dg = (root[l+1]-dp)/(two*aux[l])
      dr = mp.sqrt(dg*dg+one)
      if (mp.sign(dg) < 0):
        dr *= -1
      dg = root[m]-dp+aux[l]/(dg+dr)
      ds = one
      dc = one
      dp = zero
      mml = m-l
      for ii in range(1, mml+1):
        i = m-ii
        df = ds*aux[i]
        db = dc*aux[i]
        if (abs(df) < abs(dg)):
          ds = df/dg
          dr = mp.sqrt(ds*ds+one)
          aux[i+1] = dg*dr
          dc = one/dr
          ds = ds*dc
        else:
          dc = dg/df
          dr = mp.sqrt(dc*dc+one)
          aux[i+1] = df*dr
          ds = one/dr
          dc = dc*ds
        dg = root[i+1]-dp
        dr = (root[i]-dg)*ds+two*dc*db
        dp = ds*dr
        root[i+1] = dg+dp
        dg = dc*dr-db
        df = weight[i+1]
        weight[i+1] = ds*weight[i]+dc*df
        weight[i] = dc*weight[i]-ds*df
      root[l] = root[l]-dp
      aux[l] = dg
      aux[m] = zero
    if (j >= nloop-1):
      raise

  # scale, sort and return
  return (list(t) for t in zip(*sorted(zip(root, [beta[0]*w*w for w in weight]))))
#===============================================================================

#===============================================================================
# Generate alpha and beta recursion coefficients according to the discretized
# Stieltjes procedure for the orthogonal polynomials with an inner product
# given by the input roots and weights
#
# Code translated directly from:
#   Walter Gautschi. ACM Trans. Math. Softw. 20 (1994) 21-62
#   doi:10.1145/174603.174605
def stieltjes(n, roots, weights):
  naux = len(roots)
  if (len(weights) != naux):
    raise
  if (n > naux):
    raise

  huge = mp.mpf('1.0e+60')
  tiny = mp.mpf('1.0e-60')
  beta = [sum(weights)]
  alpha = [sum([x*w for x,w in zip(roots, weights)])/beta[0]]

  p0 = [zero for i in range(naux)]
  p1 = [zero for i in range(naux)]
  p2 = [one for i in range(naux)]
  sum0 = beta[0]
  for k in range(1, n):
    sum1 = zero
    sum2 = zero
    for m in range(naux):
      if (weights[m] == zero):
        continue
      p0[m] = p1[m]
      p1[m] = p2[m]
      p2[m] = (roots[m]-alpha[k-1])*p1[m]-beta[k-1]*p0[m]
      if ((abs(p2[m]) > huge) or (abs(sum2) > huge)):
        raise
      t = weights[m]*p2[m]*p2[m]
      sum1 = sum1+t
      sum2 = sum2+t*roots[m]
    if (abs(sum1) < tiny):
      raise
    alpha.append(sum2/sum1)
    beta.append(sum1/sum0)
    sum0 = sum1

  return alpha, beta
#===============================================================================

#===============================================================================
# Generate alpha and beta recursion coefficients according to the Lanczos
# orthogonal reduction method. This is an alternative (more robust, slower)
# to Stieltjes's procedure
#
# Code translated directly from:
#   Walter Gautschi. ACM Trans. Math. Softw. 20 (1994) 21-62
#   doi:10.1145/174603.174605
def lanczos(n, roots, weights):
  naux = len(roots)
  if (len(weights) != naux):
    raise
  if (n > naux):
    raise

  p0 = roots[:]
  p1 = [zero for i in range(naux)]
  p1[0] = weights[0]
  for i in range(1, naux):
    gam = one
    sig = zero
    t = zero
    p_i = weights[i]
    xlam = roots[i]
    for k in range(i+1):
      rho = p1[k]+p_i
      tmp = gam*rho
      tsig = sig
      if (rho <= zero):
        gam = one
        sig = zero
      else:
        gam = p1[k]/rho
        sig = p_i/rho
      tk = sig*(p0[k]-xlam)-gam*t
      p0[k] = p0[k]-(tk-t)
      t = tk
      if (sig <= zero):
        p_i = tsig*p1[k]
      else:
        p_i = t*t/sig
      tsig = sig
      p1[k] = tmp

  return p0[0:n], p1[0:n]
#===============================================================================

#===============================================================================
# Compute alpha and beta three-term recursion parameters for Rys polynomials.
# An auxiliary shifted Legendre quadrature is used, for which roots and weights
# are computed the first time this is run. Note that naux=300 is most probably
# overkill.
def rysab(n, x, naux=300, preci=mp.mpf('1e-20'), method='s'):
  # Generate roots and weights for the auxiliary quadrature (shifted Legendre)
  #   alpha(i) = 1/2
  #   beta(0)  = 1
  #   beta(i)  = 1 / (4 - 1/i^2)
  # Do this only the first time (results saved as attributes) or if
  # the "naux" or "preci" parameters change
  if ((rysab.roots_aux == None) or
      (rysab.naux != naux) or (rysab.preci != preci)):
    rysab.naux = naux
    rysab.preci = preci
    alpha_aux = [mp.mpf('0.5') for i in range(naux)]
    beta_aux = [one]
    for i in range(1, naux):
      beta_aux.append(quart/(four-one/mp.mpf(i*i)))
    rysab.roots_aux, rysab.weights_aux = gauss(alpha_aux, beta_aux, preci)

  rtsgau = [r*r for r in rysab.roots_aux]
  wtsgau = [w*mp.exp(-x*r) for r,w in zip(rtsgau, rysab.weights_aux)]
  if (method == 's'):
    return stieltjes(n, rtsgau, wtsgau)
  elif (method == 'l'):
    return lanczos(n, rtsgau, wtsgau)
  else:
    raise
rysab.roots_aux = None
#===============================================================================

#===============================================================================
# Compute roots and weights for Rys polynomials
def rysroots(n, x, preci=mp.mpf('1e-20'), **kwargs):
  kwargs.update({'preci': preci})
  alpha, beta = rysab(n, x, **kwargs)
  return tuple(gauss(alpha, beta, preci))
#===============================================================================
