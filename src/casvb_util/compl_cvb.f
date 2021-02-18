************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine compl_cvb(a,nvec,n)
c  Creates orthogonal complement.
c  On entry : A is square (NxN) and contains NVEC vectors.
c  On exit  : A is a full matrix, NVEC first vectors are untouched,
c  remaining orthonormal vectors span the orthogonal complement.
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n,n)

      i1 = mstackr_cvb(n*(nvec+n))
      i2 = mstackr_cvb(n*n)
      i3 = mstackr_cvb(n)
      call compl2_cvb(a,nvec,n,w(i1),w(i2),w(i3))
      call mfreer_cvb(i1)
      return
      end
