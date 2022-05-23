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
      subroutine schmidtd_cvb(c1,nvec1,c2,nvec2,sao,n,metr)
c  Orthogonalize nvec2 vectors in C2 on nvec1 vectors in C1.
c  C1 vectors assumed to be orthonormal.
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension c1(n,nvec1),c2(n,nvec2),sao(*)

      if(metr.eq.0)then
        call schmidtd2_cvb(c1,c1,nvec1,c2,nvec2,n)
      else
        i1 = mstackr_cvb(n*nvec1)
        call saoon_cvb(c1,w(i1),nvec1,sao,n,metr)
        call schmidtd2_cvb(c1,w(i1),nvec1,c2,nvec2,n)
        call mfreer_cvb(i1)
      endif
      return
      end
