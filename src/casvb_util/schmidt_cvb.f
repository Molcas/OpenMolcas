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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine schmidt_cvb(c,nvec,sao,n,metr)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension c(n,nvec),sao(*)

      if(metr.eq.0)then
        i1 = mstackr_cvb(nvec)
        call schmidt2_cvb(c,c,w(i1),nvec,sao,n,metr)
        call mfreer_cvb(i1)
      else
        i1 = mstackr_cvb(n*nvec)
        i2 = mstackr_cvb(nvec)
        call schmidt2_cvb(c,w(i1),w(i2),nvec,sao,n,metr)
        call mfreer_cvb(i1)
      endif
      return
      end
