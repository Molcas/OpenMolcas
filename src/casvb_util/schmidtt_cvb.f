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
      subroutine schmidtt_cvb(c,nvec,t,nt,sao,n,metr)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension c(n,nvec),t(nt,nvec),sao(*)

      if(metr.eq.0)then
        call schmidtt2_cvb(c,c,nvec,t,nt,sao,n,metr)
      else
        i1 = mstackr_cvb(n*nvec)
        call schmidtt2_cvb(c,w(i1),nvec,t,nt,sao,n,metr)
        call mfreer_cvb(i1)
      endif
      return
      end
