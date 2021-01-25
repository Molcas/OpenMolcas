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
      subroutine realz_cvb(arr,nmax,nread,ifc)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension arr(nmax)
      save nbuf
      data nbuf/100/

      i1 = mstackr_cvb(nbuf)
      nread=0
100   call fzero(w(i1),nbuf)
      nleft=nmax-nread
      call real_cvb(w(i1),min(nbuf,nleft),nread1,ifc)
      call fmove_cvb(w(i1),arr(1+nread),nread1)
      nread=nread+nread1
      if(nread1.gt.0)goto 100
      call mfreer_cvb(i1)
      return
      end
