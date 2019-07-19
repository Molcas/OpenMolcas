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
      subroutine extract_cvb(c,t,nvec1,nextract,mode,thr,s,nbf,metr)
c  MODE: 0,1 => determine NEXTRACT based on THR
c        1,3 => transform also T matrix
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension c(nbf,nvec1),t(nbf,nvec1),s(*)

      nvec=nvec1
      i1 = mstackr_cvb(nvec)
      i2 = mstackr_cvb(nbf*nvec)
      i3 = mstackr_cvb(nvec*nvec)
      i4 = mstackr_cvb(nbf*nvec)
      call svd_cvb(c,w(i1),w(i2),w(i3),nbf,nvec)

      if(mode.lt.2)then
c  Determine NEXTRACT :
        do 100 i=nvec,1,-1
100     if(w(i+i1-1).le.thr)goto 200
200     nextract=nvec-i
      endif
      call fmove_cvb(w(1+(nvec-nextract)*nbf+i2-1),c,nbf*nextract)
      if(mod(mode,2).eq.1)then
c  Apply same transformation to T :
        call mxatb_cvb(t,w(i3),nbf,nvec,nvec,w(i4))
        do 300 i=1,nvec
300     call dscal_(nbf,1.d0/w(i+i1-1),w((i-1)*nbf+i4),1)
        call fmove_cvb(w(1+(nvec-nextract)*nbf+i4-1),t,nbf*nextract)
        call schmidtt_cvb(c,nextract,t,nbf,s,nbf,metr)
      endif
      call mfreer_cvb(i1)
      return
      end
