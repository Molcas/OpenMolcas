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
      subroutine bspset2_cvb(ikcoff,nel,kbasis,need)
      implicit real*8 (a-h,o-z)
#include "frag_cvb.fh"
      dimension ikcoff(0:nel,0:nel,0:nel)

      do 100 ifrag=1,nfrag
      do 101 ion=mnion_fr(ifrag),mxion_fr(ifrag)
      nelsing=nel_fr(ifrag)-2*ion
      if(nelsing.lt.0)goto 101
      do 200 iMs=1,nMs_fr(ifrag)
      nalfsing=nalf_fr(iMs,ifrag)-ion
      if(nalfsing.lt.0)goto 200
      do 300 iS=1,nS_fr(ifrag)
      if(i2s_fr(iS,ifrag).le.nelsing.and.
     >  i2s_fr(iS,ifrag).ge.2*nalfsing-nelsing)
     >  ikcoff(nelsing,nalfsing,i2S_fr(iS,ifrag))=1
300   continue
200   continue
101   continue
100   continue
      need=0
      do 400 nel1=0,nel
      do 401 nalf1=0,nel
      do 402 i2s1=0,nel
      if(ikcoff(nel1,nalf1,i2s1).eq.1)then
        ikcoff(nel1,nalf1,i2s1)=need
        nalf1_spin=(nel1+i2s1)/2
        need=need+ifns_cvb(nel1,nalf1_spin,kbasis)*ndet_cvb(nel1,nalf1)
      endif
402   continue
401   continue
400   continue
      return
      end
