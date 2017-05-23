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
      subroutine party2_cvb(iperm,n,party)
      implicit real*8(a-h,o-z)
      dimension iperm(n)

c  Following caters for non-contiguous integers :
      ntransp=0
100   continue
      do 200 i=1,n-1
      if(iperm(i).gt.iperm(i+1))then
        ntransp=ntransp+1
        iswp=iperm(i)
        iperm(i)=iperm(i+1)
        iperm(i+1)=iswp
        do 300 j=i-1,1,-1
        if(iperm(j).gt.iperm(j+1))then
          ntransp=ntransp+1
          iswp=iperm(j)
          iperm(j)=iperm(j+1)
          iperm(j+1)=iswp
        endif
300     continue
        goto 100
      endif
200   continue
      if(mod(ntransp-n,2).eq.0)then
        party=1d0
      else
        party=-1d0
      endif
      return
      end
