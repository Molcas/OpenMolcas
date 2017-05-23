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
      integer function ioemrg_cvb(ia1,na1,ia2,na2)
      implicit real*8 (a-h,o-z)
      dimension ia1(na1),ia2(na2)

      n1=1
      n2=1
      ioe=0
100   continue
      if(n1.gt.na1)then
        ioemrg_cvb=1-2*mod(ioe,2)
        return
      elseif(n2.gt.na2)then
        ioe=ioe+(na1-n1+1)*na2
        ioemrg_cvb=1-2*mod(ioe,2)
        return
      endif
      if(ia1(n1).lt.ia2(n2))then
        ioe=ioe+n2-1
        n1=n1+1
      elseif(ia1(n1).gt.ia2(n2))then
        n2=n2+1
      elseif(ia1(n1).eq.ia2(n2))then
        ioemrg_cvb=0
        return
      endif
      goto 100
      end
      integer function ioemrg2_cvb(ia1,na1,ia2,na2,ia12)
      implicit real*8 (a-h,o-z)
      dimension ia1(na1),ia2(na2),ia12(*)

      n1=1
      n2=1
      ioe=0
      n12=1
100   continue
      if(n1.gt.na1)then
        ioemrg2_cvb=1-2*mod(ioe,2)
        do i=n2,na2
        ia12(n12)=ia2(i)
        n12=n12+1
        enddo
        return
      elseif(n2.gt.na2)then
        ioe=ioe+(na1-n1+1)*na2
        ioemrg2_cvb=1-2*mod(ioe,2)
        do i=n1,na1
        ia12(n12)=ia1(i)
        n12=n12+1
        enddo
        return
      endif
      if(ia1(n1).lt.ia2(n2))then
        ioe=ioe+n2-1
        ia12(n12)=ia1(n1)
        n1=n1+1
        n12=n12+1
      elseif(ia1(n1).gt.ia2(n2))then
        ia12(n12)=ia2(n1)
        n2=n2+1
        n12=n12+1
      elseif(ia1(n1).eq.ia2(n2))then
        ioemrg2_cvb=0
        return
      endif
      goto 100
      end
