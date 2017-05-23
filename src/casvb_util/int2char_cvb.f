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
      subroutine int2char_cvb(a,int,iform)
c  Simulates internal write :    write(a,'(i<iform>)')int
      implicit real*8(a-h,o-z)
      character*(*) a
      character*1 blank,minus,cnumb(0:9)
      save blank,minus,cnumb
      data blank/' '/,minus/'-'/
      data cnumb/'0','1','2','3','4','5','6','7','8','9'/

      la=len(a)
      if(iform.gt.la)then
        write(6,*)' Illegal call to int2char_cvb:',iform,la
        call abend_cvb()
      endif
      dum=log10(DBLE(max(abs(int),1)))
      idum=nint(dum)
      if(abs(int).ge.10**idum)then
        iamax=idum+1
      else
        iamax=idum
      endif
      if(int.lt.0)iamax=iamax+1
      if(iamax.gt.iform)then
        write(6,*)' Integer too large in int2char_cvb:',int,iform
        call abend_cvb()
      endif
      ia=0
      do 100 i=1,iform-iamax
      ia=ia+1
100   a(ia:ia)=blank
      if(int.lt.0)then
        ia=ia+1
        a(ia:ia)=minus
        iamax=iamax-1
      endif
      int2=abs(int)
      do 200 i=iamax-1,0,-1
      ia=ia+1
      numb=int2/(10**i)
      a(ia:ia)=cnumb(numb)
200   int2=int2-numb*(10**i)
      if(int.eq.0)a(iform:iform)=cnumb(0)
      return
      end
