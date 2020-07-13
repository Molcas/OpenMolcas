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
      subroutine intchk_cvb(iarr,nmax,nread,ifc,a,lflag)
      implicit real*8 (a-h,o-z)
      parameter (nkeyw=3,ncmp=4)
      character*8 keyword(nkeyw)
      character*(*) a
#include "malloc_cvb.fh"
      dimension iarr(*),ito(1)
      save keyword
      data keyword/'NONE    ','ALL     ','TO      '/

      nread=0
      lf=lflag
1000  call fstring_cvb(keyword,nkeyw,istr,ncmp,1)
      if(istr.gt.0)lf=lflag
      if(istr.eq.1)then
c 'NONE'
        nread=0
      elseif(istr.eq.2)then
c 'ALL'
        if(lflag.eq.-1)then
          nread=nmax
          do 100 i=1,nmax
          iarr(i)=i
100       continue
        else
          nread=0
          lf=1-lflag
        endif
      elseif(istr.eq.3)then
c 'TO'
        if(nread.eq.nmax)then
          write(6,'(3a)')' Too many numbers specified in ',a,
     >      ' keyword!'
          call abend_cvb()
        elseif(nread.eq.0)then
          write(6,'(3a)')' No number before ',a,' -- TO keyword!'
          call abend_cvb()
        endif
        call int_cvb(ito,1,nr,ifc)
        if(nr.eq.-1)then
          write(6,'(3a)')' No number after ',a,' -- TO keyword!'
          call abend_cvb()
        endif
        ifrom=iarr(nread)
        if(ifrom.gt.ito(1))then
          write(6,*)' From greater than to:',ifrom,ito(1)
          call abend_cvb()
        elseif(nread+ito(1)-ifrom.gt.nmax)then
          write(6,'(3a)')' Too many numbers specified in ',a,
     >      ' keyword!'
          call abend_cvb()
        endif
        do 150 i=ifrom+1,ito(1)
        nread=nread+1
        iarr(nread)=i
150     continue
      else
        call int_cvb(iarr(1+nread),nmax-nread,nr,ifc)
        if(nread.gt.0)lf=lflag
        if(nr.eq.-1)then
          write(6,'(3a)')' Too many numbers specified in ',a,
     >      ' keyword!'
          call abend_cvb()
        endif
        nread=nread+nr
      endif
      if(istr.gt.0.or.nr.gt.0)goto 1000

      if(lflag.ne.-1)lflag=lf

      do 200 i=1,nread
      if(iarr(i).lt.1.or.iarr(i).gt.nmax)then
        write(6,'(3a,i5)')' Illegal ',a,' number read!',iarr(i)
        write(6,'(a,i3)')' Must be in the range 1 --',nmax
        call abend_cvb()
      endif
200   continue

c  Sort numbers and ensure uniqueness :
      call sorti_cvb(nread,iarr)
      ncnt=1
      do 300 i=2,nread
      if(iarr(i).ne.iarr(ncnt))then
        ncnt=ncnt+1
        iarr(ncnt)=iarr(i)
      endif
300   continue
      nread=min(ncnt,nread)

      return
      end
