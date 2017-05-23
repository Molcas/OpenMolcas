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
      subroutine fout_cvb(f,a1,a2)
      implicit real*8 (a-h,o-z)
      character*(*) a1,a2
      character*15 b1
      character*46 b2
      character*12 b3
      save huge
      data huge/1d20/

      b1=a1
      b2=a2
      if(abs(f).ne.huge)then
        write(b3,'(d12.4)')f
      else
        b3='    Disabled'
      endif
      write(6,'(1x,3a)')b1,b2,b3
      return
      end
      subroutine fouti_cvb(fi,ni,a1,a2)
      implicit real*8 (a-h,o-z)
c      logical l
      character*(*) a1,a2
      character*15 b1
      character*46 b2
      character*12 b3
      dimension fi(ni)
      save huge
      data huge/1d20/
      b1=a1
      b2=a2
      b3='     ......'
      write(6,'(/,1x,3a)')b1,b2,b3
      b2=' '
c  Find IPOS : position of I index in string
      ichar0=ichar('0')
      ichar9=ichar('9')
      do ipos=15,1,-1
      if(ichar(b1(ipos:ipos)).ge.ichar0.and.
     >   ichar(b1(ipos:ipos)).le.ichar9)goto 100
      enddo
      write(6,*)' Fatal error in FOUTI!'
      call abend_cvb()
100   continue
      do i=1,ni
      if(abs(fi(i)).ne.huge)then
        write(b1(ipos:ipos),'(i1)')i
        write(b3,'(d12.4)')fi(i)
        write(6,'(1x,3a)')b1,b2,b3
      endif
      enddo
      return
      end
      subroutine foutij_cvb(fij,ni,nj,a1,a2)
      implicit real*8 (a-h,o-z)
c      logical l
      character*(*) a1,a2
      character*15 b1
      character*46 b2
      character*12 b3
      dimension fij(ni,nj)
      save huge
      data huge/1d20/
      b1=a1
      b2=a2
      b3='     ......'
      write(6,'(/,1x,3a)')b1,b2,b3
      b2=' '
c  Find IPOS/JPOS : position of I/J indices in string
      ichar0=ichar('0')
      ichar9=ichar('9')
      do jpos=15,1,-1
      if(ichar(b1(jpos:jpos)).ge.ichar0.and.
     >   ichar(b1(jpos:jpos)).le.ichar9)goto 200
      enddo
      write(6,*)' Fatal error in FOUTIJ!'
      call abend_cvb()
200   continue
      do ipos=jpos-1,1,-1
      if(ichar(b1(ipos:ipos)).ge.ichar0.and.
     >   ichar(b1(ipos:ipos)).le.ichar9)goto 300
      enddo
      write(6,*)' Fatal error in FOUTIJ!'
      call abend_cvb()
300   continue
      do j=1,nj
      do i=1,ni
      if(abs(fij(i,j)).ne.huge)then
        write(b1(ipos:ipos),'(i1)')i
        write(b1(jpos:jpos),'(i1)')j
        write(b3,'(d12.4)')fij(i,j)
        write(6,'(1x,3a)')b1,b2,b3
      endif
      enddo
      enddo
      return
      end
      subroutine iout_cvb(ii,a1,a2)
      implicit real*8 (a-h,o-z)
c      logical l
      character*(*) a1,a2
      character*15 b1
      character*46 b2
      character*12 b3
c      save huge
c      data huge/1d20/
      b1=a1
      b2=a2
      write(b3,'(i12)')ii
      write(6,'(1x,3a)')b1,b2,b3
      return
      end
      subroutine lout_cvb(l,a1,a2)
      implicit real*8 (a-h,o-z)
      logical l
      character*(*) a1,a2
      character*15 b1
      character*46 b2
      character*12 b3
c      save huge
c      data huge/1d20/
      b1=a1
      b2=a2
      if(l)then
        b3='        TRUE'
      else
        b3='       FALSE'
      endif
      write(6,'(1x,3a)')b1,b2,b3
      return
      end
