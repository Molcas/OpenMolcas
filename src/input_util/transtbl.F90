!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001-2005, Valera Veryazov                             *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************
      Subroutine TransTbl(Filename)
!
!  translate basis file names
!
      character *256 Filename, DirName, OrigName
      character *256 Line
      Integer Strnln,irecl
      External StrnLn
      logical is_error,found
      iunit=20
      iunit=isfreeunit(iunit)
      ileft=0
      i=Strnln(FileName)
      found=.false.
      do while((.not.found).and.(i.gt.1))
        if(FileName(i:i).eq.'/') then
           ileft=i
           found=.true.
        endif
        i=i-1
      enddo
!
      if(.not.found) then
        ileft=0
        i=Strnln(FileName)
        found=.false.
        do while((.not.found).and.(i.gt.1))
          if(FileName(i:i).eq.'_') then
             ileft=i
             found=.true.
          endif
          i=i-1
        enddo
      endif
      DirName=Filename(1:ileft)
      i=index(Filename,' ')
      if(i.le.0) i=len(Filename)+1
      OrigName=Filename(ileft+1:i-1)
      LenOrig=i-ileft
      call molcas_open_ext2(iunit,DirName(1:ileft)//'trans.tbl',        &
     & 'sequential','formatted',iostat,.false.,irecl,'unknown',         &
     & is_error)
!      open(iunit,file=DirName(1:ileft)//'trans.tbl',
!     &        Form='FORMATTED',IOSTAT=IOStat)
      If (IOStat.ne.0) Then
         close(iunit)
         call molcas_open_ext2(iunit,'BASLIB_trans.tbl',                &
     & 'sequential','formatted',iostat,.false.,irecl,                   &
     &  'unknown',is_error)
!         open(iunit,file='BASLIB_trans.tbl',
!     &        Form='FORMATTED',IOSTAT=IOStat)
        If (IOStat.ne.0) Then
          write(6,*) 'trans.tbl is not found'
         close(iunit)
       return
      endif
      endif
20    read(iunit,'(a)',end=30, err=30) Line
      i=index(Line,OrigName(1:LenOrig-1))
      if(i.ne.1.or.Line(LenOrig:LenOrig).ne.' ') goto 20
       i=LenOrig+1
25     if (Line(i:i).eq.' ') then
        i=i+1
        if(i.lt.len(Line))goto 25
       endif
       ib=i
       ia=index(Line(ib:),' ')
       if(ia.eq.0) ia=len(Line)+1
       FileName=DirName(1:ileft)//Line(ib:ib+ia-1)
#ifdef _DEBUGPRINT_
       write(6,*) '*** Basis set was redirected to ',FileName
#endif
30     close(iunit)
       return
      end
