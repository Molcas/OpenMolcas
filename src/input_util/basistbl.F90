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
      Subroutine BasisTbl(Label,BasDir)
!
!  a routine to translate basis set labels
!
      Character*(*) Label
      Character Temp*256, Line*256
      Character *(*) BasDir
      Integer Strnln,irecl
      External StrnLn
      Logical Exist,is_error
!     i=index(Label,'.')
!     if(i.gt.0) then
!     i=index(Label(i+1:),'.')
! return if we have a complete name
!     if(i.ne.0) return
!     endif
      Temp=BasDir//'/basis.tbl'
      Call f_Inquire(Temp,Exist)
      if(.not.Exist) Return
      jUnit=15
      iUnit=isfreeunit(jUnit)
      call molcas_Open_ext2(iunit,temp,'sequential','formatted',        &
     & iostat,.false.,irecl,'unknown',is_error)
!      open(unit=iUnit,File=Temp,form='FORMATTED',IOSTAT=IOStat)
      if(IOStat.ne.0) return
      iLast=StrnLn(Label)
!
!
!     Strip trailing dots
!
99    If (Label(iLast:iLast).eq.'.') Then
         iLast=iLast-1
         Go To 99
      End If
!     Write (*,*) 'Label(1:iLast)=',Label(1:iLast)
100     read(iUnit,'(a)',end=200,err=200) Line
        if(Line(1:1).eq.'#') goto 100
        if(Line.eq.' ') goto 100
        Call UpCase(Line)
!
!       Identify first non-blank segment in Line
!
        jlast=0
 98     If (Line(jlast+1:jlast+1).ne.' ') Then
           jLast = jLast +1
           Go To 98
        End If
        If (jlast.ne.ilast) Go To 100
!       i=index(Line,Label(1:iLast))
        i=index(Line(1:jlast),Label(1:iLast))
        if(i.ne.1) goto 100
        i=iLast+1
150      if(Line(i:i).eq.' ') then
         i=i+1
         goto 150
         endif
         ib=i
         ia=index(Line(ib:),' ')
         if(ia.eq.0) ia=len(Line)+1
#ifdef _DEBUGPRINT_
         write(6,'(3a)') Label(1:iLast),'translated to ',               &
     &      Line(ib:ib+ia-1)
#endif
        Label=Line(ib:ib+ia-1)
200     close(iUnit)
        return
        end
