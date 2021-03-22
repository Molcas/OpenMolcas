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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_RdInp()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Parse and process RPA input.

implicit none
#include "WrkSpc.fh"
#include "rpa_config.fh"
#include "rpa_data.fh"

character*9 SecNam
parameter(SecNam='RPA_RdInp')
integer mLine
parameter(mLine=100000)

logical Reduce_Prt
external Reduce_Prt
character*180 Get_Ln
external Get_Ln
integer iPrintLevel, RPA_iUHF
external iPrintLevel, RPA_iUHF

character*4 Keyword
character*180 Line

logical doTitle
logical EndOfInput
logical DebugPrint

integer ip_Integer, l_Integer
integer Lu
integer iLine
integer iUHF
integer i

! set default print level
iPrint = max(iPrintLevel(-1),0)
if (iPrint > 0 .and. iPrint < 3) then
  if (Reduce_Prt()) then
    iPrint = 0
  end if
end if

! set debug print
#ifdef _DEBUGPRINT_
DebugPrint = iPrint > 4
#else
DebugPrint = .false.
#endif
if (DebugPrint) then
  write(6,'(A,A,I3)') SecNam,': default print level is',iPrint
end if

! Set restricted (1) or unrestricted (2)
iUHF = RPA_iUHF()

! set defaults
dRPA = .true.
SOSEX = .false.
LumOrb = .false.
nTitle = 0
nFreeze(1) = 0
nFreeze(2) = 0

! open input file and find RPA section
Lu = 17
call SpoolInp(Lu)
rewind(Lu)
call RdNLst(Lu,'RPA')

! allocate dummy arrays for reading input
! (dimension should be the longest needed for reading)
l_Integer = 2
call GetMem('InpRdI','Allo','Inte',ip_Integer,l_Integer)

! parse input
doTitle = .false.
EndOfInput = .false.
iLine = 0
do while (.not. EndOfInput .and. iLine < mLine)
  iLine = iLine+1 ! eliminate risk of an infinite loop
  Line = Get_Ln(Lu)
  call StdFmt(Line,Keyword)
  if (DebugPrint) then
    write(6,'(A,A)') SecNam,' processing line:'
    write(6,'(A)') Line
    write(6,'(A,A)') 'Key=',Keyword
  end if
  !****************************
  if (Keyword(1:1) == ' ') then
    !****************************
    ! blank line: pass
    continue
    !*********************************
  else if (Keyword(1:1) == '*') then
    !*********************************
    ! comment line: pass
    continue
    !***********************************
  else if (Keyword(1:3) == 'END') then
    !***********************************
    ! end of input
    EndOfInput = .true.
    !*******************************
  else if (Keyword == 'TITL') then
    !*******************************
    ! title
    doTitle = .true.
    nTitle = nTitle+1
    if (nTitle > mTitle) then
      write(6,'(A,I8)') 'Maximum number of title lines is',mTitle
      write(6,'(A)') 'Current input line:'
      write(6,'(A)') Line
      call RPA_Warn(2,'Too many title lines in RPA input')
    else
      Line = Get_Ln(Lu)
      Title(nTitle) = Line(1:80)
    end if
    !*******************************
  else if (Keyword == 'PRIN') then
    !*******************************
    ! print level
    doTitle = .false.
    call RPA_ReadIntegerInput('PRIN',1,Lu,iWork(ip_Integer),l_Integer)
    iPrint = max(iWork(ip_Integer),0)
#   ifdef _DEBUGPRINT_
    DebugPrint = DebugPrint .or. iPrint > 4
#   endif
    !*******************************
  else if (Keyword == 'LUMO') then
    !*******************************
    ! use orbitals from InpOrb file
    doTitle = .false.
    LumOrb = .true.
    !*******************************
  else if (Keyword == 'RUNO') then
    !*******************************
    ! use orbitals from Runfile
    doTitle = .false.
    LumOrb = .false.
    !*******************************
  else if (Keyword == 'DRPA') then
    !*******************************
    ! direct RPA
    doTitle = .false.
    dRPA = .true.
    SOSEX = .false.
    !*******************************
  else if (Keyword == 'SOSE') then
    !*******************************
    ! direct RPA + SOSEX
    doTitle = .false.
    dRPA = .true.
    SOSEX = .true.
    !*******************************
  else if (Keyword == 'ALLE') then
    !*******************************
    ! All electrons correlated (no frozen)
    doTitle = .false.
    call iZero(nFro,16)
    !*******************************
  else if (Keyword == 'FREE') then
    !*******************************
    ! Freeze occupied orbitals
    doTitle = .false.
    call iZero(nFro,16)
    call RPA_ReadIntegerInput('FREE',iUHF,Lu,iWork(ip_Integer),l_Integer)
    do i=1,iUHF
      nFreeze(i) = iWork(ip_Integer-1+i)
    end do
    !*******************************
  else if (Keyword == 'DELE') then
    !*******************************
    ! Delete virtual orbitals
    doTitle = .false.
    call RPA_Warn(2,'Virtual orbital deletion not implemented yet!')
  else
    ! no keyword match
    ! => either a title line or an unknown keyword
    if (doTitle) then
      ! process as title line
      nTitle = nTitle+1
      if (nTitle > mTitle) then
        write(6,'(A,I8)') 'Maximum number of title lines is',mTitle
        write(6,'(A)') 'Current input line:'
        write(6,'(A)') Line
        call RPA_Warn(2,'Too many title lines in RPA input')
      else
        Title(nTitle) = Line(1:80)
      end if
    else
      ! unknown keyword
      write(6,'(A)') 'Offending input line:'
      write(6,'(A)') Line
      write(6,'(A,A,A)') 'Equivalent keyword input "',Keyword,'" not recognized!'
      call RPA_Warn(2,'RPA input keyword not recognized')
    end if
  end if
end do

! deallocation
call GetMem('InpRdI','Free','Inte',ip_Integer,l_Integer)

! close input file
call Close_LuSpool(Lu)

end subroutine RPA_RdInp
