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

subroutine BasisTbl(Label,BasDir)

! a routine to translate basis set labels

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
character(len=*), intent(inout) :: Label
character(len=*), intent(in) :: BasDir
character(len=256) :: Temp, Line
integer(kind=iwp) :: i, ia, ib, iLast, irecl, istatus, iUnit, jlast, jUnit
integer(kind=iwp), external :: isFreeUnit, StrnLn
logical(kind=iwp) :: Exists, is_error

!i = index(Label,'.')
!if (i > 0) then
!  i=index(Label(i+1:),'.')
!  !return if we have a complete name
!  if (i /= 0) return
!endif
Temp = BasDir//'/basis.tbl'
call f_Inquire(Temp,Exists)
if (.not. Exists) return
jUnit = 15
iUnit = isfreeunit(jUnit)
call molcas_Open_ext2(iunit,temp,'sequential','formatted',istatus,.false.,irecl,'unknown',is_error)
!open(unit=iUnit,file=Temp,form='FORMATTED',iostat=istatus)
if (istatus /= 0) return
iLast = StrnLn(Label)

! Strip trailing dots

99 if (Label(iLast:iLast) == '.') then
  iLast = iLast-1
  Go To 99
end if
!write(u6,*) 'Label(1:iLast)=',Label(1:iLast)
100 read(iUnit,'(a)',end=200,err=200) Line
if (Line(1:1) == '#') goto 100
if (Line == ' ') goto 100
call UpCase(Line)

! Identify first non-blank segment in Line

jlast = 0
98 if (Line(jlast+1:jlast+1) /= ' ') then
  jLast = jLast+1
  Go To 98
end if
if (jlast /= ilast) Go To 100
!i = index(Line,Label(1:iLast))
i = index(Line(1:jlast),Label(1:iLast))
if (i /= 1) goto 100
i = iLast+1
150 if (Line(i:i) == ' ') then
  i = i+1
  goto 150
end if
ib = i
ia = index(Line(ib:),' ')
if (ia == 0) ia = len(Line)+1
#ifdef _DEBUGPRINT_
write(u6,'(3a)') Label(1:iLast),'translated to ',Line(ib:ib+ia-1)
#endif
Label = Line(ib:ib+ia-1)
200 close(iUnit)

return

end subroutine BasisTbl
