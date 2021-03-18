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

subroutine TransTbl(Filename)

! translate basis file names

character*256 Filename, DirName, OrigName
character*256 Line
integer Strnln, irecl
external StrnLn
logical is_error, found
iunit = 20
iunit = isfreeunit(iunit)
ileft = 0
i = Strnln(FileName)
found = .false.
do while ((.not. found) .and. (i > 1))
  if (FileName(i:i) == '/') then
    ileft = i
    found = .true.
  end if
  i = i-1
end do

if (.not. found) then
  ileft = 0
  i = Strnln(FileName)
  found = .false.
  do while ((.not. found) .and. (i > 1))
    if (FileName(i:i) == '_') then
      ileft = i
      found = .true.
    end if
    i = i-1
  end do
end if
DirName = Filename(1:ileft)
i = index(Filename,' ')
if (i <= 0) i = len(Filename)+1
OrigName = Filename(ileft+1:i-1)
LenOrig = i-ileft
call molcas_open_ext2(iunit,DirName(1:ileft)//'trans.tbl','sequential','formatted',IOStat,.false.,irecl,'unknown',is_error)
!open(iunit,file=DirName(1:ileft)//'trans.tbl',form='FORMATTED',iostat=IOStat)
if (IOStat /= 0) then
  close(iunit)
  call molcas_open_ext2(iunit,'BASLIB_trans.tbl','sequential','formatted',IOStat,.false.,irecl,'unknown',is_error)
  !open(iunit,file='BASLIB_trans.tbl',form='FORMATTED',iostat=IOStat)
  if (IOStat /= 0) then
    write(6,*) 'trans.tbl is not found'
    close(iunit)
    return
  end if
end if
20 read(iunit,'(a)',end=30,err=30) Line
i = index(Line,OrigName(1:LenOrig-1))
if (i /= 1 .or. Line(LenOrig:LenOrig) /= ' ') goto 20
i = LenOrig+1
25 if (Line(i:i) == ' ') then
  i = i+1
  if (i < len(Line)) goto 25
end if
ib = i
ia = index(Line(ib:),' ')
if (ia == 0) ia = len(Line)+1
FileName = DirName(1:ileft)//Line(ib:ib+ia-1)
#ifdef _DEBUGPRINT_
write(6,*) '*** Basis set was redirected to ',FileName
#endif
30 close(iunit)

return

end subroutine TransTbl
