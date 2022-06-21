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

use Definitions, only: iwp, u6

implicit none
character(len=256), intent(inout) :: Filename
character(len=256) :: DirName, OrigName, Line
integer(kind=iwp) :: i, ia, ib, ileft, irecl, istatus, iunit, LenOrig
integer(kind=iwp), external :: isFreeUnit, StrnLn
logical(kind=iwp) :: is_error, found

iunit = isFreeUnit(20)
ileft = 0
i = StrnLn(FileName)
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
  i = StrnLn(FileName)
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
call molcas_open_ext2(iunit,DirName(1:ileft)//'trans.tbl','sequential','formatted',istatus,.false.,irecl,'unknown',is_error)
!open(iunit,file=DirName(1:ileft)//'trans.tbl',form='formatted',iostat=istatus)
if (istatus /= 0) then
  close(iunit)
  call molcas_open_ext2(iunit,'BASLIB_trans.tbl','sequential','formatted',istatus,.false.,irecl,'unknown',is_error)
  !open(iunit,file='BASLIB_trans.tbl',form='formatted',iostat=istatus)
  if (istatus /= 0) then
    write(u6,*) 'trans.tbl is not found'
    close(iunit)
    return
  end if
end if
do
  read(iunit,'(a)',iostat=istatus) Line
  if (istatus /= 0) exit
  i = index(Line,OrigName(1:LenOrig-1))
  if ((i == 1) .and. (Line(LenOrig:LenOrig) == ' ')) exit
end do
if (istatus == 0) then
  do i=LenOrig+1,len(Line)
    if (Line(i:i) /= ' ') exit
  end do
  ib = i
  ia = index(Line(ib:),' ')
  if (ia == 0) ia = len(Line)+1
  FileName = DirName(1:ileft)//Line(ib:ib+ia-1)
# ifdef _DEBUGPRINT_
  write(u6,*) '*** Basis set was redirected to ',FileName
# endif
end if
close(iunit)

return

end subroutine TransTbl
