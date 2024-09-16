!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function get_ln_quit(lunit,icritical)
!***********************************************************************
! This function replaces function getln                                *
!                                                                      *
! It reads, broadcasts, and parses an input line                       *
!                                                                      *
! Blank lines or lines containing star (*) in column 1 are skipped     *
! Lines staring with exclaimation (!) in column 1 are skipped          *
!                                                                      *
! After this routine has been called, data can be retrieved using      *
! the subroutines                                                      *
!                                                                      *
!   Get_F(icol,array,n)  (for floating point values)                   *
!   Get_I(icol,iarry,n)  (for integer values)                          *
!   Get_S(icol,strgs,n)  (for character strings)                       *
!                                                                      *
! where icol is the first non-blank work to be taken, and n is the     *
! number of data.                                                      *
!                                                                      *
!***********************************************************************

use getline_mod, only: iEnd, iGetLine, iStrt, Line, MyUnit, nCol, Quit_On_Error
use Definitions, only: iwp, u6

implicit none
character(len=180) :: get_ln_quit
integer(kind=iwp), intent(in) :: lUnit, iCritical
integer(kind=iwp) :: i, icom, istatus, j, l
character(len=256) :: filename

Quit_On_Error = .false.
myunit = lunit
do
  read(lunit,'(A)',iostat=istatus) line
  if (istatus /= 0) exit
  igetline = igetline+1
  if ((line /= ' ') .and. (line(1:1) /= '*') .and. (line(1:1) /= '!')) exit
end do
if (istatus >= 0) then
  if (istatus == 0) then
    l = len(line)
    do i=1,l
      if (ichar(line(i:i)) == 9) line(i:i) = ' '
      if (line(i:i) == ';') line(i:l) = ' '
    end do
    ncol = 0
    j = 1
    do
      icom = 0
      do i=j,l
        if (line(i:i) == ',') then
          icom = icom+1
          if (icom /= 1) exit
        else
          if (line(i:i) /= ' ') exit
        end if
      end do
      if (i > l) then
        get_ln_quit = line
        return
      end if
      do j=i,l
        if ((line(j:j) == ' ') .or. (line(j:j) == ',')) exit
      end do
      ncol = ncol+1
      istrt(ncol) = i
      iend(ncol) = j-1
    end do
  end if

  filename = ' '
  inquire(unit=lunit,name=filename)
  if (filename /= ' ') then
    write(u6,'(a,a)') 'Error reading file=',filename
  else
    write(u6,'(a,i8)') 'Error reading unit=',lunit
  end if
  write(u6,'(a)') 'Line: ',line(1:80)
  Quit_On_Error = .true.
end if
if (icritical == 0) then
  Quit_On_Error = .true.
  return
end if

filename = ' '
inquire(unit=lunit,name=filename)
if (filename /= ' ') then
  write(u6,'(a,a)') 'EOF reached for file=',filename
else
  write(u6,'(a,i8)') 'EOF reached for unit=',lunit
end if

Quit_On_Error = .true.

end function get_ln_quit
