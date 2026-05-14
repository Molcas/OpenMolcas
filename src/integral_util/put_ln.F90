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

subroutine Put_ln(In_Line)
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

use getline_mod, only: iEnd, iStrt, Line, nCol
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: In_line
integer(kind=iwp) :: i, icom, j, l

Line = In_Line
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
  if (i > l) exit

  do j=i,l
    if ((line(j:j) == ' ') .or. (line(j:j) == ',')) exit
  end do
  ncol = ncol+1
  istrt(ncol) = i
  iend(ncol) = j-1
end do

end subroutine Put_ln
