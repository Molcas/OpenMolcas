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

subroutine CHO_HEAD(STRING,LINE,LENMAX,LUNIT)
!
! Purpose: print a header.

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: STRING
character, intent(in) :: LINE
integer(kind=iwp), intent(in) :: LENMAX, LUNIT
integer(kind=iwp) :: I, LENSTR, LENTOT

LENSTR = len(STRING)
LENTOT = min(LENSTR,LENMAX-2)
if (LENTOT > 0) then
  write(LUNIT,'(//,2X,A)') STRING(1:LENTOT)
  write(LUNIT,'(2X,80A)') (LINE,I=1,LENTOT)
else
  write(LUNIT,'(//,2X,A,/)') STRING(1:LENSTR)
end if

end subroutine CHO_HEAD
