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

subroutine CHO_CLOSE(LUNIT,STAT)
!
! Purpose: close sequential unformatted fortran file.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: LUNIT
character(len=*), intent(in) :: STAT

if ((LUNIT < 1) .or. (LUNIT > 99)) then
  call CHO_QUIT('CHO_CLOSE: unit out of bounds!',104)
else
  close(LUNIT,status=STAT)
  LUNIT = -1
end if

end subroutine CHO_CLOSE
