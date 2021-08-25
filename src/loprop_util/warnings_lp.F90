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

subroutine Warnings_lp(iCode,Warning,iLength)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iCode, iLength
character(len=*) :: Warning
integer(kind=iwp) :: i

! Binary codes:
! 1: Multiple minima found
! 2: Minima not within range
! 3: Ran out of iterations

if (iLength < 25) then
  write(u6,*) 'Length of warning string must be at least 25 characters'
  call Abend()
end if

do i=1,iLength
  Warning(i:i) = ' '
end do

if (iCode == 1) then
  Warning = 'Multiple minima found'
else if (iCode == 2) then
  Warning = 'Minima not within range'
else if (iCode == 3) then
  Warning = 'Ran out of iterations'
else if (iCode == 4) then
  Warning = 'No minima found'
end if

return

end subroutine Warnings_lp
