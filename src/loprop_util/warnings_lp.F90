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
integer(kind=iwp), intent(in) :: iCode, iLength
character(len=*), intent(out) :: Warning

! Binary codes:
! 1: Multiple minima found
! 2: Minima not within range
! 3: Ran out of iterations
! 4: No minima found

if (iLength < 25) then
  write(u6,*) 'Length of warning string must be at least 25 characters'
  call Abend()
end if

select case (iCode)
  case (1)
    Warning = 'Multiple minima found'
  case (2)
    Warning = 'Minima not within range'
  case (3)
    Warning = 'Ran out of iterations'
  case (4)
    Warning = 'No minima found'
  case default
    Warning = ''
end select

return

end subroutine Warnings_lp
