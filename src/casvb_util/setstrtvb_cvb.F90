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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine setstrtvb_cvb(recn)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: recn
integer(kind=iwp) :: iadd, jadd
real(kind=wp), parameter :: recdef = 3200.2_wp
logical(kind=iwp), external :: tstfile_cvb ! ... Files/Hamiltonian available ...

if (recn /= Zero) return
if (.not. tstfile_cvb(recdef)) return
do iadd=1,99
  if (.not. tstfile_cvb(recdef+real(iadd,kind=wp))) then
    jadd = iadd-1
    recn = recdef+real(jadd,kind=wp)
    return
  end if
end do

return

end subroutine setstrtvb_cvb
