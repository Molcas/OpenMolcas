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

subroutine setsavvb_cvb(recn)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: recn
integer(kind=iwp) :: iadd
logical(kind=iwp), external :: tstfile_cvb ! ... Files/Hamiltonian available ...
real(kind=wp), parameter :: recdef = 3200.2_wp

if (recn /= Zero) return
do iadd=0,99
  if (.not. tstfile_cvb(recdef+real(iadd,kind=wp))) then
    recn = recdef+real(iadd,kind=wp)
    return
  end if
end do
recn = recdef+99.0_wp

return

end subroutine setsavvb_cvb
