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

subroutine OptRMS_Slapaf(x,y,nAt,RMS,RMSMax)

use Symmetry_Info, only: VarR, VarT
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: x(3,nAt)
real(kind=wp), intent(in) :: y(3,nAt)
real(kind=wp), intent(out) :: RMS, RMSMax
integer(kind=iwp) :: i, ixyz
real(kind=wp) :: diff, disp

! Only align if energy is rotational and translational invariant.
! (no weighting)

if (.not. (VarR .or. VarT)) then
  call Superpose(x,y,nAt,RMS,RMSMax)
else

  ! Otherwise, just compute RMS

  RMS = Zero
  RMSMax = Zero
  do i=1,nAt
    disp = Zero
    do ixyz=1,3
      diff = x(ixyz,i)-y(ixyz,i)
      RMS = RMS+diff**2
      disp = disp+diff**2
    end do
    if (sqrt(disp) > RMSMax) RMSMax = sqrt(disp)
  end do
  RMS = sqrt(RMS/real(nAt,kind=wp))
end if

return

end subroutine OptRMS_Slapaf
