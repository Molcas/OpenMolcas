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

subroutine Truncate_Grid(R,nR,nR_Eff,Radius_Max)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nR
real(kind=wp), intent(in) :: R(2,nR), Radius_Max
integer(kind=iwp), intent(inout) :: nR_Eff
integer(kind=iwp) :: i, nTmp

nTmp = nR_Eff
do i=1,nTmp
  if (R(1,i) > Radius_Max) then
    nR_Eff = i-1
    exit
  end if
end do

return

end subroutine Truncate_Grid
