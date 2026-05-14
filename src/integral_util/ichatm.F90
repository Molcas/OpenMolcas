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

function iChAtm(Coor)
! iChAtm is an integer function which will return an integer such
! that the three first bits will represent the characteristics of
! the Cartesian components. If the bit is set then the Cartesian
! component will change sign if the symmetry operator contains a
! part which operates on that particular Cartesian direction.

use Symmetry_Info, only: iChCar, iOper, nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iChAtm
real(kind=wp), intent(in) :: Coor(3)
integer(kind=iwp) :: i, iCar, j, nOper

if (nIrrep == 8) then
  nOper = 3
else if (nIrrep == 4) then
  nOper = 2
else if (nIrrep == 2) then
  nOper = 1
else
  nOper = 0
end if

! Default that none of the Cartesians will change sign.
iChAtm = 0

! Loop over the Cartesian components.

do iCar=1,3

  ! Test if component is not zero. If zero no operator will change
  ! the sign of the component.
  if (abs(Coor(iCar)) < 1.0e-12_wp) cycle

  ! Here if the Component is none zero.

  ! Loop over the group generators and check if there is an
  ! operator that will change the sign.

  ! The generators are stored in positions 1, 2, and 4.

  do i=1,nOper    ! skip the unit operator -- i=0
    j = 2**(i-1)

    ! Test if symoperation will permute component

    if (iand(iOper(j),iChCar(iCar)) /= 0) then
      iChAtm = iChAtm+2**(iCar-1)
      exit
    end if
  end do
end do

return

end function iChAtm
