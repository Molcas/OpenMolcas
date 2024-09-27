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

subroutine ChCar(iChCar,iGen,nGen)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iChCar(3)
integer(kind=iwp), intent(in) :: nGen, iGen(nGen)
integer(kind=iwp) :: i, iCar

! Generate characteristics for x, y, and z.

do iCar=1,3
  iChCar(iCar) = 0
  do i=1,nGen
    if (btest(iGen(i),iCar-1)) then
      iChCar(iCar) = 2**(iCar-1)
      exit
    end if
  end do
end do

return

end subroutine ChCar
