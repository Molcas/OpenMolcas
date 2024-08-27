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

function iChxyz(Coord,iGen,nGen)

use Constants

implicit none
integer iChxyz
integer nGen
real*8 Coord(3)
integer iGen(nGen), iChCar(3)
integer iCar

call ChCar(iChCar,iGen,nGen)

iChxyz = 0
do iCar=1,3
  if (Coord(iCar) /= Zero) iChxyz = iChxyz+iChCar(iCar)
end do

return

end function iChxyz
