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

! Routine to generate the transformation for the second moments.
subroutine M2Trans(Rotte,TD)

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: Rotte(3,3)
real(kind=wp), intent(out) :: TD(6,6)

! The transformation of x2.

TD(1,1) = Rotte(1,1)*Rotte(1,1)
TD(2,1) = Rotte(1,1)*Rotte(2,1)
TD(3,1) = Rotte(1,1)*Rotte(3,1)
TD(4,1) = Rotte(2,1)*Rotte(2,1)
TD(5,1) = Rotte(2,1)*Rotte(3,1)
TD(6,1) = Rotte(3,1)*Rotte(3,1)

! The transformation of xy.

TD(1,2) = Rotte(1,1)*Rotte(1,2)+Rotte(1,2)*Rotte(1,1)
TD(2,2) = Rotte(1,1)*Rotte(2,2)+Rotte(1,2)*Rotte(2,1)
TD(3,2) = Rotte(1,1)*Rotte(3,2)+Rotte(1,2)*Rotte(3,1)
TD(4,2) = Rotte(2,1)*Rotte(2,2)+Rotte(2,2)*Rotte(2,1)
TD(5,2) = Rotte(2,1)*Rotte(3,2)+Rotte(2,2)*Rotte(3,1)
TD(6,2) = Rotte(3,1)*Rotte(3,2)+Rotte(3,2)*Rotte(3,1)

! The transformation of xz.

TD(1,3) = Rotte(1,1)*Rotte(1,3)+Rotte(1,3)*Rotte(1,1)
TD(2,3) = Rotte(1,1)*Rotte(2,3)+Rotte(1,3)*Rotte(2,1)
TD(3,3) = Rotte(1,1)*Rotte(3,3)+Rotte(1,3)*Rotte(3,1)
TD(4,3) = Rotte(2,1)*Rotte(2,3)+Rotte(2,3)*Rotte(2,1)
TD(5,3) = Rotte(2,1)*Rotte(3,3)+Rotte(2,3)*Rotte(3,1)
TD(6,3) = Rotte(3,1)*Rotte(3,3)+Rotte(3,3)*Rotte(3,1)

! The transformation of y2.

TD(1,4) = Rotte(1,2)*Rotte(1,2)
TD(2,4) = Rotte(1,2)*Rotte(2,2)
TD(3,4) = Rotte(1,2)*Rotte(3,2)
TD(4,4) = Rotte(2,2)*Rotte(2,2)
TD(5,4) = Rotte(2,2)*Rotte(3,2)
TD(6,4) = Rotte(3,2)*Rotte(3,2)

! The transformation of yz.

TD(1,5) = Rotte(1,2)*Rotte(1,3)+Rotte(1,3)*Rotte(1,2)
TD(2,5) = Rotte(1,2)*Rotte(2,3)+Rotte(1,3)*Rotte(2,2)
TD(3,5) = Rotte(1,2)*Rotte(3,3)+Rotte(1,3)*Rotte(3,2)
TD(4,5) = Rotte(2,2)*Rotte(2,3)+Rotte(2,3)*Rotte(2,2)
TD(5,5) = Rotte(2,2)*Rotte(3,3)+Rotte(2,3)*Rotte(3,2)
TD(6,5) = Rotte(3,2)*Rotte(3,3)+Rotte(3,3)*Rotte(3,2)

! The transformation of z2.

TD(1,6) = Rotte(1,3)*Rotte(1,3)
TD(2,6) = Rotte(1,3)*Rotte(2,3)
TD(3,6) = Rotte(1,3)*Rotte(3,3)
TD(4,6) = Rotte(2,3)*Rotte(2,3)
TD(5,6) = Rotte(2,3)*Rotte(3,3)
TD(6,6) = Rotte(3,3)*Rotte(3,3)

return

end subroutine M2Trans
