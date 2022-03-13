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

subroutine OffAtom(C1,C2,C3,C4,C5)

use Constants, only: Angstrom, deg2rad
use Definitions, only: wp

implicit none
real(kind=wp) :: C1(3), C2(3), C3(3), C4(3), C5(3)
real(kind=wp) :: D(3), E(3), LD, LE, U(3), V(3)
real(kind=wp), parameter :: R = 0.2767_wp/Angstrom, Theta = 36.72_wp*deg2rad

D(1) = (C2(1)+C3(1))/2-C1(1)
D(2) = (C2(2)+C3(2))/2-C1(2)
D(3) = (C2(3)+C3(3))/2-C1(3)
LD = sqrt(D(1)**2+D(2)**2+D(3)**2)
D(1) = D(1)/LD
D(2) = D(2)/LD
D(3) = D(3)/LD

U(1) = C2(1)-C1(1)
U(2) = C2(2)-C1(2)
U(3) = C2(3)-C1(3)
V(1) = C3(1)-C1(1)
V(2) = C3(2)-C1(2)
V(3) = C3(3)-C1(3)
E(1) = U(2)*V(3)-V(2)*U(3)
E(2) = U(3)*V(1)-V(3)*U(1)
E(3) = U(1)*V(2)-V(1)*U(2)
LE = sqrt(E(1)**2+E(2)**2+E(3)**2)
E(1) = E(1)/LE
E(2) = E(2)/LE
E(3) = E(3)/LE

C4(1) = C1(1)+R*(cos(Theta)*D(1)+sin(Theta)*E(1))
C4(2) = C1(2)+R*(cos(Theta)*D(2)+sin(Theta)*E(2))
C4(3) = C1(3)+R*(cos(Theta)*D(3)+sin(Theta)*E(3))
C5(1) = C1(1)+R*(cos(Theta)*D(1)-sin(Theta)*E(1))
C5(2) = C1(2)+R*(cos(Theta)*D(2)-sin(Theta)*E(2))
C5(3) = C1(3)+R*(cos(Theta)*D(3)-sin(Theta)*E(3))

return

end subroutine OffAtom
