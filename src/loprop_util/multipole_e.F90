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

subroutine Multipole_E(q,Dipole,R,E)

use Constants, only: Zero, One, Three
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: q, Dipole, R
real(kind=wp), intent(out) :: E
real(kind=wp) :: R_Inv, R2_Inv, xDipole, xMonopole

R_Inv = One/R
R2_Inv = R_Inv*R_Inv

if (R < Zero) then
  xMonopole = -q*R2_Inv
  xDipole = -Three*Dipole*R*R2_Inv*R2_Inv
else
  xMonopole = q*R2_Inv
  xDipole = Three*Dipole*R*R2_Inv*R2_Inv
end if

E = xMonopole+xDipole

return

end subroutine Multipole_E
