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

function Multipole_Expansion(q_A,q_B,Dipole_A,Dipole_B,R_A,R_B,R,Absolute)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Multipole_Expansion
real(kind=wp), intent(in) :: q_A, q_B, Dipole_A, Dipole_B, R_A, R_B, R
logical(kind=iwp) :: Absolute
real(kind=wp) :: dR_A, dR_B, E_A, E_B

dR_A = R_A-R
dR_B = R_B-R
call Multipole_E(q_A,Dipole_A,dR_A,E_A)
call Multipole_E(q_B,Dipole_B,dR_B,E_B)

if (Absolute) then
  Multipole_Expansion = abs(E_A+E_B)
else
  Multipole_Expansion = E_A+E_B
end if

return

end function Multipole_Expansion
