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

subroutine Phi_point(iPhi,nPhi,Cos_Phi,Sin_Phi,w_Phi)

use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iPhi, nPhi
real(kind=wp), intent(out) :: Cos_Phi, Sin_Phi, w_Phi
real(kind=wp) :: q
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

q = Pi*(Two*real(iPhi,kind=wp)-One)/real(nPhi,kind=wp)
if (abs(cos(q)) > Thrs) then
  Cos_Phi = cos(q)
else
  Cos_Phi = Zero
end if
if (abs(sin(q)) > Thrs) then
  Sin_Phi = sin(q)
else
  Sin_Phi = Zero
end if
w_Phi = Two*Pi/real(nPhi,kind=wp)

return

end subroutine Phi_point
