!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine AltErr(K_Lap,Coeff,T,VV,Eps)
!-----------------------------------------------------------------------
! Function : Calculate errors in each points of alternant
!-----------------------------------------------------------------------

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: Coeff(40), T(40)
real(kind=wp), intent(out) :: VV(40), Eps
integer(kind=iwp) :: I, I_Dim, J
real(kind=wp) :: X
real(kind=wp), external :: QuadErr

I_Dim = 2*K_Lap
Eps = Zero
J = I_Dim
do I=1,I_Dim
  X = T(J)
  VV(I) = QuadErr(K_Lap,X,Coeff)
  Eps = max(Eps,abs(VV(I)))
  J = J-1
end do

return

end subroutine AltErr
