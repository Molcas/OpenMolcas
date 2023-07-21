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

implicit real*8(A-H,O-Z)
real*8 VV(40), T(40), Coeff(40)
IDim = 2*K_Lap
Eps = 0.0D+00
J = IDim
do I=1,IDim
  X = T(J)
  VV(I) = QuadErr(K_Lap,X,Coeff)
  Eps = max(Eps,abs(VV(I)))
  J = J-1
end do

return

end subroutine AltErr
