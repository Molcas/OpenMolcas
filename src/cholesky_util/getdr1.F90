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

function GetDr1(K_Lap,X,Coeff)
!-----------------------------------------------------------------------
! Function : Get 1st derivative of error function
!
!                 K_Lap
! Dr1(X) = 1/X*X - SUM Omega(I)*Alpha(I)*Exp(-Alpha(I)*X)
!                  I=1
!-----------------------------------------------------------------------

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GetDr1
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: X, Coeff(40)
integer(kind=iwp) :: I, Idx
real(kind=wp) :: Alpha, Dum, Omega

Dum = Zero
do I=1,K_Lap
  Idx = 2*I-1
  Omega = Coeff(Idx)
  Alpha = Coeff(Idx+1)
  Dum = Dum+Omega*Alpha*exp(-Alpha*X)
end do

GetDr1 = Dum-One/(X*X)

return

end function GetDr1
