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

function GetDr2(K_Lap,X,Coeff)
!-----------------------------------------------------------------------
! Function : Get 2nd derivative of error function
!
! K_Lap
!  SUM Omega(I)*Alpha(I)*Alpha(I)*Exp(-Alpha(I)*X)
!  I=1
!-----------------------------------------------------------------------

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GetDr2
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: X, Coeff(2,20)
integer(kind=iwp) :: I
real(kind=wp) :: Alpha, Dum, Omega

Dum = Zero
do I=1,K_Lap
  Omega = Coeff(1,I)
  Alpha = Coeff(2,I)
  Dum = Dum+Omega*Alpha*Alpha*exp(-Alpha*X)
end do

GetDr2 = Two/(X*X*X)-Dum

return

end function GetDr2
