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

implicit real*8(A-H,O-Z)
parameter(ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
real*8 Coeff(40)

Dum = ZERO
do I=1,K_Lap
  Idx = 2*I-1
  Omega = Coeff(Idx)
  Alpha = Coeff(Idx+1)
  Dum = Dum+Omega*Alpha*Alpha*exp(-Alpha*X)
end do

GetDr2 = TWO/(X*X*X)-Dum

return

end function GetDr2
