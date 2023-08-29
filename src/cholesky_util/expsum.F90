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

function ExpSum(X,K,Coeff)
!======================================================================*
!----------------------------------------------------------------------*
!                                                                      *
!     FUNCTION : Compute the Function in Exponential sum method.       *
!             K                                                        *
!     1/X =: SUM Omega(I)*EXP(-Alpha(I)*X)                             *
!            I=1                                                       *
!                                                                      *
!     Author(s): Akio Takatsuka (2007)                                 *
!                                                                      *
!----------------------------------------------------------------------*
!======================================================================*

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: ExpSum
real(kind=wp), intent(in) :: X, Coeff(2,20)
integer(kind=iwp), intent(in) :: K
integer(kind=iwp) :: I
real(kind=wp) :: Alpha, Dum, Omega

Dum = Zero
do I=1,K
  Omega = Coeff(1,I)
  Alpha = Coeff(2,I)
  Dum = Dum+Omega*exp(-X*Alpha)
end do
ExpSum = Dum

return

end function ExpSum
