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

function GLaSum(K,X,W,Delta)
!======================================================================*
!----------------------------------------------------------------------*
!                                                                      *
!     FUNCTION : Compute the Function in Gauss-Laguerre method         *
!             K                                                        *
!     1/X =: SUM W(I)*Func(X(I))                                       *
!            I=1                                                       *
!                                                                      *
!     Author(s): Akio Takatsuka (2007)                                 *
!                                                                      *
!----------------------------------------------------------------------*
!======================================================================*

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GLaSum
integer(kind=iwp), intent(in) :: K
real(kind=wp), intent(in) :: X(20), W(20), Delta
integer(kind=iwp) :: I
real(kind=wp) :: Dum, Val

Dum = Zero
do I=1,K
  Val = exp((-Delta+One)*X(I))
  Dum = Dum+W(I)*Val
end do
GLaSum = Dum

return

end function GLaSum
