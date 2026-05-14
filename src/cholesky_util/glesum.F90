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

function GLeSum(K,X,W,Delta)
!======================================================================*
!----------------------------------------------------------------------*
!                                                                      *
!     FUNCTION : Compute the Function in Gauss-Legendre method         *
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
real(kind=wp) :: GLeSum
integer(kind=iwp), intent(in) :: K
real(kind=wp), intent(in) :: X(20), W(20), Delta
integer(kind=iwp) :: I
real(kind=wp) :: A, Derv, Dum, Val, Varl

Dum = Zero
do I=1,K
  A = One/(One-X(I))
  Varl = X(I)*A
  Derv = Varl*A
  Val = Derv*exp(-Delta*Varl)
  Dum = Dum+W(I)*Val
end do
GLeSum = Dum

return

end function GLeSum
