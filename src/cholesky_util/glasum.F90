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

implicit real*8(A-H,O-Z)
parameter(ZERO=0.0D+00,ONE=1.0D+00)
real*8 W(20), X(20)

Dum = ZERO
do I=1,K
  Val = exp((-Delta+ONE)*X(I))
  Dum = Dum+W(I)*Val
end do
GLaSum = Dum

return

end function GLaSum
