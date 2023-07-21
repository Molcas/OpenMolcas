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

implicit real*8(A-H,O-Z)
real*8 Coeff(40)

Dum = 0.0D+00
do I=1,K
  Idx = 2*I-1
  Omega = Coeff(Idx)
  Alpha = Coeff(Idx+1)
  Dum = Dum+Omega*exp(-X*Alpha)
end do
ExpSum = Dum

return

end function ExpSum
