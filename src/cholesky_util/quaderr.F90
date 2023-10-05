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

function QuadErr(K_Lap,X,Coeff)
!-----------------------------------------------------------------------
! Function : Calculate the error between Exponensial sum and 1/x
!-----------------------------------------------------------------------

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: QuadErr
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: X, Coeff(2,20)
real(kind=wp) :: S
real(kind=wp), external :: ExpSum

S = ExpSum(X,K_Lap,Coeff)

QuadErr = S-One/X

return

end function QuadErr
