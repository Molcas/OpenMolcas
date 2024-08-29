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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine rKappa_Zeta(rKappa,Zeta,nZeta)
!***********************************************************************
!                                                                      *
! Object: modify rkappa                                                *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    :                                                         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!             September '00.                                           *
!***********************************************************************

use Constants, only: OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta
real(kind=wp), intent(inout) :: rKappa(nZeta)
real(kind=wp), intent(in) :: Zeta(nZeta)

rKappa(:) = rKappa(:)*Zeta(:)**(-OneHalf)

end subroutine rKappa_Zeta
