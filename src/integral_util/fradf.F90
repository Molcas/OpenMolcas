!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function fradf(x)
!***********************************************************************
!                                                                      *
! Object: to compute the angular contribution to the multipole integral*
!         between continuum basis functions within an R-matrix sphere  *
!         (theta integration)                                          *
!                                                                      *
!***********************************************************************

use rmat, only: expsum, l
use Definitions, only: wp

implicit none
real(kind=wp) :: fradf
real(kind=wp), intent(in) :: x

fradf = x**(l+2)*exp(-expsum*x*x)

end function fradf
