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
! Copyright (C) 1986, Bernd Artur Hess                                 *
!***********************************************************************

subroutine RELOP()
! SUBROUTINE RELOP INITIALIZES THE COMMON BLOCK USED BY
! THE RELOP PACKAGE
! V 1.0 - 12.3.86 - BERND HESS

use crelop, only: GA, IMAX
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: N
real(kind=wp), external :: GAM

IMAX = size(GA)
do N=1,IMAX
  GA(N) = GAM(N-1)
end do

return

end subroutine RELOP
