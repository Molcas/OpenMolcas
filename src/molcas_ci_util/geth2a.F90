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

function GETH2A(I,J,K,L,TUVX)

use Index_Functions, only: iTri
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GETH2A
integer(kind=iwp), intent(in) :: I, J, K, L
real(kind=wp), intent(in) :: TUVX(*)

GETH2A = TUVX(iTri(iTri(I,J),iTri(K,L)))

return

end function GETH2A
