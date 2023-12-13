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

function THETA(M,N)
! INTEGRATION OVER THETA. INCLUDES A FACTOR SIN(TH)
! FOR THE VOLUME ELEMENT

use crelop, only: GA
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: THETA
integer(kind=iwp), intent(in) :: M, N

if (mod(N,2) == 1) then
  THETA = Zero
else
  THETA = GA(M+2)*GA(N+1)/GA(M+N+3)
end if

return

end function THETA
