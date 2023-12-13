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

function ICUNP(ICSPCK,L)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ICUNP
integer(kind=iwp), intent(in) :: ICSPCK(*), L
integer(kind=iwp) :: INTW, IPOW

INTW = ICSPCK((L+14)/15)
IPOW = 2**(28-2*mod(L-1,15))
ICUNP = mod(INTW/IPOW,4)

return

end function ICUNP
