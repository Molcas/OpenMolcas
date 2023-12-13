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

function OVLMP(NP,ZP,NQ,ZQ)
!     <NP,ZP|NQ,ZQ>
! Note that both NP and NQ should be either even or odd, but not mixed.

use AMatrix, only: DFAC
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: OVLMP
integer(kind=iwp), intent(in) :: NP, NQ
real(kind=wp), intent(in) :: ZP, ZQ
integer(kind=iwp) :: IT11, IT22, IT33
real(kind=wp) :: RT1, ZAV

IT11 = 2*NP
IT22 = 2*NQ
RT1 = sqrt(ZP**(IT11+1)*ZQ**(IT22+1))/(DFAC(IT11)*DFAC(IT22))
RT1 = sqrt(RT1)
IT33 = NP+NQ
ZAV = Half*(ZP+ZQ)
OVLMP = RT1*DFAC(IT33)/sqrt(ZAV**(IT33+1))

return

end function OVLMP
