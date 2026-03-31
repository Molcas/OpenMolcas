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

subroutine MULMAT(NSS,XMATR,XMATI,ee,Z)

use Constants, only: Zero, cZero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSS
real(kind=wp), intent(in) :: XMATR(NSS,NSS), XMATI(NSS,NSS)
real(kind=wp), intent(out) :: ee
complex(kind=wp), intent(out) :: Z(NSS,NSS)
integer(kind=iwp) :: ISS, JSS

ee = Zero
Z(:,:) = cZero

do ISS=1,NSS
  do JSS=1,NSS
    ee = ee+XMATR(ISS,JSS)**2+XMATI(ISS,JSS)**2
    Z(ISS,JSS) = Z(ISS,JSS)+cmplx(XMATR(ISS,JSS),XMATI(ISS,JSS),kind=wp)
  end do
end do

end subroutine MULMAT
