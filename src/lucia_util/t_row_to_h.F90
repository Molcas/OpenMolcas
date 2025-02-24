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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

subroutine T_ROW_TO_H(T,H,K,TKK)
! Set H integrals
!
!    Column K : H(P,K) = T(P,K)/T(K,K), P /= K
!    Other Columns     = 0
! - and return T_{kk} in TKK
!
! Jeppe Olsen, Jan 98
! For rotation of CI vectors

use lucia_data, only: IBSO, ISMFSO, NTOOB, NTOOBS, PGINT1A
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: T(NTOOB**2)
real(kind=wp), intent(out) :: H(NTOOB**2), TKK
integer(kind=iwp), intent(in) :: K
integer(kind=iwp) :: IOFF, KOFF, KREL, KSM, NK

KSM = ISMFSO(K)
KOFF = IBSO(KSM)
KREL = K-KOFF+1
NK = NTOOBS(KSM)

H(:) = Zero

IOFF = PGINT1A(1)%A(KSM)
H(IOFF+(KREL-1)*NK:IOFF+KREL*NK-1) = T(IOFF+(KREL-1)*NK:IOFF+KREL*NK-1)
TKK = H(IOFF-1+(KREL-1)*NK+KREL)
if (TKK /= Zero) then
  H(IOFF+(KREL-1)*NK:IOFF+KREL*NK-1) = H(IOFF+(KREL-1)*NK:IOFF+KREL*NK-1)/TKK
  !H(IOFF-1+(K-1)*NK+K) = H(IOFF-1+(K-1)*NK+K)-One
  H(IOFF-1+(KREL-1)*NK+KREL) = Zero
!else
!  TKK = One
end if

end subroutine T_ROW_TO_H
