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

use GLBBAS, only: PGINT1A
use lucia_data, only: IBSO, ISMFSO, NTOOB, NTOOBS
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: T(*), H(*), TKK
integer(kind=iwp) :: K
integer(kind=iwp) :: IOFF, KOFF, KREL, KSM, NK
real(kind=wp) :: FAC
integer(kind=iwp), external :: IFRMR

KSM = ISMFSO(K)
KOFF = IBSO(KSM)
KREL = K-KOFF+1
NK = NTOOBS(KSM)

call SETVEC(H,Zero,NTOOB**2)

IOFF = IFRMR(PGINT1A(1)%A,1,KSM)
call COPVEC(T(IOFF+(KREL-1)*NK),H(IOFF+(KREL-1)*NK),NK)
TKK = H(IOFF-1+(KREL-1)*NK+KREL)
if (TKK /= Zero) then
  FAC = One/TKK
  call SCALVE(H(IOFF+(KREL-1)*NK),FAC,NK)
  !H(IOFF-1+(K-1)*NK+K) = H(IOFF-1+(K-1)*NK+K)-One
  H(IOFF-1+(KREL-1)*NK+KREL) = Zero
else
  !TKK = One
  TKK = Zero
end if

end subroutine T_ROW_TO_H
