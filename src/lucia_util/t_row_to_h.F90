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

use Constants, only: Zero, One
use GLBBAS, only: PGINT1A
use lucia_data, only: NTOOB, IBSO, ISMFSO, NTOOBS

implicit none
integer K
real*8 TKK
! Input (in blocked form)
real*8 T(*)
! Output (also in blocked form)
real*8 H(*)

integer, external :: IFRMR
integer KSM, KOFF, KREL, NK, IOFF
real*8 FAC

KSM = ISMFSO(K)
KOFF = IBSO(KSM)
KREL = K-KOFF+1
NK = NTOOBS(KSM)

call SETVEC(H,ZERO,NTOOB**2)

IOFF = IFRMR(PGINT1A(1)%I,1,KSM)
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
