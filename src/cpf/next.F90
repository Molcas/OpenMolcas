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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine NEXT(P,DPS,CN)

use cpf_global, only: IADDP, IPRINT, ITPUL, Lu_CI, NCONF
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: P(*), DPS(*)
real(kind=wp), intent(in) :: CN(*)
integer(kind=iwp) :: I, IAD, IIN, ITM, J
real(kind=wp) :: CTOT

IAD = IADDP(1)
call dDAFILE(Lu_CI,2,P,NCONF,IAD)
ITM = ITPUL-1
do I=1,ITM
  IIN = I+1
  CTOT = Zero
  do J=IIN,ITPUL
    CTOT = CTOT+CN(J)
  end do
  IAD = IADDP(I+1)
  call dDAFILE(Lu_CI,2,DPS,NCONF,IAD)
  P(1:NCONF) = P(1:NCONF)+CTOT*DPS(1:NCONF)
end do
if (IPRINT >= 15) write(u6,19) (P(I),I=1,NCONF)

return

19 format(6X,'C(NEXT)',5F10.6)

end subroutine NEXT
