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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Pot_nuc(CCoor,pot,nGrid)

use Basis_Info
use Center_Info

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 CCoor(3,nGrid), pot(nGrid)
real*8 C(3), TC(3)
integer iStabM(0:7), iDCRT(0:7)

! compute nuclear contribution to potential

do iGrid=1,nGrid
  pot(iGrid) = 0d0
end do

!hjw is this always correct?
istabm(0) = 0
nstabm = 1

kdc = 0
if (nCnttp > 0) kdc = -dbsc(1)%nCntr ! to make sure we start at 0
do kCnttp=1,nCnttp
  kdc = kdc+dbsc(kCnttp)%nCntr
  if (dbsc(kCnttp)%Charge == Zero) cycle

  do kCnt=1,dbsc(kCnttp)%nCntr

    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      do iGrid=1,nGrid
        r12 = sqrt((TC(1)-CCoor(1,iGrid))**2+(TC(2)-CCoor(2,iGrid))**2+(TC(3)-CCoor(3,iGrid))**2)
        if (r12 > 1.d-8) pot(iGrid) = pot(iGrid)+dbsc(kCnttp)%Charge*fact/r12
      end do

    end do
  end do
end do

return

end subroutine Pot_nuc
