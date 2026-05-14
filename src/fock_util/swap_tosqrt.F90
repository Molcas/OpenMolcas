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
! Copyright (C) Francesco Aquilante                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine swap_tosqrt(irc,iLoc,nRS,JSYM,XLT,Xab)

use Cholesky, only: iBas, iiBstR, iRS2F, nnBstR
use Symmetry_Info, only: Mul
use Data_Structures, only: NDSBA_Type
use Definitions, only: wp, iwp

#include "intent.fh"
implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iLoc, nRS, JSYM
type(NDSBA_Type), intent(_OUT_) :: XLT
real(kind=wp), intent(in) :: Xab(nRS)
integer(kind=iwp) :: iag, ias, ibg, ibs, iSyma, iSymb, jRab, kRab
integer(kind=iwp), external :: cho_isao

!                                                                      *
!***********************************************************************
!                                                                      *
if (JSYM /= 1) then ! NON TOTAL-SYMMETRIC

  do jRab=1,nnBstR(jSym,iLoc)

    kRab = iiBstr(jSym,iLoc)+jRab ! already in 1st red set

    iag = iRS2F(1,kRab) ! global address
    ibg = iRS2F(2,kRab)

    iSyma = cho_isao(iag) ! symmetry block
    iSymb = Mul(jSym,iSyma) ! sym(a) > sym(b)

    ias = iag-ibas(iSyma)
    ibs = ibg-ibas(iSymb)

    XLT%SB(iSyma,iSymb)%A2(ias,ibs) = sqrt(abs(Xab(kRab)))

  end do ! jRab loop

else if (JSYM == 1) then

  do jRab=1,nnBstR(jSym,iLoc)

    kRab = iiBstr(jSym,iLoc)+jRab ! already in 1st red set

    iag = iRS2F(1,kRab) ! global address
    ibg = iRS2F(2,kRab)

    iSyma = cho_isao(iag) ! sym(a)=sym(b)

    ias = iag-ibas(iSyma) ! address within that symm block
    ibs = ibg-ibas(iSyma)

    XLT%SB(iSyma,iSyma)%A2(ias,ibs) = sqrt(abs(Xab(kRab)))
    XLT%SB(iSyma,iSyma)%A2(ibs,ias) = sqrt(abs(Xab(kRab)))

  end do ! jRab loop

end if

irc = 0

return

end subroutine swap_tosqrt
