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

subroutine CtlDns(iDCRR,iDCRS,iDCRT,jOp)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iDCRR, iDCRS, iDCRT
integer(kind=iwp), intent(out) :: jOp(6)
integer(kind=iwp) :: iR, iRT, iRTS, iS, iT, iTS
integer(kind=iwp), external :: NrOpr

! Djl. Some care has to be taken here. Assume that there
! are two operators, T and S which generates the center
! pairs A,T(B) and A,S(B). If these pairs are symmetry
! related we will only

! Dij
iR = iDCRR
jOp(1) = NrOpr(iR)+1
! Dkl
iS = iDCRS
jOp(2) = NrOpr(iS)+1
! Dik
iT = iDCRT
jOp(3) = NrOpr(iT)+1
! Dil
iTS = ieor(iT,iS)
jOp(4) = NrOpr(iTS)+1
! Djk
iRT = ieor(iR,iT)
jOp(5) = NrOpr(iRT)+1
! Djl
iRTS = ieor(iRT,iS)
jOp(6) = NrOpr(iRTS)+1

return

end subroutine CtlDns
