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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_WrVec_Par(irc,Vec,NVT,myRankSP,SP)
!
! PARALLEL VERSION
!     Write the vectors in blocks.

use Cholesky, only: LuTmp, nnBstRSh, nSym
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: NVT(*), myRankSP(*), SP
real(kind=wp), intent(_IN_) :: Vec(*)
integer(kind=iwp) :: iAdr, iSP, iSym, j, kV, lTot
integer(kind=iwp), parameter :: iOpt = 1

irc = 0

iSP = myRankSP(SP)
kV = 1
do iSym=1,nSym
  lTot = nnBstRSh(iSym,iSP,2)*NVT(iSym)
  if (lTot > 0) then
    iAdr = 0
    do j=1,SP-1
      iAdr = iAdr+nnBstRSh(iSym,myRankSP(j),2)*NVT(iSym)
    end do
    call DDAFile(LuTmp(iSym),iOpt,Vec(kV),lTot,iAdr)
    kV = kV+lTot
  end if
end do

end subroutine Cho_XCV_WrVec_Par
