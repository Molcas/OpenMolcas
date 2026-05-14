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

subroutine Cho_XCV_WrVec_Ser(irc,Vec,iSP)
!
! SERIAL VERSION
!     Simply write the partial vectors to disk at the appropriate
!     addresses on the vector files.

use Cholesky, only: iiBstRSh, LuCho, nnBstR, nnBstRSh, nSym, NumCho
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iSP
real(kind=wp), intent(_IN_) :: Vec(*)
integer(kind=iwp) :: iAdr, iAdr0, iSym, j, kV, lTot
integer(kind=iwp), parameter :: iOpt = 1

irc = 0

kV = 1
do iSym=1,nSym
  lTot = nnBstRSh(iSym,iSP,2)
  if (lTot > 0) then
    iAdr0 = iiBstRSh(iSym,iSP,2)
    do J=1,NumCho(iSym)
      iAdr = iAdr0+nnBstR(iSym,2)*(J-1)
      call DDAFile(LuCho(iSym),iOpt,Vec(kV),lTot,iAdr)
      kV = kV+lTot
    end do
  end if
end do

end subroutine Cho_XCV_WrVec_Ser
