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

subroutine Cho_XCV_RdVec_(irc,Vec,myRankSP,n_myRankSP,NVT,J1,J2,iSym)
!
! Read the vector blocks.

use ChoSwp, only: nnBstRSh

implicit none
integer irc
real*8 Vec(*)
integer n_myRankSP
integer myRankSP(n_myRankSP)
integer NVT
integer J1, J2, iSym
#include "cholesky.fh"
integer iOpt
parameter(iOpt=2)
integer kV, n, i
integer lTot, iAdr, iAdr0
integer iSP

irc = 0

n = J2-J1+1
iAdr0 = 0
kV = 1
do i=1,n_myRankSP
  iSP = myRankSP(i)
  lTot = nnBstRSh(iSym,iSP,2)*n
  if (lTot > 0) then
    iAdr = iAdr0+nnBstRSh(iSym,iSP,2)*(J1-1)
    call DDAFile(LuTmp(iSym),iOpt,Vec(kV),lTot,iAdr)
    kV = kV+lTot
  end if
  iAdr0 = iAdr0+nnBstRSh(iSym,iSP,2)*NVT
end do

end subroutine Cho_XCV_RdVec_
