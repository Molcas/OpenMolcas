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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_AmpDiag(irc,Diag,EOcc,EVir)
!
! Jonas Bostrom, Jan. 2010.
!
! Purpose: Construct diagonal for decomposition of amplitude
!          vectors.

use Symmetry_Info, only: Mul
use Cholesky, only: nSym
use ChoMP2, only: iMoMo, iOcc, iVir, nMoMo, nOcc, nVir
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: iA, iAI, iI, iSym, iSymA, iSymI, iVecType, kD0, kD1, kD2
real(kind=wp) :: DE, Ei

irc = 0

! Initialization.
! ---------------

iVecType = 6
kD0 = 0

! Construct Diagonal.
! -------------------

do iSym=1,nSym
  do iSymI=1,nSym
    iSymA = Mul(iSymI,iSym)
    kD1 = kD0+iMoMo(iSymA,iSymI,iVecType)
    do iI=1,nOcc(iSymI)
      kD2 = kD1+nVir(iSymA)*(iI-1)
      Ei = EOcc(iOcc(iSymI)+iI)
      do iA=1,nVir(iSymA)
        iAI = kD2+iA
        DE = Two*(EVir(iVir(iSymA)+iA)-Ei)
        Diag(iAI) = Diag(iAI)/DE
      end do
    end do
  end do
  kD0 = kD0+nMoMo(iSym,iVecType)
end do

return

end subroutine ChoMP2g_AmpDiag
