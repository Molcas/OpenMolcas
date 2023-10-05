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

subroutine Cho_1VecRd_SP(Vec,lVec,jVec,iSym,LstSP,nSP,iRedC,iLoc)
!
! Purpose: read vector jVec, sym. iSym, from disk. Read only
!          components from shell pairs in list LstSP. Use scratch
!          location iLoc to set index arrays. On input, iRedC
!          specifies the reduced set available at location iLoc:
!          specify -1 if not set (or unknown). On exit, iRedC
!          identifies the reduced set for which indices are
!          available at location iLoc. NOTE: only WA files!!

use Cholesky, only: Cho_AdrVec, iiBstRSh, InfVec, LuCho, LuPri, nnBstRSh, NumCho
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lVec, jVec, iSym, nSP, LstSP(nSP), iLoc
real(kind=wp), intent(out) :: Vec(lVec)
integer(kind=iwp), intent(inout) :: iRedC
integer(kind=iwp) :: iAdr, iAdr0, iOpt, irc, iRed, iShlAB, iSP, kV, lTot
character(len=*), parameter :: SecNam = 'Cho_1VecRd_SP'
integer(kind=iwp), external :: Cho_P_LocalSP

! Return if no vectors are available on disk.
! -------------------------------------------

if (NumCho(iSym) < 1) return

! Check that vector storage mode is word-addressable (WA).
! --------------------------------------------------------

if (Cho_AdrVec /= 1) then
  write(Lupri,*) SecNam,': WA address mode is required!'
  write(Lupri,*) 'Cho_AdrVec is: ',Cho_AdrVec,' (should be 1)'
  call Cho_Quit('WA address mode is required in '//SecNam,104)
end if

! Get reduced set of this vector.
! -------------------------------

if ((jVec > 0) .and. (jVec <= NumCho(iSym))) then
  iRed = InfVec(jVec,2,iSym)
else
  call Cho_Quit('Red. set error in '//SecNam,104)
  iRed = -999999
end if

! Set reduced set (if needed).
! ----------------------------

if (iRedC /= iRed) then
  call Cho_X_SetRed(irc,iLoc,iRed)
  if (irc /= 0) then
    write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
    call Cho_Quit('Error in '//SecNam,104)
  end if
  iRedC = iRed
end if

! Read vector elements.
! ---------------------

iAdr0 = InfVec(jVec,3,iSym)
kV = 1
do iSP=1,nSP
  iShlAB = Cho_P_LocalSP(LstSP(iSP))
  iOpt = 2
  lTot = nnBstRSh(iSym,iShlAB,iLoc)
  iAdr = iAdr0+iiBstRSh(iSym,iShlAB,iLoc)
  call dDAFile(LuCho(iSym),iOpt,Vec(kV),lTot,iAdr)
  kV = kV+lTot
end do

end subroutine Cho_1VecRd_SP
