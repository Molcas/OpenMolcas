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

subroutine Cho_SetShP2RS_2(irc,iLoc,iShlAB,nAB)
!
! Purpose: set mapping from shell pair iShlAB to reduced shell pair
!          Based on reduced set stored at location iLoc = 2 or 3.
!          If a non-zero code (irc) is returned, nothing has been
!          set!!

use Index_Functions, only: nTri_Elem
use Cholesky, only: iiBstR, iiBstRSh, IndRed, iOff_Batch, iShP2RS, iSP2F, MySP, nBstSh, nnBstRSh, nSym
#ifdef _DEBUGPRINT_
use Cholesky, only: IndRSh
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iLoc, iShlAB, nAB(*)
integer(kind=iwp) :: iAB, iAB1, iAB2, iShlA, iShlB, iSym, jAB, kAB, l_iShP2RS, lTst, NumAB
#ifdef _DEBUGPRINT_
character(len=*), parameter :: SecNam = 'Cho_SetShP2RS_2'
#endif

! Check allocations.
! ------------------

call Cho_InvPck(iSP2F(mySP(iShlAB)),iShlA,iShlB,.true.)
if (iShlA == iShlB) then
  NumAB = nTri_Elem(nBstSh(iShlA))
else
  NumAB = nBstSh(iShlA)*nBstSh(iShlB)
end if
lTst = 2*NumAB
l_iShP2RS = 0
if (allocated(iShP2RS)) l_iShP2RS = size(iShP2RS)
if ((l_iShP2RS < 1) .or. (l_iShP2RS < lTst)) then
  irc = 102
  return
end if

! Check iLoc.
! -----------

if ((iLoc < 2) .or. (iLoc > 3)) then
  irc = 104
  return
end if

! Set mapping array.
! iShP2RS(1,AB) = index in current reduced set, symmetry reduced.
! iShP2RS(2,AB) = symmetry block.
! Zeros are returned if the element AB is not a member of the
! current reduced set.
! ---------------------------------------------------------------

iShP2RS(:,1:NumAB) = 0

do iSym=1,nSym
  if (nAB(iSym) > 0) then
    iAB1 = iiBstRSh(iSym,iShlAB,iLoc)+1
    iAB2 = iAB1+nnBstRSh(iSym,iShlAB,iLoc)-1
    do iAB=iAB1,iAB2
      jAB = IndRed(iiBstR(iSym,iLoc)+iAB,iLoc) ! addr in 1st rs
      kAB = IndRed(jAB,1) ! addr in full shell pair
#     ifdef _DEBUGPRINT_
      nErr = 0
      if (IndRSh(jAB) /= iSP2F(mySP(iShlAB))) then
        write(Lupri,*) SecNam,': inconsistent shell pairs!'
        write(Lupri,*) SecNam,': from input: ',iSP2F(mySP(iShlAB)),'  from IndRsh: ',IndRSh(jAB)
        nErr = nErr+1
      end if
      if ((kAB < 1) .or. (kAB > NumAB)) then
        write(Lupri,*) SecNam,': shell pair address error!'
        write(Lupri,*) SecNam,': kAB = ',kAB
        write(Lupri,*) SecNam,': min and max allowed: 1 ',NumAB
        nErr = nErr+1
      end if
      if (nErr /= 0) then
        write(Lupri,*) SecNam,': Shell A, B, AB: ',iShlA,iShlB,iShlAB
        write(Lupri,*) SecNam,': iLoc: ',iLoc
        write(Lupri,*) SecNam,': symmetry block: ',iSym
        write(Lupri,*) SecNam,': red. set address, first and current: ',jAB,iiBstR(iSym,iLoc)+iAB
        call Cho_Quit('Error detected in '//SecNam,104)
      end if
#     endif
      iShP2RS(1,kAB) = iOff_Batch(iSym,iShlAB)+iAB-iAB1+1
      iShP2RS(2,kAB) = iSym
    end do
  end if
end do

! Set return code 0: all ok!
! --------------------------

irc = 0

end subroutine Cho_SetShP2RS_2
