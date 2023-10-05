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

subroutine Cho_SetShP2Q_2(irc,iLoc,iShlAB,nAB)
!
! Purpose: set mapping from shell pair iShlAB to qualified
!          columns within current reduced set (stored at location
!          iLoc = 2 or 3).
!          If a non-zero code (irc) is returned, nothing has been
!          set!!

use Index_Functions, only: nTri_Elem
use Cholesky, only: IndRed, IndRSh, iQuAB, iShP2Q, iSP2F, nBstSh, nQual, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc, nAB(nSym)
integer(kind=iwp), intent(in) :: iLoc, iShlAB
integer(kind=iwp) :: iAB, iQ, iShlA, iShlB, iSym, jAB, kAB, kShlAB, l_iShP2Q, lTst, NumAB

! Check allocations.
! ------------------

call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
if (iShlA == iShlB) then
  NumAB = nTri_Elem(nBstSh(iShlA))
else
  NumAB = nBstSh(iShlA)*nBstSh(iShlB)
end if
lTst = 2*NumAB
l_iShP2Q = 0
if (allocated(iShP2Q)) l_iShP2Q = size(iShP2Q)
if ((l_iShP2Q < 1) .or. (l_iShP2Q < lTst)) then
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
! iShP2Q(1,AB) = index among qualified, symmetry reduced.
! iShP2Q(2,AB) = symmetry block.
! Zeros are returned if the element AB is not qualified.
! -------------------------------------------------------

iShP2Q(:,1:NumAB) = 0
nAB(:) = 0

do iSym=1,nSym
  do iQ=1,nQual(iSym)
    iAB = iQuAB(iQ,iSym) ! addr in curr. red. set
    jAB = IndRed(iAB,iLoc)  ! addr in 1st red. set
    kShlAB = IndRSh(jAB)  ! shell pair (full)
    if (kShlAB == iSP2F(iShlAB)) then
      kAB = IndRed(jAB,1) ! addr in full shell pair
      iShP2Q(1,kAB) = iQ
      iShP2Q(2,kAB) = iSym
      nAB(iSym) = nAB(iSym)+1
    end if
  end do
end do

! Set return code 0: all ok!
! --------------------------

irc = 0

end subroutine Cho_SetShP2Q_2
