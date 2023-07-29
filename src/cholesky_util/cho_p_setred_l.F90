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

subroutine Cho_P_SetRed_L()
!
! Purpose: set next local reduced set. The next global reduced set
!          must be available at (global) index array location 2.

use Cholesky, only: iiBstR, iiBstR_G, iiBstRSh, iiBstRSh_G, iL2G, IndRed, IndRed_G, LuPri, MySP, nnBstR, nnBstRSh, nnBstRSh_G, &
                    nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, i0, iC, irc, iShlAB, iSP, iSym, j, k, k0, l, ll, nDim
character(len=*), parameter :: SecNam = 'Cho_P_SetRed_L'

! Copy current local reduced set (at location 2) to location 3.
! -------------------------------------------------------------

call Cho_X_RSCopy(irc,2,3)
if (irc /= 0) then
  write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
  call Cho_Quit('Error in '//SecNam,104)
end if

! Re-initialize local reduced set indices at location 2.
! ------------------------------------------------------

nDim = nSym*nnShl
IndRed(:,2) = 0
call iZero(iiBstRSh(:,:,2),nDim)
call iZero(nnBstRSh(:,:,2),nDim)
call iZero(iiBstR(:,2),nSym)
call iZero(nnBstR(:,2),nSym)
nnBstRT(2) = 0

! Set local nnBstRSh counter at location 2.
! -----------------------------------------

do iSP=1,nnShl
  iShlAB = mySP(iSP)
  do iSym=1,nSym
    nnBstRSh(iSym,iSP,2) = nnBstRSh_G(iSym,iShlAB,2)
  end do
end do

! Set remaining reduced set indices (excl. IndRed), location 2.
! -------------------------------------------------------------

call Cho_SetRedInd(2)

! Set local IndRed to point to local 1st reduced set.
! ---------------------------------------------------

iC = 0
do iSym=1,nSym
  do iSP=1,nnShl
    iShlAB = mySP(iSP)
    i0 = iiBstR_G(iSym,2)+iiBstRSh_G(iSym,iShlAB,2)
    k0 = iiBstR(iSym,3)+iiBstRSh(iSym,iSP,3)
    do i=1,nnBstRSh_G(iSym,iShlAB,2)
      j = IndRed_G(i0+i,2) ! addr in global rs1
      iC = iC+1
      k = 0
      do while (k < nnBstRSh(iSym,iSP,3))
        k = k+1
        ll = IndRed(k0+k,3)
        l = iL2G(ll)
        if (l == j) then
          IndRed(iC,2) = ll
          k = nnBstRSh(iSym,iSP,3) ! break while loop
        end if
      end do
    end do
  end do
end do

end subroutine Cho_P_SetRed_L
