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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************
!  Cho_X_Bookmark
!
!> @brief
!>   Get number of Cholesky vectors needed to achieve a given integral target accuracy \p Thr
!> @author Thomas Bondo Pedersen, August 2012
!>
!> @details
!> Return number of Cholesky vectors in each symmetry block needed to
!> achieve an integral accuracy of \p Thr (&le; decomposition threshold).
!> The actual integral accuracy in each symmetry block is returned in
!> array \p delta (for 1C-CD, the accuracy is judged by 1-center
!> diagonal elements only!). Note that \p mSym may be smaller than \p nSym
!> (so that one may find a bookmark in irrep 1 only, for example).
!> Available for full CD and 1C-CD, but *not RI*.
!>
!> On exit, return codes signify:
!>
!> - \p irc = ``-1``: bookmarks not available
!> - \p irc =  ``0``: all OK
!> - \p irc =  ``1``: illegal input
!>
!> @param[in]  Thr   Target integral accuracy
!> @param[in]  mSym  Number of irreps
!> @param[out] nVec  Number of Cholesky vectors
!> @param[out] delta Actual integral accuracy
!> @param[out] irc   Return code
!***********************************************************************

subroutine Cho_X_Bookmark(Thr,mSym,nVec,delta,irc)

use Cholesky, only: BkmThr, BkmVec, Cho_Real_Par, nRow_BkmThr, nSym, ThrCom
#ifdef _DEBUGPRINT_
use Cholesky, only: NumCho
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Thr
integer(kind=iwp), intent(in) :: mSym
integer(kind=iwp), intent(out) :: nVec(mSym), irc
real(kind=wp), intent(out) :: delta(mSym)
integer(kind=iwp) :: iRS, iSym, l, n
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: BkmScr(:)
logical(kind=iwp), parameter :: DebugPrint = .false.
character(len=*), parameter :: SecNam = 'Cho_X_Bookmark'

! Set return code.
! ----------------

irc = 0

! Check that bookmarks are available.
! -----------------------------------

if ((.not. allocated(BkmVec)) .or. (.not. allocated(BkmThr))) then
  irc = -1
  return
end if

! Check input.
! ------------

if ((sign(One,Thr) < Zero) .or. (Thr < ThrCom) .or. (mSym < 1) .or. (mSym > nSym)) then
  irc = 1
  return
end if

! Debug print.
! ------------

if (DebugPrint) then
  call Cho_Head(SecNam//': Bookmarks (nVec,delta)','-',80,u6)
  do iSym=1,nSym
    write(u6,'(A,I2,A)') 'Symmetry block',iSym,' Bookmarks (nVec,delta)'
    write(u6,'(5(1X,A,I6,A,ES15.8,A))') ('(',BkmVec(iRS,iSym),',',BkmThr(iRS,iSym),')',iRS=1,nRow_BkmThr)
  end do
end if

! Find largest accuracy below Thr.
! --------------------------------

do iSym=1,mSym
  iRS = 0
  Found = .false.
  do while ((iRS < nRow_BkmThr) .and. (.not. Found))
    iRS = iRS+1
    Found = BkmThr(iRS,iSym) <= Thr
  end do
  if (.not. Found) then
    call Cho_Quit('Bug detected in '//SecNam,104)
  else
    nVec(iSym) = BkmVec(iRS,iSym)
    delta(iSym) = BkmThr(iRS,iSym)
  end if
end do

! Parallel run: take into account the distribution of vectors
! across nodes. Set nVec to the number of vectors needed on this
! node (process).
if (Cho_Real_Par) then
  l = nVec(1)
  do iSym=2,mSym
    l = max(l,nVec(iSym))
  end do
  call mma_Allocate(BkmScr,l,Label='BkmScr')
  do iSym=1,mSym
    call Cho_P_Distrib_Vec(1,nVec(iSym),BkmScr,n)
    nVec(iSym) = n
  end do
  call mma_deallocate(BkmScr)
end if

#ifdef _DEBUGPRINT_
! Self test
do iSym=1,mSym
  if (nVec(iSym) > NumCho(iSym)) call Cho_Quit('Error in '//SecNam,104)
end do
#endif

end subroutine Cho_X_Bookmark
