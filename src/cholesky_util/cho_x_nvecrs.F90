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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_nVecRS
!
!> @brief
!>   Find first vector and number of vectors in reduced set \p iRed, sym. block \p iSym
!> @author Thomas Bondo Pedersen
!>
!> @details
!> This routine finds the first vector and number of
!> vectors in reduced set \p iRed, sym. block \p iSym.
!> Note that \p iVec=\p nVec = ``0`` may be returned---this is
!> perfectly acceptable: a given reduced set may be
!> empty. However, if negative numbers are returned
!> (\p iVec < ``0`` and \p nVec < ``0``), an error has ocurred. This
!> should be tested by the caller!!
!>
!> @note
!> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
!>
!> @param[in]  iRed Reduced set
!> @param[in]  iSym Symmetry block (1-8)
!> @param[out] iVec First vector in red. set \p iRed, sym. \p iSym
!> @param[out] nVec Number of vectors in red. set \p iRed, sym. \p iSym
!***********************************************************************

subroutine Cho_X_nVecRS(iRed,iSym,iVec,nVec)

use Cholesky, only: InfVec, MaxVec, nSym, NumCho
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iRed, iSym
integer(kind=iwp), intent(out) :: iVec, nVec
integer(kind=iwp) :: irc, jRed, jVec, LastRed
logical(kind=iwp) :: Found
character(len=*), parameter :: SecNam = 'Cho_X_nVecRS'

! Check input.
! ------------

irc = 0
if ((iSym < 1) .or. (iSym > nSym)) irc = -1
if ((NumCho(iSym) < 0) .or. (NumCho(iSym) > MaxVec)) irc = -2
if (NumCho(iSym) == 0) then
  iVec = 0
  nVec = 0
  return
end if
LastRed = InfVec(NumCho(iSym),2,iSym)
if (LastRed < 1) irc = -3
if (iRed < 1) irc = -4
if (irc /= 0) then
  iVec = irc
  nVec = irc
  return
end if
if (iRed > LastRed) then
  iVec = 0
  nVec = 0
  return
end if
nVec = 0

! Find first vector in reduced set iRed.
! --------------------------------------

Found = .false.
jVec = 0
do while ((jVec < NumCho(iSym)) .and. (.not. Found))
  jVec = jVec+1
  jRed = InfVec(jVec,2,iSym)
  if (jRed == iRed) then
    iVec = jVec
    Found = .true.
  else if (jRed > iRed) then
    jVec = NumCho(iSym)  ! break loop
  end if
end do

! No first vector <=> 0 vectors in reduced set iRed.
! --------------------------------------------------

if (.not. Found) then
  iVec = 0
  nVec = 0
  return
end if

! Count number of vectors in reduced set iRed.
! --------------------------------------------

nVec = 1
jVec = iVec
do while (jVec < NumCho(iSym))
  jVec = jVec+1
  jRed = InfVec(jVec,2,iSym)
  if (jRed == iRed) then
    nVec = nVec+1
  else
    jVec = NumCho(iSym) ! break loop
  end if
end do

#ifdef _DEBUGPRINT_
! Debug: print result.
! --------------------

write(u6,*) SecNam,': there are ',nVec,' vectors in reduced set ',iRed,' (sym. block ',iSym,')'
write(u6,*) SecNam,': first vector is: ',iVec
#endif

end subroutine Cho_X_nVecRS
