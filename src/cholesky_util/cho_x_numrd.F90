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
!  Cho_X_NumRd
!
!> @brief
!>   Return the number of Cholesky vectors that may be read into \p Mem words of memory.
!> @author Thomas Bondo Pedersen
!>
!> @details
!> The count starts at vector \p iVec1 of symmetry \p iSym
!> (this is needed since the vectors are stored in
!> different reduces sets).
!> On exit, \p Cho_X_NumRd is negative if some error has
!> occurred (``-1``, ``-2``, and ``-3`` signify errors in input
!> variables, ``-4`` indicates an error in \p Cho_X_SetRed).
!>
!> @note
!> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
!>
!> @param[in]     iVec1 First vector
!> @param[in]     iSym  Symmetry
!> @param[in,out] iRedC Reduced set in core (location ``3``); ``0`` (or ``-1``) if unknown or undefined
!> @param[in]     Mem   Memory available for read
!***********************************************************************

function Cho_X_NumRd(iVec1,iSym,iRedC,Mem)

use Cholesky, only: InfVec, MaxVec, nDimRS, nnBstR, nSym, NumCho
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_X_NumRd
integer(kind=iwp), intent(in) :: iVec1, iSym, Mem
integer(kind=iwp), intent(inout) :: iRedC
integer(kind=iwp) :: iLoc, irc, iRed, iVec, Need, NumRd

if ((iSym < 1) .or. (iSym > nSym)) then
  Cho_X_NumRd = -1
else if (NumCho(iSym) < 1) then
  Cho_X_NumRd = 0
else if (NumCho(iSym) > MaxVec) then
  Cho_X_NumRd = -2
else if ((iVec1 < 1) .or. (iVec1 > NumCho(iSym))) then
  Cho_X_NumRd = -3
else if (Mem < 1) then
  Cho_X_NumRd = 0
else
  NumRd = 0
  Need = 0
  iVec = iVec1-1
  if (.not. allocated(nDimRS)) then
    iLoc = 3
    do while ((iVec < NumCho(iSym)) .and. (Need < Mem))
      iVec = iVec+1
      iRed = InfVec(iVec,2,iSym)
      if (iRed /= iRedC) then
        irc = 0
        call Cho_X_SetRed(irc,iLoc,iRed)
        if (irc /= 0) then
          Cho_X_NumRd = -4
          iRedC = -1
          return
        end if
        iRedC = iRed
      end if
      Need = Need+nnBstR(iSym,iLoc)
      if (Need <= Mem) NumRd = NumRd+1
    end do
  else
    do while ((iVec < NumCho(iSym)) .and. (Need < Mem))
      iVec = iVec+1
      iRed = InfVec(iVec,2,iSym)
      Need = Need+nDimRS(iSym,iRed)
      if (Need <= Mem) NumRd = NumRd+1
    end do
  end if
  Cho_X_NumRd = NumRd
end if

end function Cho_X_NumRd
