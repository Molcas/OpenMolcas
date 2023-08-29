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
!  Cho_X_Final
!
!> @brief
!>   Finalize Cholesky utilities
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Deallocates memory, closes files, etc., as initialized
!> by ::Cho_X_Init. On exit, \p irc = ``0`` signals successful finalization.
!>
!> @param[out] irc Return code
!***********************************************************************

subroutine Cho_X_Final(irc)

use Cholesky, only: BkmThr, BkmVec, ChoIniCheck, MySP, nCol_BkmThr, nCol_BkmVec, nRow_BkmThr, nRow_BkmVec
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: ChoIsIni

! Set  error code.
! ----------------

irc = 0

! Read initialization integer flag from runfile.
! ----------------------------------------------

call Get_iScalar('ChoIni',ChoIsIni)

! Finalize if needed.
! -------------------

if (ChoIsIni == ChoIniCheck) then

  ! Close files.
  ! ------------

  call Cho_OpenVR(2,2)

  ! Deallocate vector buffer.
  ! -------------------------

  call Cho_VecBuf_Final()

  ! Deallocate memory.
  ! ------------------

  call Cho_X_Dealloc(irc)
  if (irc == 0) then

    if (allocated(MySP)) call mma_deallocate(MySP)
    if (allocated(BkmVec)) then
      call mma_deallocate(BkmVec)
      nRow_BkmVec = 0
      nCol_BkmVec = 0
    end if
    if (allocated(BkmThr)) then
      call mma_deallocate(BkmThr)
      nRow_BkmThr = 0
      nCol_BkmThr = 0
    end if

  end if

  ! Reset initialization integer on runfile to "not set".
  ! -----------------------------------------------------

  ChoIsIni = ChoIniCheck+1
  call Put_iScalar('ChoIni',ChoIsIni)

end if

end subroutine Cho_X_Final
