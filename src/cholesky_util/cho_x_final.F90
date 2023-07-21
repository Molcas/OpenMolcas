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

use ChoArr, only: MySP
use ChoBkm, only: BkmVec, BkmThr, nRow_BkmVec, nCol_BkmVec, nRow_BkmThr, nCol_BkmThr
use ChoIni
use stdalloc, only: mma_deallocate

implicit none
integer irc
character*11 SecNam
parameter(SecNam='Cho_X_Final')
integer ChoIsIni

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
  if (irc /= 0) Go To 1

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

  ! Reset initialization integer on runfile to "not set".
  ! -----------------------------------------------------

1 ChoIsIni = ChoIniCheck+1
  call Put_iScalar('ChoIni',ChoIsIni)

end if

end subroutine Cho_X_Final
