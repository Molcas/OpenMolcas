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
!  Cho_X_ReoVec
!
!> @brief
!>   Reorder Cholesky vectors on disk
!> @author Thomas Bondo Pedersen
!>
!> @details
!> @note
!> You should consider using on-the-fly read/reorder routine ::Cho_X_getVFull instead!!
!>
!> This routine makes it possible to reorder Cholesky vectors
!> from reduced storage to full storage in exactly the same
!> manner as when specifying the keywork ``REORder`` in the
!> Cholesky input section (in Seward). If the vectors have
!> already been reordered (checked through runfile), the routine
!> returns immediately.
!> The reordered vectors
!> are stored in the files ``CHFVnm`` where ``n`` is the symmetry of
!> the first AO index, ``m`` that of the second. The resulting files
!> are split (if needed).
!>
!> @note
!> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
!>
!> @param[out] irc return code
!***********************************************************************

subroutine Cho_X_ReoVec(irc)

use Cholesky, only: nnBstRT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: iReo, l_Wrk
integer(kind=iwp), allocatable :: iRS2F(:,:)
real(kind=wp), allocatable :: Wrk(:)

irc = 0

call Get_iScalar('Cholesky Reorder',iReo)
if (iReo == 0) then
  call mma_allocate(iRS2F,3,nnBstRT(1),Label='iRS2F')
  call mma_maxDBLE(l_Wrk)
  call mma_allocate(Wrk,l_Wrk,Label='Wrk')
  call Cho_ReoVec(iRS2F,3,nnBstRT(1),Wrk,l_Wrk)
  call mma_deallocate(Wrk)
  call mma_deallocate(iRS2F)
  iReo = 1
  call Put_iScalar('Cholesky Reorder',iReo)
end if

end subroutine Cho_X_ReoVec
