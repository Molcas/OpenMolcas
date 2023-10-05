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

subroutine Cho_P_WrDiag()
!
! Purpose: store global diagonal on disk (Parallel only).
!          NB: on exit, initial global diagonal is stored!

use Cholesky, only: Cho_Real_Par, Diag_G, nnBstRT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
real(kind=wp), allocatable :: Diag_L(:)

if (Cho_Real_Par) then
  call mma_allocate(Diag_L,nnBstRT(1),Label='Diag_L')
  call Cho_IODiag(Diag_L,2)
  call Cho_P_SyncDiag(Diag_L,1)
  call Cho_P_IndxSwp()
  call Cho_IODiag(Diag_G,1)
  call Cho_P_IndxSwp()
  call mma_deallocate(Diag_L)
end if

end subroutine Cho_P_WrDiag
