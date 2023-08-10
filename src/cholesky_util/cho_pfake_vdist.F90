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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_PFake_VDist()
!
! Author: Thomas Bondo Pedersen, April 2007.
! Purpose: fake parallel distribution of vectors.

use Para_Info, only: Is_Real_Par, nProcs
use Cholesky, only: Cho_Fake_Par, InfVec, MaxVec, nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iRedC, iSym, iV, l_Wrk, nRead, nV
integer(kind=iwp), allocatable :: IDV(:), InfV(:,:)
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'Cho_PFake_VDist'

! Return if nothing to do.
! ------------------------

if ((nProcs == 1) .or. (.not. Is_Real_Par())) return ! serial run
if (.not. CHO_FAKE_PAR) return ! true parallel run

! Memory allocation.
! ------------------

call mma_allocate(InfV,2,MaxVec+1,Label='InfV')
call mma_allocate(IDV,MaxVec,Label='IDV')
call mma_maxDBLE(l_Wrk)
call mma_allocate(Wrk,l_Wrk,Label='Wrk')

! Distribute vectors in each symmetry:
! read from LuCho files; put into LuCho_G files.
! ----------------------------------------------

iRedC = -1
do iSym=1,nSym
  InfV(:,:) = 0
  nV = 0
  call Cho_Distrib_Vec(1,NumCho(iSym),IDV,nV)
  iV = 0
  do while (iV < nV)
    nRead = 0
    call Cho_PFake_GetVec(Wrk,size(Wrk),IDV(iV+1),nV-iV,InfV(:,iV+1),iSym,nRead,iRedC)
    if (nRead < 1) call Cho_Quit('Insufficient memory in '//SecNam,101)
    call Cho_PFake_PutVec(Wrk,InfV,nRead,iSym,iV+1)
    iV = iV+nRead
  end do
# ifdef _DEBUGPRINT_
  if (iV /= nV) call Cho_Quit('Logical error in '//SecNam,103)
# endif
  InfVec(1:nV,3,iSym) = InfV(2,1:nV)
end do

! Deallocation.
! -------------

call mma_deallocate(Wrk)
call mma_deallocate(IDV)
call mma_deallocate(InfV)

end subroutine Cho_PFake_VDist
