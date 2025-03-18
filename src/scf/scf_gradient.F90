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

subroutine SCF_Gradient()
! Compute and store the (local) gradient at the current orbitals and iteration

use LnkLst, only: LLlGrd, PutVec
use InfSCF, only: CMO, Iter, FockMO, kOV, mOV, nBO, nBT, nD, nOO, OneHam, Ovrlp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
real(kind=wp), allocatable :: GrdOO(:,:), GrdOV(:)

call mma_allocate(GrdOO,nOO,nD,Label='GrdOO')
call mma_allocate(GrdOV,mOV,Label='GrdOV')

! Compute the last gradient
call Mk_FockMO(OneHam,Ovrlp,nBT,CMO,nBO,FockMO,nOO,nD,Iter)
call EGrad(OneHam,Ovrlp,nBT,CMO,nBO,GrdOO,nOO,nD,Iter)
call vOO2OV(GrdOO,nOO,GrdOV,mOV,nD,kOV)

! Write gradient to linked list
call PutVec(GrdOV,mOV,Iter,'OVWR',LLlGrd)

call mma_deallocate(GrdOO)
call mma_deallocate(GrdOV)

end subroutine SCF_Gradient
