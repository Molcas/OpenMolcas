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

subroutine Make_Fluctuating_Charges(nAtoms,iANr,nij,nPert,rMP,nElem,EC,Alpha)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, iANr(nAtoms), nij, nPert, nElem
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1)
real(kind=wp), intent(in) :: EC(3,nij), Alpha
real(kind=wp), allocatable :: A(:,:), AInv(:,:), dQ(:), Lambda(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! For the A-matrix

call mma_allocate(AInv,nAtoms,nAtoms,label='AInv')
call mma_allocate(A,nAtoms,nAtoms,label='A')

call Build_AMatrix(nAtoms,iANr,A,AInv,EC,nij,Alpha)

call mma_deallocate(A)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the fluctuating charges

call mma_allocate(Lambda,nAtoms,label='Lambda')
call mma_allocate(dQ,nAtoms,label='dQ')

call Fluctuating(AInv,nAtoms,Lambda,dQ,nij,nPert,iANr,rMP,nElem,EC,Alpha)

call mma_deallocate(Lambda)
call mma_deallocate(dQ)
call mma_deallocate(AInv)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Make_Fluctuating_Charges
