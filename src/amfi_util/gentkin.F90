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

subroutine gentkin(L,TKIN,nprims,exponents,rootOVLPinv)
!bs subroutine to generate the kinetic energy

use AMFI_global, only: MxprimL
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, nprims
real(kind=wp), intent(out) :: TKIN(nprims,nprims)
real(kind=wp), intent(in) :: exponents(nprims), rootOVLPinv(MxprimL,MxprimL)
integer(kind=iwp) :: irun1, irun2
real(kind=wp), allocatable :: dummy(:,:), dummy2(:,:)
real(kind=wp), external :: Tkinet

call mma_allocate(dummy,nprims,nprims,label='dummy')
call mma_allocate(dummy2,nprims,nprims,label='dummy2')

!bs build the symmetric matrix
do irun2=1,nprims
  do irun1=1,irun2
    dummy(irun1,irun2) = Tkinet(L,exponents(irun1),exponents(irun2))
    dummy(irun2,irun1) = dummy(irun1,irun2)
  end do
end do
!bs now transform by rootOVLPinv*dummy*rootOVLPinv
call dgemm_('N','N',nprims,nprims,nprims,One,dummy,nprims,rootOVLPinv,MxprimL,Zero,dummy2,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,rootOVLPinv,MxprimL,dummy2,nprims,Zero,Tkin,nprims)

call mma_deallocate(dummy)
call mma_deallocate(dummy2)

return

end subroutine gentkin
