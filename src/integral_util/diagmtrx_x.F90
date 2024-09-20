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

subroutine DiagMtrx_x(H,nH,iNeg)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(inout) :: H(nH,nH)
integer(kind=iwp), intent(out) :: iNeg
integer(kind=iwp) :: i, ij, j
real(kind=wp) :: SumHii, Temp
real(kind=wp), allocatable :: Diag(:,:), EVal(:), EVec(:,:), HU(:,:)

call mma_allocate(EVal,nTri_Elem(nH),label='EVal')
call mma_allocate(EVec,nH,nH,label='EVec')

! Copy elements for H

SumHii = Zero
ij = 0
do i=1,nH
  do j=1,i
    ij = ij+1
    EVal(ij) = H(i,j)
  end do
  SumHii = SumHii+H(i,i)
end do
!write(u6,*) ' SumHii=',SumHii

! Set up a unit matrix

call unitmat(EVec,nH)

! Compute eigenvalues and eigenvectors

call Jacob(EVal,EVec,nH,nH)
call Jacord(EVal,EVec,nH,nH)

! Print out the result

iNeg = 0
do i=1,nH
  if (EVal(nTri_Elem(i)) < Zero) iNeg = iNeg+1
end do

call mma_allocate(Diag,nH,nH,label='Diag')
call mma_allocate(HU,nH,nH,label='HU')

Diag(:,:) = Zero
do i=1,nH
  temp = EVal(nTri_Elem(i))
  !write(u6,'(A,G10.4)') 'Hii=',temp
  Diag(i,i) = max(abs(temp),1.0e-15_wp)
end do

call DGEMM_('N','N',nH,nH,nH,One,EVec,nH,Diag,nH,Zero,HU,nH)
call DGEMM_('N','T',nH,nH,nH,One,HU,nH,EVec,nH,Zero,H,nH)

call mma_deallocate(HU)
call mma_deallocate(Diag)
call mma_deallocate(EVec)
call mma_deallocate(EVal)

return

end subroutine DiagMtrx_x
