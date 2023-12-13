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

subroutine Print_EigenValues(H,nH)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(in) :: H(nTri_Elem(nH))
integer(kind=iwp) :: i
real(kind=wp), allocatable :: EVal(:), EVec(:,:)

call mma_allocate(EVal,nTri_Elem(nH),label='EVal')
call mma_allocate(EVec,nH,nH,label='EVec')

! Copy elements for H

EVal(:) = H

! Set up a unit matrix

call unitmat(EVec,nH)

! Compute eigenvalues and eigenvectors

call Jacob(EVal,EVec,nH,nH)
call Jacord(EVal,EVec,nH,nH)

! Print out the result

write(u6,*)
write(u6,*) 'Eigenvalues of the matrix'
write(u6,*)
write(u6,'(10F15.8)') (EVal(nTri_Elem(i)),i=1,nH)

call mma_deallocate(EVal)
call mma_deallocate(EVec)

return

end subroutine Print_EigenValues
