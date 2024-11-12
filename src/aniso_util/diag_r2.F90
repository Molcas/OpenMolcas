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

subroutine DIAG_R2(MATRIX,N,INFO,W,Z)
! THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A Real SQUARE
! MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
! ARE DIRECTED INTO W1 AND THE Real EIGENVECTORS ARE WRITTEN TO Z1.

use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: MATRIX(N,N)
integer(kind=iwp), intent(out) :: INFO
real(kind=wp), intent(out) :: W(N), Z(N,N)
integer(kind=iwp) :: I, J
real(kind=wp), allocatable :: AP(:), WORK(:)

! initializations
INFO = 0
if (N < 1) return

W(:) = Zero
Z(:,:) = Zero

call mma_allocate(AP,nTri_Elem(N),'AP')
call mma_allocate(WORK,3*N,'WORK')
WORK(:) = Zero

do j=1,N
  do i=1,j
    AP(iTri(i,j)) = MATRIX(i,j)
  end do
end do
! diagonalize:
call DSPEV_('V','U',N,AP,W,Z,N,WORK,INFO)

call mma_deallocate(AP)
call mma_deallocate(WORK)

return

end subroutine DIAG_R2
