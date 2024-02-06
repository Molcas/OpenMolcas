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

use Constants, only: Zero

implicit none
#include "stdalloc.fh"
integer, intent(in) :: N
integer, intent(out) :: INFO
real(kind=8), intent(in) :: MATRIX(N,N)
real(kind=8), intent(out) :: W(N), Z(N,N)
! local variables:
integer :: I, J
real(kind=8), allocatable :: AP(:)   !(N*(N+1)/2)
real(kind=8), allocatable :: WORK(:) !(3*N)
real(kind=8), allocatable :: W1(:)   !(N)
real(kind=8), allocatable :: Z1(:,:) !(N,N)

! initializations
INFO = 0
if (N < 1) return

call dcopy_(N,[Zero],0,W,1)
call dcopy_(N*N,[Zero],0,Z,1)

call mma_allocate(AP,(N*(N+1)/2),'AP')
call mma_allocate(WORK,3*N,'WORK')
call mma_allocate(W1,N,'W1')
call mma_allocate(Z1,N,N,'Z1')
call dcopy_(N*(N+1)/2,[Zero],0,AP,1)
call dcopy_(3*N,[Zero],0,WORK,1)
call dcopy_(N,[Zero],0,W1,1)
call dcopy_(N*N,[Zero],0,Z1,1)

do j=1,N
  do i=1,j
    AP(i+(j-1)*j/2) = MATRIX(i,j)
  end do
end do
! diagonalize:
call DSPEV_('V','U',N,AP,W1,Z1,N,WORK,INFO)

call dcopy_(N,W1,1,W,1)
call dcopy_(N*N,Z1,1,Z,1)

call mma_deallocate(AP)
call mma_deallocate(WORK)
call mma_deallocate(W1)
call mma_deallocate(Z1)

return

end subroutine DIAG_R2
