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

subroutine DIAG_R2_RASSI(MATRIX,NBTOT,INFO,W1,Z1)
! THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A REAL SQUARE
! MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
! ARE DIRECTED INTO W1 AND THE REAL EIGENVECTORS ARE WRITTEN TO Z1.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NBTOT, INFO
real(kind=wp) :: MATRIX(NBTOT,NBTOT), W1(NBTOT), Z1(NBTOT,NBTOT)
integer(kind=iwp) :: I, J
real(kind=wp), allocatable :: AP(:), WORK(:), W(:), Z(:,:)

call mma_allocate(AP,Nbtot*(Nbtot+1)/2,Label='AP')
call mma_allocate(WORK,3*Nbtot,Label='WORK')
call mma_allocate(W,Nbtot,Label='W')
call mma_allocate(Z,Nbtot,Nbtot,Label='Z')

! initializations
INFO = 0

do j=1,Nbtot
  do i=1,j
    AP(i+(j-1)*j/2) = MATRIX(i,j)
  end do
end do

call dspev_('V','U',Nbtot,AP,W,Z,Nbtot,WORK,INFO)
do I=1,Nbtot
  W1(I) = W(I)
  do J=1,Nbtot
    Z1(J,I) = Z(J,I)
  end do
end do

call mma_deallocate(AP)
call mma_deallocate(WORK)
call mma_deallocate(W)
call mma_deallocate(Z)

end subroutine DIAG_R2_RASSI
