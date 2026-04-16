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

subroutine DIAG_R2_RASSI(MATRIX,NBTOT,INFO,W,Z)
! THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A REAL SQUARE
! MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
! ARE DIRECTED INTO W AND THE REAL EIGENVECTORS ARE WRITTEN TO Z.

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBTOT
real(kind=wp), intent(in) :: MATRIX(NBTOT,NBTOT)
integer(kind=iwp), intent(out) :: INFO
real(kind=wp), intent(out) :: W(NBTOT), Z(NBTOT,NBTOT)
integer(kind=iwp) :: J
real(kind=wp), allocatable :: AP(:), WORK(:)

call mma_allocate(AP,nTri_Elem(Nbtot),Label='AP')
call mma_allocate(WORK,3*Nbtot,Label='WORK')

! initializations
INFO = 0

do j=1,Nbtot
  AP(nTri_Elem(j-1)+1:nTri_Elem(j)) = MATRIX(1:j,j)
end do

call dspev_('V','U',Nbtot,AP,W,Z,Nbtot,WORK,INFO)

call mma_deallocate(AP)
call mma_deallocate(WORK)

end subroutine DIAG_R2_RASSI
