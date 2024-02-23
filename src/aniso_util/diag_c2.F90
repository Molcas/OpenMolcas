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

subroutine diag_c2(matrix,n,info,w,z)
! this routine performs the diagonalization of a Complex square
! matrix with the dimension nbtot. the eigenvalues of the diagonalization
! are directed into w1 and the Complex eigenvectors are written to z1.

use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: matrix(n,n)
integer(kind=iwp), intent(out) :: info
real(kind=wp), intent(out) :: w(n)
complex(kind=wp), intent(out) :: z(n,n)
integer(kind=iwp) :: i, j
real(kind=wp) :: RM
real(kind=wp), allocatable :: rwork(:)
complex(kind=wp), allocatable :: ap(:), work(:)
real(kind=wp), external :: dznrm2_

info = 0
W(:) = Zero
Z(:,:) = cZero

RM = dznrm2_(n*n,matrix,1)

if (RM > Zero) then
  call mma_allocate(ap,nTri_Elem(n),'ap')
  call mma_allocate(work,2*n-1,'work')
  call mma_allocate(rwork,3*n-2,'rwork')
  work(:) = cZero
  rwork(:) = Zero

  do j=1,n
    do i=1,j
      ap(iTri(i,j)) = matrix(i,j)
    end do
  end do
  ! diagonalize:
  call zhpev_('v','u',n,ap,w,z,n,work,rwork,info)

  call mma_deallocate(rwork)
  call mma_deallocate(ap)
  call mma_deallocate(work)
else
  ! return dummy results:
  do i=1,n
    z(i,i) = cOne
  end do
end if

return

end subroutine diag_c2
