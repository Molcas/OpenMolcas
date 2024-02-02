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
!   this routine performs the diagonalization of a Complex square
!   matrix with the dimension nbtot. the eigenvalues of the diagonalization
!   are directed into w1 and the Complex eigenvectors are written to z1.

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer :: info, i, j, n
complex(kind=8), intent(in) :: matrix(n,n)
complex(kind=8), intent(out) :: z(n,n)
real(kind=8), intent(out) :: w(n)
! local variables:
real(kind=8), allocatable :: rwork(:) !rwork(3*n-2)
real(kind=8), allocatable :: w1(:)    !w1(n)
complex(kind=8), allocatable :: ap(:)    !ap(n*(n+1)/2)
complex(kind=8), allocatable :: work(:)  !work(2*n-1)
complex(kind=8), allocatable :: z1(:,:)  !z1(n,n)
real(kind=8), external :: dznrm2_
real(kind=8) :: RM

info = 0
call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,Z,1)
call dcopy_(N,[0.0_wp],0,W,1)

RM = 0.0_wp
RM = dznrm2_(n*n,matrix,1)

if (RM > 0.0_wp) then
  call mma_allocate(ap,(n*(n+1)/2),'ap')
  call mma_allocate(work,(2*n-1),'work')
  call mma_allocate(z1,n,n,'work')
  call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,AP,1)
  call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,work,1)
  call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,Z1,1)

  call mma_allocate(rwork,(3*n-2),'rwork')
  call mma_allocate(w1,n,'w1')
  call dcopy_(3*N-2,[0.0_wp],0,rwork,1)
  call dcopy_(N,[0.0_wp],0,W1,1)

  do j=1,n
    do i=1,j
      ap(i+(j-1)*j/2) = matrix(i,j)
    end do
  end do
  ! diagonalize:
  call zhpev_('v','u',n,ap,w1,z1,n,work,rwork,info)
  ! save results:
  call dcopy_(N,W1,1,W,1)
  call zcopy_(N*N,Z1,1,Z,1)

  call mma_deallocate(rwork)
  call mma_deallocate(w1)
  call mma_deallocate(ap)
  call mma_deallocate(work)
  call mma_deallocate(z1)
else
  ! return dummy results:
  do i=1,n
    w(i) = 0.0_wp
    z(i,i) = (1.0_wp,0.0_wp)
  end do
end if

return

end subroutine diag_c2
