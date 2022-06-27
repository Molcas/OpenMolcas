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

subroutine x2c_makx(m,n,f,s,x)
! Make X matrix from m-dimensional (m=2n) Fock(f) and Overlap(s) matrix
!
! X is the relation(transfer) matrix of Large--Small component coefficients
! of electron solutions (positive energy solutions)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m, n
real(kind=wp), intent(in) :: f(m,m), s(m,m)
real(kind=wp), intent(out) :: x(n,n)
integer(kind=iwp) :: i, j, k, lwork, info
real(kind=wp), allocatable, target :: tmp(:), w(:), tF(:,:), tS(:,:)
real(kind=wp), pointer :: pF(:), pS(:)

lwork = 8*m
call mma_allocate(tF,m,m,label='TmpF')
call mma_allocate(tS,m,m,label='TmpS')
call mma_allocate(w,m,label='Eig')
call mma_allocate(tmp,lwork,label='Work')

! Copy Fock and Overlap matrix to temp arrays

tF(:,:) = f(:,:)
tS(:,:) = s(:,:)

! Diagonalization of Fock matrix with given overlap matrix

call dsygv_(1,'V','L',m,tF,m,tS,m,w,tmp,lwork,info)

! Calculate the X matrix from electron solutions

pF(1:m*m) => tF
pS(1:m*m) => tS
k = 1
do i=n+1,2*n
  do j=1,n
    ! store large component coefficients of electron solutions (matrix A)
    pF(k) = tF(j,i)
    ! store small component coefficients of electron solutions (matrix B)
    pS(k) = tF(n+j,i)
    k = k+1
  end do
end do
nullify(pF,pS)
! compute X=BA^{-1}
call XDR_dmatinv(tF,n)
call dmxma(n,'N','N',tS,tF,x,One)

! Free temp memories

call mma_deallocate(tF)
call mma_deallocate(tS)
call mma_deallocate(w)
call mma_deallocate(tmp)

return

end subroutine x2c_makx
