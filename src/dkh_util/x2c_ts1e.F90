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

subroutine x2c_ts1e(n,s,t,v,w,ul,us,clight)
! Evaluate the X2C Hamiltonian matrix and store the transform matrices

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
! w   aka pVp
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: s(n,n), t(n,n), v(n,n), w(n,n)
real(kind=wp), intent(out) :: ul(n,n), us(n,n)
real(kind=wp), intent(in) :: clight
integer(kind=iwp) :: m, i, j
real(kind=wp) :: c_2c2, c_4c2, c_2c
real(kind=wp), allocatable :: F(:,:), St(:,:), X(:,:), A(:,:), B(:,:), C(:,:), SS(:,:)

! Construct the full dimensional (2*n) Fock matrix and overlap matrix

m = n+n
c_2c = clight*Two
c_2c2 = clight*clight*Two
c_4c2 = c_2c2*Two
! rescale, in order to make X matrix close to unit matrix
w(:,:) = w(:,:)/c_4c2
call mma_allocate(F,m,m,label='TmpF')
call mma_allocate(St,m,m,label='TmpS')
St(:,:) = Zero
do i=1,n
  do j=1,n
    St(j,i) = s(j,i)
    St(n+j,n+i) = t(j,i)/c_2c2
    F(i,j) = v(j,i)
    F(n+i,j) = t(j,i)
    F(i,n+j) = t(j,i)
    F(n+i,n+j) = w(j,i)-t(j,i)
  end do
end do

! Call diagonalization routine to obtain the X matrix

call mma_allocate(X,n,n,label='TmpX')
call x2c_makx(m,n,F,St,X)

! Calculate transformed Hamiltonian matrix

call mma_allocate(A,n,n,label='TmpA')
call mma_allocate(B,n,n,label='TmpB')
call mma_allocate(C,n,n,label='TmpC')
call mma_allocate(SS,n,n,label='TmpC')
call dmxma(n,'C','N',X,t,A,One)
call dmxma(n,'N','N',t,X,B,One)
call dmxma(n,'N','N',A,X,C,One)
! X-projected overlap matrix
SS(:,:) = s(:,:)+C(:,:)/c_2c2
! X-projected kinetic matrix
t(:,:) = A(:,:)+B(:,:)-C(:,:)
call XDR_dmatsqrt(s,n)
call dmxma(n,'C','N',s,SS,A,One)
call dmxma(n,'N','N',A,s,B,One)
call XDR_dmatsqrt(B,n)
call dmxma(n,'N','N',s,B,C,One)
call XDR_dmatinv(s,n)
! renormalization matrix, also the upper part of transformation matrix
call dmxma(n,'N','N',C,s,ul,One)
! lower part of the transformation matrix
call dmxma(n,'N','N',X,ul,us,One)

! Apply transformation to kinetic and potential matrices

call dmxma(n,'C','N',ul,t,A,One)
call dmxma(n,'N','N',A,ul,t,One)
call dmxma(n,'C','N',ul,v,A,One)
call dmxma(n,'N','N',A,ul,v,One)
call dmxma(n,'C','N',us,w,A,One)
call dmxma(n,'N','N',A,us,w,One)
v(:,:) = v(:,:)+t(:,:)+w(:,:)
! since pVp was rescaled, the lower part of transformation matrix also need to be rescaled
us(:,:) = us(:,:)/c_2c

! Free temp memories

call mma_deallocate(F)
call mma_deallocate(St)
call mma_deallocate(X)
call mma_deallocate(A)
call mma_deallocate(B)
call mma_deallocate(C)
call mma_deallocate(SS)

return

end subroutine x2c_ts1e
