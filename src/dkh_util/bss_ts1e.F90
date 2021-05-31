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

subroutine bss_ts1e(n,s,t,v,w,ul,us,clight)
! Evaluate the BSS Hamiltonian matrix and store the transform matrices

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
! w   aka pVp
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: s(n,n), t(n,n), w(n,n), clight
real(kind=wp), intent(inout) :: v(n,n)
real(kind=wp), intent(out) :: ul(n,n), us(n,n)
integer(kind=iwp) :: m, i, j, lwork, info
real(kind=wp), allocatable :: Tr(:,:), Bk(:,:), EL(:,:), ES(:,:), OL(:,:), OS(:,:), Ep(:), E0(:), KC(:,:), tF(:,:), W1(:), tmp(:), &
                              A(:,:), B(:,:), X(:,:), M1(:,:,:)

! Transform to the free-particle Foldy-Wothuysen picture

call mma_allocate(Tr,n,n,label='Tr') !IFG
call mma_allocate(Bk,n,n,label='Back')
call mma_allocate(EL,n,n,label='mEL')
call mma_allocate(ES,n,n,label='mES')
call mma_allocate(OL,n,n,label='mOL')
call mma_allocate(OS,n,n,label='mOS')
call mma_allocate(Ep,n,label='Ep')
call mma_allocate(E0,n,label='E0')
call mma_allocate(KC,n,3,label='KC') !IFG

call XDR_fpFW(n,s,t,v,w,Tr,Bk,EL,ES,OL,OS,Ep,E0,KC(:,1),KC(:,2),KC(:,3),clight)

! Diagonalize to get the X matrix

do i=1,n
  ! add (diagonal) kinetic matrix in fpFW picture
  EL(i,i) = EL(i,i)+E0(i)
  ES(i,i) = ES(i,i)-Ep(i)-clight*clight
end do
m = n+n
lwork = 8*m
call mma_allocate(tF,m,m,label='Fock')
call mma_allocate(W1,m,label='Eig')
call mma_allocate(tmp,lwork,label='Tmp')
do i=1,n
  do j=1,n
    tF(j,i) = EL(j,i)
    tF(j,i+n) = OL(j,i)
    tF(j+n,i) = OS(j,i)
    tF(j+n,i+n) = ES(j,i)
  end do
end do
call dsyev_('V','L',m,tF,m,W1,tmp,lwork,info)
call mma_deallocate(W1)
call mma_deallocate(tmp)
call mma_allocate(A,n,n,label='tmpA')
call mma_allocate(B,n,n,label='tmpB')
call mma_allocate(X,n,n,label='tmpX')
do i=1,n
  do j=1,n
    A(j,i) = tF(j,i+n)
    B(j,i) = tF(j+n,i+n)
  end do
end do
call XDR_dmatinv(A,n)
call dmxma(n,'N','N',B,A,X,One)

! Apply decoupling transformation to the Fock matrix

call dmxma(n,'C','N',X,OS,A,One)
call dmxma(n,'C','N',X,ES,B,One)
call dmxma(n,'N','N',B,X,ES,One)
call dmxma(n,'N','N',OL,X,B,One)
EL(:,:) = EL(:,:)+A(:,:)+B(:,:)+ES(:,:)
call dmxma(n,'C','N',X,X,A,One)
do i=1,n
  A(i,i) = A(i,i)+One
end do
! renormalization matrix
call XDR_dmatsqrt(A,n)
! decoupled electron Fock matrix in moment space (eigenfunction of T matrix)
call dmxma(n,'C','N',A,EL,B,One)
call dmxma(n,'N','N',B,A,v,One)

! Back transform to non-orthogonal basis picture

call dmxma(n,'C','N',Bk,v,B,One)
call dmxma(n,'N','N',B,Bk,v,One)

! Calculate transform matrices in non-orthogonal basis space

call dmxma(n,'N','N',X,A,B,One)
! A/B is the upper/lower part of transformation matrix in fpFW picture
call mma_allocate(M1,n,n,4,label='TmpM') !IFG
call XDR_mkutls(n,A,B,Tr,Bk,KC(:,1),KC(:,2),KC(:,3),ul,us,M1(:,:,1),M1(:,:,2),M1(:,:,3),M1(:,:,4))
call mma_deallocate(M1)

! Free temp memories

call mma_deallocate(Tr)
call mma_deallocate(Bk)
call mma_deallocate(EL)
call mma_deallocate(ES)
call mma_deallocate(OL)
call mma_deallocate(OS)
call mma_deallocate(Ep)
call mma_deallocate(E0)
call mma_deallocate(KC)
call mma_deallocate(tF)
call mma_deallocate(A)
call mma_deallocate(B)
call mma_deallocate(X)

return

end subroutine bss_ts1e
