!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  trg_SVD
!
!> @brief  Transform a gradient to a rotated reference
!> @author Ignacio Fdez. Galvan
!>
!> @details
!> Transforms a gradient computed with respect to the rotation parameters
!> at the "current" orbitals to the rotation parameters at the "reference"
!> orbitals.
!>
!> @see ::trg_series
!>
!> @param[in]     Nv   Number of virtual orbitals
!> @param[in]     No   Number of occupied orbitals
!> @param[in]     X    Rotation parameters from reference to current
!> @param[in,out] G    Gradient at the current orbitals, will be replaced with the transformed gradient
!***********************************************************************

subroutine trg_svd(Nv,No,X,G)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: Nv, No
real(kind=wp), intent(in) :: X(Nv,No)
real(kind=wp), intent(inout) :: G(Nv,No)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: aux
real(kind=wp), allocatable :: D(:,:,:), G1(:,:), P(:,:), Q(:,:), R(:,:), s(:), tmp(:,:), W(:,:), Z(:,:)
real(kind=wp), parameter :: thrs = sqrt(epsilon(thrs))

! This applies the transformation
!
!   G_trans = R [(Z+Z^T) (*) D1 + (Z-Z^T) (*) D2 ] Q^T +
!             + (R sinc(s) R^T) G (I - Q Q^T)
!
! where
!
!   X = R s Q^T  (SVD of X)
!   Z = R^T G Q
!   D1_ij = sinc(s_i-s_j)/2     D2_ij = sinc(s_i+s_j)/2
!   (*) denotes an element-wise product
!   sinc(x) = sin(x)/x ; sinc(0) = 1
!
! for the special case where G and X are both antisymmetric with a
! single nonzero off-diagonal block.

k = min(No,Nv)
if (k < 1) return

call mma_allocate(Q,Nv,k,label='Q')
call mma_allocate(R,k,No,label='R')
call mma_allocate(s,k,label='s')

! Do the SVD, note the R matrix is already transposed
call mma_allocate(G1,Nv,No,label='G1')
G1(:,:) = X(:,:)
call large_svd(Nv,No,G1,Q,R,s)

call mma_allocate(D,k,k,2,label='D')
do j=1,k
  do i=1,j
    aux = abs(s(i)-s(j))
    if (aux < thrs) then
      D(i,j,1) = Half
    else
      D(i,j,1) = Half*sin(aux)/aux
    end if
    aux = abs(s(i)+s(j))
    if (aux < thrs) then
      D(i,j,2) = Half
    else
      D(i,j,2) = Half*sin(aux)/aux
    end if
    D(j,i,:) = D(i,j,:)
  end do
end do

! Z = Q^T G R
call mma_allocate(Z,k,k,label='Z')
call mma_allocate(tmp,k,k,label='tmp')
if (No <= Nv) then
  call dgemm_('T','N',No,No,Nv,One,Q,Nv,G,Nv,Zero,tmp,No)
  call dgemm_('N','T',No,No,No,One,tmp,No,R,No,Zero,Z,No)
else
  call dgemm_('N','T',Nv,Nv,No,One,G,Nv,R,Nv,Zero,tmp,Nv)
  call dgemm_('T','N',Nv,Nv,Nv,One,Q,Nv,tmp,Nv,Zero,Z,Nv)
end if
! C = (Z+Z^T) (*) D1 + (Z-Z^T) (*) D2
!   = Z (*) (D1+D2) + Z.T (*) (D1-D2)
tmp(:,:) = Z(:,:)*(D(:,:,1)+D(:,:,2))+transpose(Z(:,:))*(D(:,:,1)-D(:,:,2))
call mma_deallocate(D)
! G1 = Q C R^T
if (No <= Nv) then
  call dgemm_('N','N',No,No,No,One,tmp,No,R,No,Zero,Z,No)
  call dgemm_('N','N',Nv,No,No,One,Q,Nv,Z,No,Zero,G1,Nv)
else
  call dgemm_('N','N',Nv,Nv,Nv,One,Q,Nv,tmp,Nv,Zero,Z,Nv)
  call dgemm_('N','N',Nv,No,Nv,One,Z,Nv,R,Nv,Zero,G1,Nv)
end if
call mma_deallocate(Z)

call mma_allocate(W,k,k,label='W')
call mma_allocate(P,No+Nv-k,No+Nv-k,label='P')
call mma_allocate(Z,Nv,No,label='Z')
if (No <= Nv) then
  ! W = R sinc(s) R^T
  do i=1,k
    if (s(i) < thrs) then
      tmp(i,:) = R(i,:)
    else
      tmp(i,:) = sin(s(i))/s(i)*R(i,:)
    end if
  end do
  call dgemm_('T','N',No,No,No,One,R,No,tmp,No,Zero,W,No)
  ! P = I - Q Q^T
  call dgemm_('N','T',Nv,Nv,No,-One,Q,Nv,Q,Nv,Zero,P,Nv)
  do i=1,Nv
    P(i,i) = P(i,i)+One
  end do
  ! G = P G W
  call dgemm_('N','N',Nv,No,No,One,G,Nv,W,No,Zero,Z,Nv)
  call dgemm_('N','N',Nv,No,Nv,One,P,Nv,Z,Nv,Zero,G,Nv)
else
  ! W = Q sinc(s) Q^T
  do i=1,k
    if (s(i) < thrs) then
      tmp(:,i) = Q(:,i)
    else
      tmp(:,i) = sin(s(i))/s(i)*Q(:,i)
    end if
  end do
  call dgemm_('N','T',Nv,Nv,Nv,One,tmp,Nv,Q,Nv,Zero,W,Nv)
  ! P = I - R R^T
  call dgemm_('T','N',No,No,Nv,-One,R,Nv,R,Nv,Zero,P,No)
  do i=1,No
    P(i,i) = P(i,i)+One
  end do
  ! G = W G P
  call dgemm_('N','N',Nv,No,Nv,One,W,Nv,G,Nv,Zero,Z,Nv)
  call dgemm_('N','N',Nv,No,No,One,Z,Nv,P,No,Zero,G,Nv)
end if
! G = G + G1
G(:,:) = G(:,:)+G1(:,:)

call mma_deallocate(Q)
call mma_deallocate(R)
call mma_deallocate(s)
call mma_deallocate(W)
call mma_deallocate(P)
call mma_deallocate(Z)
call mma_deallocate(G1)
call mma_deallocate(tmp)

end subroutine trg_svd
