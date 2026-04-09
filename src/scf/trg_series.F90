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
!  trg_series
!
!> @brief  Transform a gradient to a rotated reference, using a series expansion
!> @author Ignacio Fdez. Galvan
!>
!> @details
!> Transforms a gradient computed with respect to the rotation parameters
!> at the "current" orbitals to the rotation parameters at the "reference"
!> orbitals.
!>
!> @see ::trg_svd
!>
!> @param[in]     Nv   Number of virtual orbitals
!> @param[in]     No   Number of occupied orbitals
!> @param[in]     X    Rotation parameters from reference to current
!> @param[in,out] G    Gradient at the current orbitals, will be replaced with the transformed gradient
!***********************************************************************

subroutine trg_series(Nv,No,X,G)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: Nv, No
real(kind=wp), intent(in) :: X(Nv,No)
real(kind=wp), intent(inout) :: G(Nv,No)
integer(kind=iwp) :: i, j, k, term
real(kind=wp) :: aux
real(kind=wp), allocatable :: M(:,:), P(:,:), Q(:,:), tmp(:,:)
integer(kind=iwp), parameter :: maxterm = 30
real(kind=wp), parameter :: thrs = epsilon(thrs)

! This applies the transformation
!
!   G_trans = G + [X,G]/2! + [X,[X,G]]/3! + [X,[X,[X,G]]]/4! + ...
!
! for the special case where G and X are both antisymmetric with a
! single nonzero off-diagonal block. Only the odd terms are relevant.

k = min(No,Nv)
if (k < 1) return

call mma_allocate(P,k,k,label='P')
call mma_allocate(M,k,k,label='M')
call mma_allocate(Q,Nv,No,label='Q')
call mma_allocate(tmp,Nv,No,label='tmp')
if (No <= Nv) then
  ! P = X^T X
  call dgemm_('T','N',No,No,Nv,One,X,Nv,X,Nv,Zero,P,No)
else
  ! P = X X^T
  call dgemm_('N','T',Nv,Nv,No,One,X,Nv,X,Nv,Zero,P,Nv)
end if

Q(:,:) = G(:,:)
do term=1,maxterm
  if (No <= Nv) then
    ! M = Q^T X
    call dgemm_('T','N',No,No,Nv,One,Q,Nv,X,Nv,Zero,M,No)
  else
    ! M = X Q^T
    call dgemm_('N','T',Nv,Nv,No,One,X,Nv,Q,Nv,Zero,M,Nv)
  end if
  ! M = 2M - M^T = M + (M - M^T)
  do j=1,k
    do i=1,j-1
      aux = M(j,i)-M(i,j)
      M(j,i) = M(j,i)+aux
      M(i,j) = M(i,j)-aux
    end do
  end do
  if (No <= Nv) then
    ! Q = X M - Q P
    call dgemm_('N','N',Nv,No,No,One,Q,Nv,P,No,Zero,tmp,Nv)
    call dgemm_('N','N',Nv,No,No,One,X,Nv,M,No,Zero,Q,Nv)
  else
    ! Q = M X - P Q
    call dgemm_('N','N',Nv,No,Nv,One,P,Nv,Q,Nv,Zero,tmp,Nv)
    call dgemm_('N','N',Nv,No,Nv,One,M,Nv,X,Nv,Zero,Q,Nv)
  end if
  ! G = G + Q/(2k+1)!
  Q(:,:) = (Q(:,:)-tmp(:,:))/real(2*term*(2*term+1),kind=wp)
  G(:,:) = G(:,:)+Q(:,:)
  ! Series is converged
  if (maxval(abs(Q)) < thrs) exit
end do
! Not converged
if (term > maxterm) call Abend()

call mma_deallocate(P)
call mma_deallocate(M)
call mma_deallocate(Q)
call mma_deallocate(tmp)

end subroutine trg_series
