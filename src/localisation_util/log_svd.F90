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
!  Log_SVD
!
!> @brief Compute the antisymmetric matrix corresponding to a rotation matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the antisymmetric blocked matrix \f$ A \f$ such that \f$ \exp(A) \f$ produces an equivalent rotation
!> to the input rotation matrix (equivalent in the sense that the occupied and virtual orbitals will span the same space,
!> but the full rotation will not be identical). It also returns the redundant occupied-occupied and virtual-virtual
!> rotations to make them identical.
!>
!> @param[in]     N        Size of the square matrix
!> @param[in]     no       Size of the first diagonal zero block (number of occupied orbitals)
!> @param[in,out] A        Rotation matrix, it is replaced by the antisymmetric matrix with nonzero off-diagonal block,
!>                         and with the redundant rotations in the diagonal blocks
!***********************************************************************

subroutine Log_SVD(N,no,A)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, no
real(kind=wp), intent(inout) :: A(N,N)
integer(kind=iwp) :: i, nv
real(kind=wp), allocatable :: B(:,:), s(:), tmp(:), V(:,:), W(:,:)
real(kind=wp), parameter :: thrs = sqrt(epsilon(thrs))

! At the end we will get a matrix:
!
!   ( Uo -- )
!   ( X  Uv )
!
! such that the original A is equal to
!
!        ( 0  -X^T ) ( Uo 0  )
! A = exp( X    0  ) ( 0  Uv )

nv = N-no
if (no*nv < 1) return

! Pick up the occupied-occupied block
call mma_allocate(B,no,no,label='B')
B(:,:) = A(1:no,1:no)

call mma_allocate(V,no,no,label='V')
call mma_allocate(W,no,no,label='W')
call mma_allocate(s,no,label='s')

! Do the SVD, note the W matrix is already transposed
!   A(oo) = V s W^T
call large_svd(no,no,B,V,W,s)

! Redundant occupied-occupied rotation
!   U(o) = V W^T
call dgemm_('N','N',no,no,no,One,V,no,W,no,Zero,A(1,1),N)

! Redundant virtual-virtual rotation
!   U(v) = A(vv) - A(vo) W (1/(s+1)) V^T A(ov)
call mma_allocate(tmp,no*max(no,nv),label='tmp')
do i=1,no
  tmp((i-1)*no+1:i*no) = V(:,i)/(s(i)+One)
end do
! B = V (1/(s+1)) W^T
call dgemm_('N','N',no,no,no,One,tmp,no,W,no,Zero,B,no)
! U(v) = A(vv) - A(vo) B^T A(ov)
call dgemm_('T','N',no,nv,no,One,B,no,A(1,no+1),N,Zero,tmp,no)
call dgemm_('N','N',nv,nv,no,-One,A(no+1,1),N,tmp,no,One,A(no+1,no+1),N)

! pick up the virtual-occupied block
tmp(1:no*nv) = pack(A(no+1:,1:no),.true.)

! X = A(vo) (V scos(s) W^T)^T
! scos(x) = acos(x)/sqrt(1-x^2) ; scos(1) = 1
do i=1,no
  if (One-s(i) > thrs) V(:,i) = acos(s(i))/sqrt(One-s(i)**2)*V(:,i)
end do
call dgemm_('N','N',no,no,no,One,V,no,W,no,Zero,B,no)
call dgemm_('N','T',nv,no,no,One,tmp,nv,B,no,Zero,A(no+1,1),N)

call mma_deallocate(B)
call mma_deallocate(V)
call mma_deallocate(W)
call mma_deallocate(s)
call mma_deallocate(tmp)

end subroutine Log_SVD
