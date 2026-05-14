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
!  Exp_SVD
!
!> @brief Compute the exponential of a single-block antisymmetric matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix \f$ A \f$ through its singular value decomposition.
!> The exponential of \f$ A \f$ is an orthogonal matrix.
!> This is specialized for the case where there is a single nonzero off-diagonal block in \f$ A \f$
!>
!> @see ::exp_series, ::exp_series2, ::exp_schur, ::exp_eig
!>
!> @param[in]     N        Size of the square matrix
!> @param[in]     no       Size of the first diagonal zero block (number of occupied orbitals)
!> @param[in,out] A        Antisymmetric real matrix, it is replaced by its exponential
!> @param[out]    maxtheta Maximum rotation angle (singular values)
!***********************************************************************

subroutine Exp_SVD(N,no,A,maxtheta)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, no
real(kind=wp), intent(inout) :: A(N,N)
real(kind=wp), intent(out) :: maxtheta
integer(kind=iwp) :: i, j, k, nv
real(kind=wp), allocatable :: Q(:,:), R(:,:), s(:), sx(:), tmp(:), X(:,:)

maxtheta = Zero
if (N < 1) return
nv = N-no
if (no*nv == 0) then
  call unitmat(A,N)
  return
end if

! Pick up the nonzero block
!     ( 0  -X^T )
! A = ( X    0  )
call mma_allocate(X,nv,no,label='X')
X(:,:) = A(no+1:,1:no)

k = min(no,nv)
call mma_allocate(Q,nv,k,label='Q')
call mma_allocate(R,k,no,label='R')
call mma_allocate(s,k,label='s')

! Do the SVD, note the R matrix is already transposed
call large_svd(nv,no,X,Q,R,s)
call mma_deallocate(X)
maxtheta = s(1)

call mma_allocate(sx,k,label='sx')
call mma_allocate(tmp,no*nv,label='tmp')

! Build the rotation matrix
! Based on J. Chem. Phys. 101 (1994) 3862 doi:10.1063/1.467504
! but using SVD of X (or, actually, of -X^T in the paper) instead of EVD of X X^T

sx(:) = cos(s(:))
! Choose the option with smaller matrix multiplications
if (no <= nv) then
  ! OO block: R cos(s) R^T
  j = 0
  do i=1,no
    tmp(j+1:j+k) = sx(:)*R(:,i)
    j = j+k
  end do
  call dgemm_('T','N',no,no,k,One,tmp,k,R,k,Zero,A(1,1),N)
  ! VV block: I + Q (cos(s)-I) Q^T
  sx(:) = sx(:)-One
  j = 0
  do i=1,k
    tmp(j+1:j+nv) = sx(i)*Q(:,i)
    j = j+nv
  end do
  call dgemm_('N','T',nv,nv,k,One,tmp,nv,Q,nv,Zero,A(no+1,no+1),N)
  do i=1,nv
    A(no+i,no+i) = A(no+i,no+i)+One
  end do
else
  ! VV block: Q cos(s) Q^T
  j = 0
  do i=1,k
    tmp(j+1:j+nv) = sx(i)*Q(:,i)
    j = j+nv
  end do
  call dgemm_('N','T',nv,nv,k,One,tmp,nv,Q,k,Zero,A(no+1,no+1),N)
  ! OO block: I + R (cos(s)-I) R^T
  sx(:) = sx(:)-One
  j = 0
  do i=1,no
    tmp(j+1:j+k) = sx(:)*R(:,i)
    j = j+k
  end do
  call dgemm_('T','N',no,no,k,One,tmp,k,R,k,Zero,A(1,1),N)
  do i=1,no
    A(i,i) = A(i,i)+One
  end do
end if
! VO block: Q sin(s) R^T
j = 0
do i=1,k
  tmp(j+1:j+nv) = sin(s(i))*Q(:,i)
  j = j+nv
end do
call dgemm_('N','N',nv,no,k,One,tmp,nv,R,k,Zero,A(no+1,1),N)
! OV block: -VO^T
A(1:no,no+1:) = -transpose(A(no+1:,1:no))

call mma_deallocate(R)
call mma_deallocate(Q)
call mma_deallocate(s)
call mma_deallocate(sx)
call mma_deallocate(tmp)

end subroutine Exp_SVD
