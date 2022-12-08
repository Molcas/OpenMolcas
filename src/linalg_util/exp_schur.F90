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
! Copyright (C) 2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Exp_Schur
!
!> @brief Compute the exponential of an antisymmetric matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix \f$ X \f$ through its Schur decomposition.
!> The exponential of \f$ X \f$ is an orthogonal matrix.
!> The Schur form of \f$ X \f$ is antisymmetric real, block-diagonal, with \f$ 2\times 2\f$ or \f$ 1\times 1\f$ diagonal blocks.
!> If \f$ X = Z T Z^T \f$, then \f$ \exp(X) = Z \exp(T) Z^T \f$, and \f$ exp(T) \f$ is trivial to compute since each block can
!> be treated separately: \f$ (0, \lambda) \to (\cos(\lambda), \pm\sin(\lambda) \f$.
!>
!> @param[in]     N        Size of the square matrix
!> @param[in,out] X        Antisymmetric real matrix, it is replaced by its exponential
!> @param[out]    maxtheta Maximum rotation angle (absolute value in the Schur form)
!***********************************************************************

subroutine Exp_Schur(N,X,maxtheta)

!define _USE_LAPACK_
#ifdef _USE_LAPACK_
use sorting_funcs, only: geq_r
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: X(N,N)
integer(kind=iwp) :: i, info
real(kind=wp) :: c, maxtheta, s, theta
real(kind=wp), allocatable :: tmp(:,:), vs(:,:), work(:)
#ifdef _USE_LAPACK_
integer(kind=iwp) :: j, lwork, sdim
real(kind=wp) :: dwrk(1)
logical(kind=iwp) :: bwork(1), pair
real(kind=wp), allocatable :: wri(:,:)
#endif
real(kind=wp), parameter :: thrs = epsilon(thrs)

maxtheta = Zero
if (N < 1) return

#ifdef _USE_LAPACK_

! Ensure the matrix is antisymmetric
do i=1,N
  do j=i+1,N
    X(i,j) = -X(j,i)
  end do
end do

call mma_allocate(wri,N,2,label='wri')
call mma_allocate(vs,N,N,label='vs')

! Schur factorization of X: X = Z T Z^T, exp(X) = Z exp(T) Z^T
! (this, unfortunately, does not take advantage of the antisymmetry)
lwork = -1
call dgees_('V','N',geq_r,N,X,N,sdim,wri(:,1),wri(:,2),vs,N,dwrk,lwork,bwork,info)
lwork = int(dwrk(1),kind=iwp)
call mma_allocate(work,lwork,label='work')
call dgees_('V','N',geq_r,N,X,N,sdim,wri(:,1),wri(:,2),vs,N,work,lwork,bwork,info)
call mma_deallocate(work)
call mma_deallocate(wri)

if (info /= 0) call abend()

call mma_allocate(tmp,N,N,label='tmp')

! Compute tmp = Z exp(T)
!
! tmp will contain diagonal blocks that are either 0, and exp(0) = 1
! or of the form  C = [ 0 x], and exp(C) = [ cos(x) sin(x)]
!                     [-x 0]               [-sin(x) cos(x)]
i = 1
do while (i <= N)
  if (i == N) then
    pair = .false.
  else
    pair = abs(X(i+1,i)) > thrs
  end if
  if (pair) then
    theta = sign(sqrt(-X(i+1,i)*X(i,i+1)),X(i,i+1))
    maxtheta = max(maxtheta, abs(theta))
    c = cos(theta)
    s = sin(theta)
    tmp(:,i) = c*vs(:,i)-s*vs(:,i+1)
    tmp(:,i+1) = s*vs(:,i)+c*vs(:,i+1)
    i = i+2
  else
    tmp(:,i) = vs(:,i)
    i = i+1
  end if
end do

! Compute exp(X) = tmp Z^T
call dgemm_('N','T',N,N,N,One,tmp,N,vs,N,Zero,X,N)

call mma_deallocate(tmp)
call mma_deallocate(vs)

#else

call mma_allocate(work,N,label='work')
call mma_allocate(vs,N,N,label='vs')

! Since dgemm does not work in-place, and we want
! the result in-place, there must be a copy at some point
vs(:,:) = X
call schur_skew(N,vs,work,info)
if (info /= 0) call abend()

call mma_allocate(tmp,N,N,label='tmp')

i = 1
do while (i <= N)
  theta = work(i)
  if ((i < N) .and. (abs(theta) > thrs)) then
    maxtheta = max(maxtheta, abs(theta))
    c = cos(theta)
    s = sin(theta)
    tmp(:,i) = c*vs(:,i)-s*vs(:,i+1)
    tmp(:,i+1) = s*vs(:,i)+c*vs(:,i+1)
    i = i+2
  else
    tmp(:,i) = vs(:,i)
    i = i+1
  end if
end do

call mma_deallocate(work)

! Compute exp(X) = tmp Z^T
call dgemm_('N','T',N,N,N,One,tmp,N,vs,N,Zero,X,N)

call mma_deallocate(vs)
call mma_deallocate(tmp)

#endif

end subroutine Exp_Schur
