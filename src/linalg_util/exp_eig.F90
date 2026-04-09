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
!  Exp_eig
!
!> @brief Compute the exponential of an antisymmetric matrix
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix \f$ X \f$ through the eigenvalue decomposition of \f$ P = X X^T \f$
!> The exponential of \f$ X \f$ is an orthogonal matrix.
!> The product \f$ P = -X^2 \f$ is symmetric and positive (semi-)definite. Its eigenvalues are the squares of the eigenvalues
!> of \f$ X \f$. The exponential series can be split in two subseries for the odd and even elements, which correspond to
!> \f$ \cos(\sqrt(P)) \f$ and \f$ P \sin(\sqrt(P))/P \f$.
!>
!> @see ::exp_schur, ::exp_series, ::exp_series2, ::exp_svd
!>
!> @param[in]     N        Size of the square matrix
!> @param[in,out] X        Antisymmetric real matrix, it is replaced by its exponential
!> @param[out]    maxtheta Maximum rotation angle
!***********************************************************************

subroutine Exp_eig(N,X,maxtheta)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: X(N,N)
real(kind=wp), intent(out) :: maxtheta
integer(kind=iwp) :: i, j, k
real(kind=wp), allocatable :: P(:,:), Ptri(:), tmp(:,:)
real(kind=wp), parameter :: thrs = sqrt(epsilon(thrs))

maxtheta = Zero
if (N < 1) return

call mma_allocate(P,N,N,label='P')
call mma_allocate(Ptri,nTri_Elem(N),label='Ptri')

! Ensure the matrix is antisymmetric
do i=1,N
  do j=i+1,N
    X(i,j) = -X(j,i)
  end do
end do
! P = X X^T = -X X
call dgemm_('N','T',N,N,N,One,X,N,X,N,Zero,P,N)
k = 0
do j=1,N
  do i=1,j
    k = k+1
    Ptri(k) = P(i,j)
  end do
end do

! Eigenvalue decomposition of P
!   P = U s U^T
! Eigenvectors returned in P
call NIdiag_New(Ptri,P,N,N)
! Put sqrt(s) in the first N elements of Ptri
do i=1,N
  Ptri(i) = sqrt(max(Zero,Ptri(nTri_Elem(i))))
  maxtheta = max(maxtheta,Ptri(i))
end do

! exp(X) = [X U sinc(sqrt(s)) + U cos(sqrt(s))] U^T
call mma_allocate(tmp,N,N,label='tmp')
call dgemm_('N','N',N,N,N,One,X,N,P,N,Zero,tmp,N)
do i=1,N
  if (Ptri(i) < thrs) then
    tmp(:,i) = tmp(:,i)+P(:,i)
  else
    tmp(:,i) = sin(Ptri(i))/Ptri(i)*tmp(:,i)+cos(Ptri(i))*P(:,i)
  end if
end do
call dgemm_('N','T',N,N,N,One,tmp,N,P,N,Zero,X,N)

call mma_deallocate(P)
call mma_deallocate(Ptri)
call mma_deallocate(tmp)

end subroutine Exp_eig
