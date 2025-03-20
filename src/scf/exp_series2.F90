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
!  exp_series2
!
!> @brief  Compute the exponential of the U matrix
!> @author Ignacio Fdez. Galvan
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix U through
!> a Taylor expansion of the K matrix.
!>
!> @see ::exp_series, ::exp_svd, ::exp_schur, ::exp_eig
!>
!> @param[in]     N        Size of the square matrix
!> @param[in]     No       Number of Occupied Orbitals
!> @param[in,out] U        U matrix is replaced by its exponential
!***********************************************************************

subroutine exp_series2(N,No,U)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: N, No
real(kind=wp), intent(inout) :: U(N,N)
integer(kind=iwp) :: i, k, Nv
real(kind=wp) :: f1, f2, f3, mval
real(kind=wp), allocatable :: cosP(:,:), mcosP(:,:), sincP(:,:), P(:,:), tmp(:,:), Uaa(:,:), Uab(:,:), Ubb(:,:)
integer(kind=iwp), parameter :: maxorder=50
real(kind=wp), parameter :: thrs = 1.0e-16_wp

if (N < 1) return
Nv = N-No
if (No*Nv < 1) then
  call unitmat(U,N)
  return
end if

! initially:
!     ( 0  0 )        ( 0  -X^T )
! U = ( X  0 ) -> A = ( X    0  )
!
! P = X^T X
!
!          (  cos(P^(1/2))         -sinc(P^(1/2))     )
! exp(A) = ( X sinc(P^(1/2))  I - X mcos(P^(1/2)) X^T )
!
! mcos(x) = (1-cos(x))/x^2

k = min(No,Nv)
call mma_allocate(P,k,k,label='P')
call mma_allocate(cosP,k,k,label='cosP')
call mma_allocate(mcosP,k,k,label='mcosP')
call mma_allocate(sincP,k,k,label='sincP')
call mma_allocate(Uaa,k,k,label='Uaa')
call mma_allocate(Uab,k,k,label='Uab')
call mma_allocate(Ubb,k,k,label='Ubb')
call mma_allocate(tmp,k,k,label='tmp')

! Most of the time No < Nv, but just in case, we select the smaller matrices to multiply
! We include the minus sign here to take care of the alternating sign in the Taylor expansions
if (No <= Nv) then
  ! P = -X^T X
  call dgemm_('T','N',No,No,Nv,-One,U(No+1,1),N,U(No+1,1),N,Zero,P,No)
else
  ! P = -X X^T
  call dgemm_('N','T',Nv,Nv,No,-One,U(No+1,1),N,U(No+1,1),N,Zero,P,Nv)
end if
call unitmat(cosP,k)  ! cos(P^(1/2))
call unitmat(sincP,k) ! sinc(P^(1/2)) = sin(P^(1/2))*P^(-1/2)
mcosP(:,:) = Zero     ! (I-cos(P^(1/2)))*P^(-1)
call unitmat(Uaa,k)
call unitmat(Uab,k)
call unitmat(Ubb,k)
Ubb(:,:) = Ubb(:,:)*Half
do i=2,maxorder,2
  f1 = One/real((i-1)*i,kind=wp)
  call dgemm_('N','N',k,k,k,f1,Uaa,k,P,k,Zero,tmp,k)
  Uaa(:,:) = tmp(:,:)
  cosP(:,:) = cosP(:,:)+Uaa(:,:)
  f2 = One/real(i*(i+1),kind=wp)
  call dgemm_('N','N',k,k,k,f2,Uab,k,P,k,Zero,tmp,k)
  Uab(:,:) = tmp(:,:)
  sincP(:,:) = sincP(:,:)+Uab(:,:)
  mval = maxval(abs(Uaa(:,:)))
  mcosP(:,:) = mcosP(:,:)+Ubb(:,:)
  if (mval < thrs) exit
  ! This is done after checking for convergence, to ensure the order of expansion is consistent
  f3 = One/real((i+1)*(i+2),kind=wp)
  call dgemm_('N','N',k,k,k,f3,Ubb,k,P,k,Zero,tmp,k)
  Ubb(:,:) = tmp(:,:)
end do
call mma_deallocate(Uaa)
call mma_deallocate(Uab)
call mma_deallocate(Ubb)
call mma_deallocate(tmp)
if (No <= Nv) then
  call mma_allocate(tmp,Nv,No,label='tmp')
  ! OO block   cos(P^(1/2))
  U(1:No,1:No) = cosP(:,:)
  ! VV block   I - X (I-cos(P^(1/2)))*P^(-1) X^T
  call dgemm_('N','N',Nv,No,No,-One,U(No+1,1),N,mcosP,No,Zero,tmp,Nv)
  call dgemm_('N','T',Nv,Nv,No,One,tmp,Nv,U(No+1,1),N,Zero,U(No+1,No+1),N)
  do i=1,Nv
    U(No+i,No+i) = U(No+i,No+i)+One
  end do
  ! VO block   X sinc(P^(1/2))
  call dgemm_('N','N',Nv,No,No,One,U(No+1,1),N,sincP,No,Zero,tmp,Nv)
  U(No+1:,1:No) = tmp(:,:)
  ! OV block   -VO^T
  U(1:No,No+1:) = -transpose(tmp(:,:))
else
  call mma_allocate(tmp,No,Nv,label='tmp')
  ! VV block   cos(P^(1/2))
  U(No+1:,No+1:) = cosP(:,:)
  ! OO block   I - X (I-cos(P^(1/2)))*P^(-1) X^T
  call dgemm_('T','N',No,Nv,Nv,-One,U(No+1,1),N,mcosP,Nv,Zero,tmp,No)
  call dgemm_('N','N',No,No,Nv,One,tmp,No,U(No+1,1),N,Zero,U(1,1),N)
  do i=1,No
    U(i,i) = U(i,i)+One
  end do
  ! OV block   -X^T sinc(P^(1/2))
  call dgemm_('T','N',No,Nv,Nv,-One,U(No+1,1),N,sincP,Nv,Zero,tmp,No)
  U(1:No,No+1:) = tmp(:,:)
  ! VO block   -OV^T
  U(No+1:,1:No) = -transpose(tmp(:,:))
end if

call mma_deallocate(P)
call mma_deallocate(cosP)
call mma_deallocate(sincP)
call mma_deallocate(mcosP)
call mma_deallocate(tmp)

end subroutine exp_series2
