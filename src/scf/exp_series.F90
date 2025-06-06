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
! Copyright (C) 2022, Danjo De Chavez                                  *
!***********************************************************************
!  exp_series
!
!> @brief  Compute the exponential of the U matrix
!> @author Danjo De Chavez
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix U through
!> a modified Taylor expansion of the K matrix. This takes advantage of the
!> subblocks and trend in the powers of K matrix. \cite Sei2022-JCTC-18-4164
!>
!> \note
!> Some equations in the reference are wrong.
!>
!> @see ::exp_series2, ::exp_svd, ::exp_schur, ::exp_eig
!>
!> @param[in]     N        Size of the square matrix
!> @param[in]     No       Number of Occupied Orbitals
!> @param[in,out] U        U matrix is replaced by its exponential
!***********************************************************************

subroutine exp_series(N,No,U)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: N, No
real(kind=wp), intent(inout) :: U(N,N)
integer(kind=iwp) :: cnt, Nv
real(kind=wp) :: factor, ithrsh, ithrshoo, ithrshvo, ithrshvv
real(kind=wp), allocatable :: Koo(:,:), Kvo(:,:), Kvv(:,:), theta(:,:), Uoo(:,:), Uvo(:,:), Uvv(:,:), xUoo(:,:), xUvo(:,:), &
                              xUvv(:,:)
real(kind=wp), parameter :: thrsh = 1.0e-20_wp

if (N < 1) return
Nv = N-No

call mma_allocate(theta,Nv,No,label='theta')

call mma_allocate(Koo,No,No,label='Koo')
call mma_allocate(Kvv,Nv,Nv,label='Kvv')
call mma_allocate(Kvo,Nv,No,label='Kvo')

call mma_allocate(Uoo,No,No,label='Uoo')
call mma_allocate(Uvv,Nv,Nv,label='Uvv')
call mma_allocate(Uvo,Nv,No,label='Uvo')

call mma_allocate(xUoo,No,No,label='xUoo')
call mma_allocate(xUvv,Nv,Nv,label='xUvv')
call mma_allocate(xUvo,Nv,No,label='xUvo')

theta(:,:) = U(No+1:N,:No)

U(:,:) = Zero

cnt = 1
factor = One

ithrsh = 2.0e-16_wp

! Initialization
! Taylor expansion terms to n=1

Uvo(:,:) = Zero

call unitmat(Uoo,No)

xUoo(:,:) = Zero
xUvv(:,:) = Zero
xUvo(:,:) = Zero

call unitmat(Uvv,Nv)

Kvo(:,:) = theta
Uvo(:,:) = theta

! Main Loop
! Taylor expansion terms from n=2 to convergence

do while (thrsh < ithrsh)
  cnt = cnt+1
  factor = factor*cnt

  if (mod(cnt,2) == 0) then
    call dgemm_('T','N',No,No,Nv,One,-theta,Nv,Kvo,Nv,Zero,Koo,No)
    call dgemm_('N','T',Nv,Nv,No,One,Kvo,Nv,-theta,Nv,Zero,Kvv,Nv)
    Uoo(:,:) = Uoo+Koo/factor
    Uvv(:,:) = Uvv+Kvv/factor

  else
    call dgemm_('N','N',Nv,No,No,One,theta,Nv,Koo,No,Zero,Kvo,Nv)
    Uvo(:,:) = Uvo+Kvo/factor

    ithrshoo = maxval(abs(Uoo-xUoo)/(abs(Uoo)+thrsh))
    ithrshvv = maxval(abs(Uvv-xUvv)/(abs(Uvv)+thrsh))
    ithrshvo = maxval(abs(Uvo-xUvo)/(abs(Uvo)+thrsh))
    ithrsh = max(ithrshoo,ithrshvv,ithrshvo)

    xUoo(:,:) = Uoo
    xUvv(:,:) = Uvv
    xUvo(:,:) = Uvo

  end if
end do

U(:No,:No) = Uoo
U(No+1:N,:No) = Uvo
U(:No,No+1:N) = -transpose(Uvo)
U(No+1:N,No+1:N) = Uvv

call mma_deallocate(Koo)
call mma_deallocate(Kvv)
call mma_deallocate(Kvo)

call mma_deallocate(Uoo)
call mma_deallocate(Uvv)
call mma_deallocate(Uvo)

call mma_deallocate(xUoo)
call mma_deallocate(xUvv)
call mma_deallocate(xUvo)

call mma_deallocate(theta)

return

end subroutine exp_series
