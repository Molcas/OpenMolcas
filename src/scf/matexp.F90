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
!  matexp
!
!> @brief  Compute the exponential of the U matrix
!> @author Danjo De Chavez (2022)
!>
!> @details
!> Computes the exponential of an antisymmetric real matrix U through
!> a modified Taylor expansion of the K matrix. This takes advantage of the
!> subblocks and trend in the powers of K matrix. \cite Sei2022-JCTC-18-4164
!>
!> \note
!> Some equations in the reference are wrong.
!>
!> @param[in]     N        Size of the square matrix
!> @param[in]     No       Number of Occupied Orbitals
!> @param[in,out] U        U matrix is replaced by its exponential
!***********************************************************************

subroutine matexp(N,No,U)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: N, No
real(kind=wp), intent(inout)  :: U(N,N)

integer(kind=iwp)  :: Nv
integer(kind=iwp)  :: count, i

real(kind=wp), parameter :: thrsh = 1.0D-20
real(kind=wp) :: ithrsh, ithrshoo, ithrshvv, ithrshvo, factor

real(kind=wp), allocatable :: Uoo(:,:), xUoo(:,:), Koo(:,:)
real(kind=wp), allocatable :: Uvv(:,:), xUvv(:,:), Kvv(:,:)
real(kind=wp), allocatable :: Uvo(:,:), xUvo(:,:), Kvo(:,:)
real(kind=wp), allocatable :: Uov(:,:)
real(kind=wp), allocatable :: theta(:,:)

if (N < 1) return
Nv = N-No

call mma_allocate(theta,Nv,No,label='theta')

call mma_allocate(Koo,No,No,label='Koo')
call mma_allocate(Kvv,Nv,Nv,label='Kvv')
call mma_allocate(Kvo,Nv,No,label='Kvo')

call mma_allocate(Uoo,No,No,label='Uoo')
call mma_allocate(Uvv,Nv,Nv,label='Uvv')
call mma_allocate(Uov,No,Nv,label='Uov')
call mma_allocate(Uvo,Nv,No,label='Uvo')

call mma_allocate(xUoo,No,No,label='xUoo')
call mma_allocate(xUvv,Nv,Nv,label='xUvv')
call mma_allocate(xUvo,Nv,No,label='xUvo')

theta(:,:) = U(No+1:N,:No)

U(:,:) = Zero

count  = 1
factor = One

ithrsh = 2.0E-16_wp

! Initialization
! Taylor expansion terms to n=1

Uov(:,:) = Zero
Uvo(:,:) = Zero

Uoo(:,:) = Zero

xUoo(:,:) = Zero
xUvv(:,:) = Zero
xUvo(:,:) = Zero

do i=1,No
  Uoo(i,i) = One
end do

Uvv(:,:) = Zero

do i=1,Nv
  Uvv(i,i) = One
end do

Kvo(:,:) = theta
Uvo(:,:) = theta

! Main Loop
! Taylor expansion terms from n=2 to convergence

do while (thrsh < ithrsh)
  count = count+1

  if (mod(count,2)==0) then
    call dgemm_('T','N',No,No,Nv,One,-theta,Nv,Kvo,Nv,Zero,Koo,No)
    call dgemm_('N','T',Nv,Nv,No,One,Kvo,Nv,-theta,Nv,Zero,Kvv,Nv)
    factor = factor*count
    Uoo(:,:) = Uoo + Koo/factor
    Uvv(:,:) = Uvv + Kvv/factor

  else
    call dgemm_('N','N',Nv,No,No,One,theta,Nv,Koo,No,Zero,Kvo,Nv)
    factor = factor*count
    Uvo(:,:) = Uvo + Kvo/factor

    ithrshoo = maxval(abs(Uoo-xUoo)/(abs(Uoo)+thrsh))
    ithrshvv = maxval(abs(Uvv-xUvv)/(abs(Uvv)+thrsh))
    ithrshvo = maxval(abs(Uvo-xUvo)/(abs(Uvo)+thrsh))
    ithrsh   = max(ithrshoo, ithrshvv, ithrshvo)

    xUoo(:,:) = Uoo
    xUvv(:,:) = Uvv
    xUvo(:,:) = Uvo

  end if
end do

Uov(:,:) = -transpose(Uvo)

U(:No,:No)       = Uoo
U(No+1:N,:No)    = Uvo
U(:No,No+1:N)    = Uov
U(No+1:N,No+1:N) = Uvv

call mma_deallocate(Koo)
call mma_deallocate(Kvv)
call mma_deallocate(Kvo)

call mma_deallocate(Uoo)
call mma_deallocate(Uvv)
call mma_deallocate(Uov)
call mma_deallocate(Uvo)

call mma_deallocate(xUoo)
call mma_deallocate(xUvv)
call mma_deallocate(xUvo)

call mma_deallocate(theta)

return

end subroutine matexp
