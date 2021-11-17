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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine Optimize(ipow,var,coef,x,energy,Hess,nterm,nvar,ndata)
!  Purpose:
!    Optimize structure given the polynomial fit from PolFit.
!
!  Input:
!    ipow      : Two dimensional integer array - terms in polynomial.
!    var       : Real two dimensional array - coordinates for which we know the energy.
!    coef      : Real two dimensional array - coefficients for each term specified by ipow.
!
!  Output:
!    x         : Real array - geometry in minimum.
!    energy    : Real - energy in minimum.
!    Hess      : Real two dimensional array - Hessian in minimum.
!
!  Calls:
!    Dool_MULA
!    funcval
!    gradient
!    Hessian
!    thirdDer
!    fourthDer
!    ShiftHess
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar), ndata
real(kind=wp), intent(in) :: var(ndata,nvar), coef(nterm,1)
real(kind=wp), intent(out) :: x(nvar), energy, Hess(nterm,nvar)
integer(kind=iwp) :: i, iseed, iter, iterm, iv, ivar, j
real(kind=wp) :: delta_norm, det, fval, rand_number, rscale, rsum
logical(kind=iwp) :: shift
real(kind=wp), allocatable :: delta(:,:), grad(:), var_intervals(:,:), xmin(:)
integer(kind=iwp), parameter :: maxiter = 100
real(kind=wp), parameter :: delta_max = One
real(kind=r8), external :: Random_Molcas

! Initialize.

call mma_allocate(var_intervals,nvar,2,label='var_intervals')
call mma_allocate(grad,nvar,label='grad')
call mma_allocate(delta,nvar,1,label='delta')
call mma_allocate(xmin,nvar,label='xmin')

! Determine intervals for random starting point in optimization.
var_intervals(:,1) = var(1,:)
var_intervals(:,2) = var(1,:)
do ivar=1,nvar
  do iterm=2,nterm
    if (var(iterm,ivar) < var_intervals(ivar,1)) then
      var_intervals(ivar,1) = var(iterm,ivar)
    else if (var(iterm,ivar) > var_intervals(ivar,2)) then
      var_intervals(ivar,2) = var(iterm,ivar)
    end if
  end do
end do

! Optimize using Newton-Raphson.
iseed = 12345
energy = 1.0e3_wp
!call Srand(Half)
do i=1,1
  do j=1,nvar
    !rand_number = rand()
    ! Random_molcas takes an integer seed 0 < iseed < 2**31
    ! and produces a pseudorandom number 0 < z < 1.0
    rand_number = Random_molcas(iseed)
    iseed = nint(rand_number*Two**31)
    x(j) = var_intervals(j,1)+(var_intervals(j,2)-var_intervals(j,1))*rand_number
  end do
  call gradient(x,coef,ipow,grad,nterm,nvar)
  call Hessian(x,coef,ipow,Hess,nterm,nvar)
  call ShiftHess(Hess,shift,nterm,nvar)
  delta(:,1) = -grad
  call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
  !call calcNorm(delta(:,1),delta_norm)
  rsum = Zero
  do iv=1,nvar
    rsum = rsum+delta(iv,1)**2
  end do
  delta_norm = sqrt(rsum)

  rscale = One
  if (delta_norm > delta_max) rscale = delta_max/delta_norm
  x(:) = x(:)+delta(:,1)*rscale
  iter = 0
  do while ((delta_norm > 1.0e-12_wp) .and. (iter <= maxiter))
    iter = iter+1
    call gradient(x,coef,ipow,grad,nterm,nvar)
    call Hessian(x,coef,ipow,Hess,nterm,nvar)
    call ShiftHess(Hess,shift,nterm,nvar)
    delta(:,1) = -grad
    call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
    !call calcNorm(delta(:,1), delta_norm)
    rsum = Zero
    do iv=1,nvar
      rsum = rsum+delta(iv,1)**2
    end do
    delta_norm = sqrt(rsum)

    rscale = One
    if (delta_norm > delta_max) rscale = delta_max/delta_norm
    x(:) = x(:)+delta(:,1)*rscale
  end do
  if (iter >= maxiter) write(u6,*) 'WARNING!! No convergence in Optimize'
  call funcval(x,coef,ipow,fval,nterm,nvar)
  if (fval < energy) then
    energy = fval
    xmin(:) = x
  end if
end do
x(:) = xmin

call mma_deallocate(var_intervals)
call mma_deallocate(grad)
call mma_deallocate(delta)
call mma_deallocate(xmin)

end subroutine Optimize
