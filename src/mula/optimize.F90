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
!    var       : Real*8 two dimensional array - coordinates
!                for which we know the energy.
!    coef      : Real*8 two dimensional array - coefficients
!                for each term specified by ipow.
!
!  Output:
!    x         : Real*8 array - geometry in minimum.
!    energy    : Real*8 - energy in minimum.
!    Hess      : Real*8 two dimensional array - Hessian in
!                minimum.
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

use Constants, only: Zero, One, Two
use Definitions, only: wp, u6

!use Linalg
implicit real*8(a-h,o-z)
parameter(maxiter=100)
parameter(delta_max=One)
integer ipow(nterm,nvar)
real*8 var(ndata,nvar)
real*8 coef(nterm,1)
real*8 x(nvar)
real*8 xmin(nvar)
real*8 grad(nvar)
real*8 delta(nvar,1)
real*8 Hess(nterm,nvar)
real*8 var_intervals(nvar,2)
logical shift

! Initialize.

! Determine intervals for random starting point in optimization.
do iv=1,nvar
  var_intervals(iv,1) = var(1,iv)
  var_intervals(iv,2) = var(1,iv)
end do
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
  do iv=1,nvar
    delta(iv,1) = -grad(iv)
  end do
  call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
  !call calcNorm(delta(:,1),delta_norm)
  sum = Zero
  do iv=1,nvar
    sum = sum+delta(iv,1)**2
  end do
  delta_norm = sqrt(sum)

  scale = One
  if (delta_norm > delta_max) scale = delta_max/delta_norm
  do iv=1,nvar
    x(iv) = x(iv)+delta(iv,1)*scale
  end do
  iter = 0
  do while ((delta_norm > 1.0e-12_wp) .and. (iter <= maxiter))
    iter = iter+1
    call gradient(x,coef,ipow,grad,nterm,nvar)
    call Hessian(x,coef,ipow,Hess,nterm,nvar)
    call ShiftHess(Hess,shift,nterm,nvar)
    do iv=1,nvar
      delta(iv,1) = -grad(iv)
    end do
    call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
    !call calcNorm(delta(:,1), delta_norm)
    sum = Zero
    do iv=1,nvar
      sum = sum+delta(iv,1)**2
    end do
    delta_norm = sqrt(sum)

    scale = One
    if (delta_norm > delta_max) scale = delta_max/delta_norm
    do iv=1,nvar
      x(iv) = x(iv)+delta(iv,1)*scale
    end do
  end do
  if (iter >= maxiter) write(u6,*) 'WARNING!! No convergence in Optimize'
  call funcval(x,coef,ipow,fval,nterm,nvar)
  if (fval < energy) then
    energy = fval
    do iv=1,nvar
      xmin(iv) = x(iv)
    end do
  end if
end do
do iv=1,nvar
  x(iv) = xmin(iv)
end do

end subroutine Optimize
