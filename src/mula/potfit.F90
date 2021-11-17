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

!module OptMod

!  Contains:
!    PotFit       (ipow,var,yin,coef,x,energy,Hess,D3,D4,trfName,
!                  stand_dev,max_err,find_minimum)
!    Optimize     (ipow,var,coef,x,energy,Hess,D3,D4)
!    var_to_qvar  (var,qvar,ref,qref,alpha,trfName)
!    qvar_to_var  (var,x,grad,Hess,D3,D4,ref,qref,trfName,alpha)
!    funcval      (x,coef,ipow) Result(sum)
!    gradient     (x,coef,ipow) Result(grad)
!    Hessian      (x,coef,ipow) Result(Hess)
!    thirdDer     (x,coef,ipow) Result(D3)
!    fourthDer    (x,coef,ipow) Result(D4)
!    ShiftHess    (Hess,shift)
!
!  Uses:
!    Constants
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!contains

subroutine PotFit(nterm,nvar,ndata,ipow,var,yin,coef,x,nOsc,energy,grad,Hess,D3,D4,trfName,stand_dev,max_err,find_minimum, &
                  max_term,use_weight,l_Hess_1,l_Hess_2,l_D)
!  Purpose:
!    Fit polynomial to data and then depending upun the variable
!    find_minimum, we either calculate the energy, gradient, Hessian,
!    and qubic and quartic force constants in x or we optimize the
!    structure before we perform the calculations.
!
!  Input:
!    ipow         : Two dimensional integer array - terms in polynomial.
!    var          : Real two dimensional array - coordinates for which we know the energy.
!    yin          : Real array - energy in points given in var.
!    trfName      : Character array - transformation associated with each internal coordinate.
!    find_minimum : Logical variable - whether or not to optimize structure.
!    max_term     : Integer - highest power if term in polynomial fit.
!    use_weight   : Logical variable - whether or not to use weights to make the minimum more important in the fit.
!
!  Output:
!    coef         : Real two dimensional array - coefficients for each term specified by ipow.
!    x            : Real array - geometry in internal coordinates.
!    energy       : Real - energy in minimum.
!    Hess         : Real two dimensional array - Hessian in minimum.
!    D3           : Real three dimensional array - cubic force constants.
!    D4           : Real four dimensional array - quartic force constants.
!    stand_dev    : Real variable - standard deviation of fitted values in Subroutine PolFit.
!    max_err      : Real variable - maximum error of fitted values in Subroutine PolFit.
!
!  Uses:
!    LinAlg
!    VibMod
!    RandomMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use RandomMod, only: Ranmar, Rmarin
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, iPow(nterm,nvar), ndata, nOsc, max_term, l_Hess_1, l_Hess_2, l_D
real(kind=wp), intent(in) :: var(ndata,nvar), yin(ndata)
real(kind=wp), intent(out) :: coef(nterm,1), x(nosc), energy, grad(nosc), Hess(l_Hess_1,l_Hess_2), D3(l_D,l_D,l_D), &
                              D4(l_D,l_D,l_D,l_D), stand_dev, max_err
character(len=80), intent(in) :: trfName(nvar)
logical(kind=iwp), intent(in) :: find_minimum, use_weight
integer(kind=iwp) :: ifit, iPoints, iRandNum, istart, istop, ivar, ivarMax, length, nPoints
real(kind=wp) :: alphaMax, stand_dev_best
logical(kind=iwp) :: fit_alpha
character(len=32) :: Inline, trfcode
real(kind=wp), allocatable :: alpha(:), alpha0(:), alphaRad(:), alphaStart(:), alphaStop(:), alphaTemp(:), diff_vec(:), qref(:), &
                              qvar(:,:), RandVec(:), ref(:)

! Initialize.
coef(:,:) = Zero
call mma_allocate(qvar,ndata,nvar,label='qvar')
call mma_allocate(ref,nvar,label='ref')
call mma_allocate(qref,nvar,label='qref')

! Check if a fit of alpha is wanted.
fit_alpha = .false.
call mma_allocate(alphaStart,nvar,label='alphaStart')
call mma_allocate(alphaStop,nvar,label='alphaStop')
call mma_allocate(alpha,nvar,label='alpha')
call mma_allocate(alpha0,nvar,label='alpha')
call mma_allocate(alphaTemp,nvar,label='alpha')
call mma_allocate(alphaRad,nvar,label='alpha')
alphaStart(:) = -One
alphaStop(:) = -One
do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ifit = index(trfcode,'FIT')
  if (ifit > 0) then
    fit_alpha = .true.
    Inline = trfName(ivar)(1:32)
    istart = index(Inline,'ALPHA=')
    istart = istart+6
    Inline = trfCode(istart:len(trfCode))
    istop = len(Inline)
    read(Inline(1:istop),*) alphaStart(ivar),alphaStop(ivar)
    alphaRad(ivar) = (alphaStop(ivar)-alphaStart(ivar))*Half
    alpha0(ivar) = alphaStart(ivar)+alphaRad(ivar)
  end if
end do

call mma_allocate(diff_vec,ndata,label='diff_vec')

if (fit_alpha) then

  ! Determine which alpha radius is the largest.
  alphaMax = alphaRad(1)
  ivarMax = 1
  do ivar=2,nvar
    if (alphaRad(ivar) > alphaMax) then
      alphaMax = alphaRad(ivar)
      ivarMax = ivar
    end if
  end do

  ! Initialize random number calculation.
  nPoints = 20000
  call RmarIn(1802,9373)

  ! Start with the (nRandom+1)'th set of random numbers.
  !nRandom = 4
  !do i=1,nRandom
  !  call Ranmar(temp,2)
  !end do

  ! Store nPoints random numbers in vector RandVec.
  length = nvar*nPoints
  call mma_allocate(RandVec,length,label='RandVec')

  call Ranmar(RandVec,length)

  ! Simulated annealing.
  stand_dev_best = 100.0_wp
  iRandNum = 1
  do while (alphaRad(ivarMax) > 0.1_wp)
    alphaTemp(:) = alpha0
    do iPoints=1,1000
      do ivar=1,nvar
        if (alphaStart(ivar) > Zero) then
          alpha(ivar) = alpha0(ivar)+alphaRad(ivar)*(Two*RandVec(iRandNum+ivar-1)-One)
        end if
      end do
      iRandNum = iRandNum+nvar
      call var_to_qvar(var,qvar,ref,qref,alpha,trfName,ndata,nvar)
      call PolFit(ipow,nvar,qvar,yin,ndata,coef,nterm,stand_dev,max_err,diff_vec,use_weight)
      if (stand_dev < stand_dev_best) then
        stand_dev_best = stand_dev
        alphaTemp(:) = alpha
      end if
    end do
    alpha0(:) = alphaTemp
    write(u6,*) alpha0(:),stand_dev_best
    alphaRad(:) = alphaRad*Half
  end do
  alpha(:) = alpha0
  write(u6,*) alpha(1),alpha(2)
  call mma_deallocate(RandVec)

end if

call mma_deallocate(alphaStart)
call mma_deallocate(alphaStop)
call mma_deallocate(alpha0)
call mma_deallocate(alphaTemp)
call mma_deallocate(alphaRad)

! Transform coordinates to get better numerical fit.
call var_to_qvar(var,qvar,ref,qref,alpha,trfName,ndata,nvar)

! Fit polynomial to energies.
call PolFit(ipow,nvar,qvar,yin,ndata,coef,nterm,stand_dev,max_err,diff_vec,use_weight)

!write(u6,*) ndata

!do i=1,ndata
!  write(u6,'(4f15.8,es15.6)') (var(i,j),j=1,nvar),yin(i),diff_vec(i)
!end do

call mma_deallocate(diff_vec)

if (find_minimum) then
  call Optimize(ipow,qvar,coef,x,energy,Hess,nterm,nvar,ndata)
else
  call x_to_qvar(x,ref,qref,alpha,trfName,nvar)
end if

! Calculate gradient and quadratic, cubic and quartic force constants.
call funcval(x,coef,ipow,energy,nterm,nvar)
call gradient(x,coef,ipow,grad,nterm,nvar)
call Hessian(x,coef,ipow,Hess,nterm,nvar)
if (max_term > 2) then
  call thirdDer(x,coef,ipow,D3,nterm,nvar)
else
  D3(:,:,:) = Zero
end if
if (max_term > 3) then
  call fourthDer(x,coef,ipow,D4,nterm,nvar)
else
  D4(:,:,:,:) = Zero
end if

! Transform back to original coordinates.
call qvar_to_var(x,grad,Hess,D3,D4,ref,qref,trfName,alpha,nvar)

call mma_deallocate(alpha)
call mma_deallocate(qvar)
call mma_deallocate(ref)
call mma_deallocate(qref)

end subroutine PotFit
