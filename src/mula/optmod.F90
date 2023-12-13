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
!    PotFit       (ipow,var,yin,coef,x,energy,Hess,D3,D4,trfName,stand_dev,max_err,find_minimum)
!    Optimize     (ipow,var,coef,x,energy,Hess,D3,D4)
!    var_to_qvar  (var,qvar,ref,qref,alpha,trfName)
!    x_to_qvar    (x,ref,qref,alpha,trfName)
!    qvar_to_var  (var,x,grad,Hess,D3,D4,ref,qref,trfName,alpha)
!    funcval      (x,coef,ipow) Result(sum)
!    gradient     (x,coef,ipow) Result(grad)
!    Hessian      (x,coef,ipow) Result(Hess)
!    thirdDer     (x,coef,ipow) Result(D3)
!    fourthDer    (x,coef,ipow) Result(D4)
!    ShiftHess    (Hess,shift)
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
!    ShiftHess
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

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
real(kind=wp), external :: Random_Molcas

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
    !call calcNorm(delta(:,1),delta_norm)
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

subroutine var_to_qvar(var,qvar,ref,qref,alpha,trfName,ndata,nvar)
!  Purpose:
!    Transform coordinates given in input using tranformation specified in input.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndata, nvar
real(kind=wp), intent(in) :: var(ndata,nvar)
real(kind=wp), intent(out) :: qvar(ndata,nvar), ref(nvar), qref(nvar)
real(kind=wp), intent(inout) :: alpha(nvar)
character(len=80), intent(in) :: trfName(nvar)
integer(kind=iwp) :: ia, ic, id, idata, ie, ifit, ir, is, istart, istop, ivar, ix
real(kind=wp) :: angsc, rsum, v
character(len=32) :: Inline, trfCode
real(kind=wp), allocatable :: par(:,:)

call mma_allocate(par,ndata,nvar,label='par')

do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ix = index(trfcode,'AS IT IS')
  ia = index(trfcode,'-AVG')
  ie = index(trfcode,'EXP')
  ir = index(trfcode,'RAD')
  id = index(trfcode,'DEG')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  angsc = One
  if (ir > 0) angsc = One
  if (id > 0) angsc = deg2rad
  do idata=1,ndata
    v = var(idata,ivar)
    if (ic > 0) then
      par(idata,ivar) = cos(angsc*v)
    else if (is > 0) then
      par(idata,ivar) = sin(angsc*v)
    else if ((ix > 0) .or. (ie > 0)) then
      par(idata,ivar) = v
    else
      write(u6,*) ' TRFCODE ERROR.'
      call abend()
    end if
  end do
  ! Calculate refrence value.
  rsum = Zero
  do idata=1,ndata
    rsum = rsum+var(idata,ivar)
  end do
  ref(ivar) = rsum/ndata
  if (ic > 0) then
    ref(ivar) = cos(angsc*ref(ivar))
  else if (is > 0) then
    ref(ivar) = sin(angsc*ref(ivar))
  end if

  ! Subtract reference value if requested.
  if (ia > 0) then
    do idata=1,ndata
      v = par(idata,ivar)
      if ((ic > 0) .or. (is > 0)) then
        par(idata,ivar) = ref(ivar)-v
      else
        par(idata,ivar) = v-ref(ivar)
      end if
    end do
  end if
end do

! Transform coordinates.
do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ie = index(trfcode,'EXP')
  ifit = index(trfcode,'FIT')
  if ((ie > 0) .and. (ifit == 0)) then
    Inline = trfName(ivar)(1:32)
    istart = index(Inline,'ALPHA=')
    istart = istart+6
    Inline = trfCode(istart:len(trfCode))
    istop = index(Inline,' ')
    istop = istop-1
    read(Inline(1:istop),*) alpha(ivar)
  end if
  do idata=1,ndata
    if (ie > 0) then
      qvar(idata,ivar) = One-exp(-alpha(ivar)*par(idata,ivar))
    else
      qvar(idata,ivar) = par(idata,ivar)
    end if
  end do
end do

call mma_deallocate(par)

! Calculate reference value of transformed coordinates.
do ivar=1,nvar
  rsum = Zero
  do idata=1,ndata
    rsum = rsum+qvar(idata,ivar)
  end do
  qref(ivar) = rsum/ndata
end do

! Subtract reference value from transformed coordinates.
do ivar=1,nvar
  do idata=1,ndata
    qvar(idata,ivar) = qvar(idata,ivar)-qref(ivar)
  end do
end do

end subroutine var_to_qvar

subroutine x_to_qvar(x,ref,qref,alpha,trfName,nDimX)
!  Purpose:
!    Transform the coordinates of a given point.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDimX
real(kind=wp), intent(inout) :: x(nDimX), alpha(nDimx)
real(kind=wp), intent(in) :: ref(nDimX), qref(nDimX)
character(len=80), intent(in) :: trfName(nDimX)
integer(kind=iwp) :: ia, ic, ie, is, istart, istop, ivar
real(kind=wp) :: v
character(len=32) :: Inline, trfCode
real(kind=wp), allocatable :: par(:)

call mma_allocate(par,nDimX,label='par')

do ivar=1,nDimX
  trfcode = trfName(ivar)(1:32)
  ia = index(trfcode,'-AVG')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  if (ic > 0) then
    par(ivar) = cos(x(ivar))
  else if (is > 0) then
    par(ivar) = sin(x(ivar))
  else
    par(ivar) = x(ivar)
  end if
  if (ia > 0) then
    v = par(ivar)
    if ((ic > 0) .or. (is > 0)) then
      par(ivar) = ref(ivar)-v
    else
      par(ivar) = v-ref(ivar)
    end if
  end if
end do

! Transform coordinates.
do ivar=1,nDimX
  trfcode = trfName(ivar)(1:32)
  ie = index(trfcode,'EXP')
  if (ie > 0) then
    Inline = trfName(ivar)(1:32)
    istart = index(Inline,'ALPHA=')
    istart = istart+6
    Inline = trfCode(istart:len(trfCode))
    istop = index(Inline,' ')
    istop = istop-1
    read(Inline(1:istop),*) alpha(ivar)
  end if
  if (ie > 0) then
    x(ivar) = One-exp(-alpha(ivar)*par(ivar))
  else
    x(ivar) = par(ivar)
  end if
end do

call mma_deallocate(par)

! Subtract reference value from transformed coordinates.
do ivar=1,nDimX
  x(ivar) = x(ivar)-qref(ivar)
end do

end subroutine x_to_qvar

subroutine qvar_to_var(x,grad,Hess,D3,D4,ref,qref,trfName,alpha,nvar)
!  Purpose:
!    Tranform coordinates, gradient, Hessian, third derivatives and
!    fourth derivates back to the coordinates originally specified.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Four, Six
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nvar
real(kind=wp), intent(inout) :: x(nvar), grad(nvar), Hess(nvar,nvar), D3(nvar,nvar,nvar), D4(nvar,nvar,nvar,nvar)
real(kind=wp), intent(in) :: ref(nvar), qref(nvar), alpha(nvar)
character(len=80), intent(in) :: trfName(nvar)
integer(kind=iwp) :: i, ia, ic, ie, is, ivar, ix, j, k, l
character(len=32) :: trfCode
logical(kind=iwp) :: ijEq, jkEq, klEq
real(kind=wp), allocatable :: D3trans(:,:,:), D4trans(:,:,:,:), q(:), s(:), t(:), Temp(:,:), tempgrad(:), u(:), v(:)

call mma_allocate(s,nvar,label='s')
call mma_allocate(t,nvar,label='t')
call mma_allocate(u,nvar,label='u')
call mma_allocate(v,nvar,label='v')
call mma_allocate(q,nvar,label='q')
call mma_allocate(tempgrad,nvar,label='tempgrad')
call mma_allocate(Temp,nvar,nvar,label='Temp')
call mma_allocate(D3trans,nvar,nvar,nvar,label='D3trans')
call mma_allocate(D4trans,nvar,nvar,nvar,nvar,label='D4trans')

x(:) = x+qref
do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ix = index(trfcode,'AS IT IS')
  ia = index(trfcode,'-AVG')
  ie = index(trfcode,'EXP')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  if (ic > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)-ref(ivar)
    end if
    q(ivar) = acos(-x(ivar))
    t(ivar) = sin(q(ivar))
    u(ivar) = cos(q(ivar))
    v(ivar) = -sin(q(ivar))
    s(ivar) = -cos(q(ivar))
  else if (is > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)-ref(ivar)
    end if
    q(ivar) = asin(-x(ivar))
    t(ivar) = -cos(q(ivar))
    u(ivar) = sin(q(ivar))
    v(ivar) = cos(q(ivar))
    s(ivar) = -sin(q(ivar))
  else if (ie > 0) then
    q(ivar) = -log(One-x(ivar))/alpha(ivar)
    t(ivar) = alpha(ivar)*(One-x(ivar))
    u(ivar) = -alpha(ivar)**2*(One-x(ivar))
    v(ivar) = alpha(ivar)**3*(One-x(ivar))
    s(ivar) = -alpha(ivar)**4*(One-x(ivar))
    if (ia > 0) then
      q(ivar) = q(ivar)+ref(ivar)
    end if
  else if (ix > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)+ref(ivar)
    end if
    q(ivar) = x(ivar)
    t(ivar) = One
    u(ivar) = Zero
    v(ivar) = Zero
    s(ivar) = Zero
  else
    write(u6,*) ' TRFCODE ERROR.'
    call abend()
  end if
end do

! Transform gradient.
do i=1,nvar
  tempgrad(i) = t(i)*grad(i)
end do

! Transform Hessian.
do i=1,nvar
  do j=i,nvar
    if (i == j) then
      Temp(i,i) = t(i)**2*Hess(i,i)+u(i)*grad(i)
    else
      Temp(i,j) = t(i)*t(j)*Hess(i,j)
    end if
  end do
end do
do i=1,nvar
  do j=i,nvar
    Temp(j,i) = Temp(i,j)
  end do
end do

! Transform third derivatives.
do i=1,nvar
  do j=i,nvar
    ijEq = (i == j)
    do k=j,nvar
      jkEq = (j == k)
      if (ijEq .and. jkEq) then
        D3trans(i,i,i) = t(i)*t(i)*t(i)*D3(i,i,i)+Three*u(i)*t(i)*Hess(i,i)+v(i)*grad(i)
      else if (ijEq .and. (.not. jkEq)) then
        D3trans(i,i,k) = t(i)*t(i)*t(k)*D3(i,i,k)+u(i)*t(k)*Hess(i,k)
      else if ((.not. ijEq) .and. jkEq) then
        D3trans(i,j,j) = t(i)*t(j)*t(j)*D3(i,j,j)+t(i)*u(j)*Hess(i,j)
      else if ((.not. ijEq) .and. (.not. jkEq)) then
        D3trans(i,j,k) = t(i)*t(j)*t(k)*D3(i,j,k)
      end if
    end do
  end do
end do
do i=1,nvar
  do j=i,nvar
    do k=j,nvar
      D3trans(i,k,j) = D3trans(i,j,k)
      D3trans(j,i,k) = D3trans(i,j,k)
      D3trans(j,k,i) = D3trans(i,j,k)
      D3trans(k,i,j) = D3trans(i,j,k)
      D3trans(k,j,i) = D3trans(i,j,k)
    end do
  end do
end do

! Transform fourth derivatives.
do i=1,nvar
  do j=i,nvar
    ijEq = (i == j)
    do k=j,nvar
      jkEq = (j == k)
      do l=k,nvar
        klEq = (k == l)
        if (ijEq .and. jkEq .and. klEq) then
          D4trans(i,i,i,i) = t(i)**4*D4(i,i,i,i)+Six*u(i)*t(i)**2*D3(i,i,i)+Four*v(i)*t(i)*Hess(i,i)+Three*u(i)*u(i)*Hess(i,i)+ &
                             s(i)*grad(i)
        else if (ijEq .and. jkEq .and. (.not. klEq)) then
          D4trans(i,i,i,l) = (t(i)*t(i)*t(i)*D4(i,i,i,l)+Three*u(i)*t(i)*D3(i,i,l)+v(i)*Hess(i,l))*t(l)
        else if ((.not. ijEq) .and. jkEq .and. klEq) then
          D4trans(i,j,j,j) = t(i)*(t(j)*t(j)*t(j)*D4(i,j,j,j)+Three*u(j)*t(j)*D3(i,j,j)+v(j)*Hess(i,j))
        else if ((.not. ijEq) .and. jkEq .and. (.not. klEq)) then
          D4trans(i,j,j,l) = t(i)*(t(j)*t(j)*D4(i,j,j,l)+u(j)*D3(i,j,l))*t(l)
        else if (ijEq .and. (.not. jkEq) .and. (.not. klEq)) then
          D4trans(i,i,k,l) = (t(i)*t(i)*D4(i,i,k,l)+u(i)*D3(i,k,l))*t(k)*t(l)
        else if ((.not. ijEq) .and. (.not. jkEq) .and. klEq) then
          D4trans(i,j,k,k) = (t(k)*t(k)*D4(i,j,k,k)+u(k)*D3(i,j,k))*t(i)*t(j)
        else if (ijEq .and. klEq) then
          D4trans(i,i,k,k) = t(i)*t(i)*t(k)*t(k)*D4(i,i,k,k)+t(i)*t(i)*u(k)*D3(i,i,k)+u(i)*t(k)*t(k)*D3(i,k,k)+u(i)*u(k)*Hess(i,k)
        else if ((.not. ijEq) .and. (.not. jkEq) .and. (.not. klEq)) then
          D4trans(i,j,k,l) = t(i)*t(j)*t(k)*t(l)*D4(i,j,k,l)
        end if
      end do
    end do
  end do
end do
do i=1,nvar
  do j=i,nvar
    do k=j,nvar
      do l=k,nvar
        D4trans(i,k,j,l) = D4trans(i,j,k,l)
        D4trans(i,j,l,k) = D4trans(i,j,k,l)
        D4trans(i,l,k,j) = D4trans(i,j,k,l)
        D4trans(i,l,j,k) = D4trans(i,j,k,l)
        D4trans(i,k,l,j) = D4trans(i,j,k,l)

        D4trans(j,i,k,l) = D4trans(i,j,k,l)
        D4trans(j,k,i,l) = D4trans(i,j,k,l)
        D4trans(j,i,l,k) = D4trans(i,j,k,l)
        D4trans(j,l,k,i) = D4trans(i,j,k,l)
        D4trans(j,l,i,k) = D4trans(i,j,k,l)
        D4trans(j,k,l,i) = D4trans(i,j,k,l)

        D4trans(k,i,j,l) = D4trans(i,j,k,l)
        D4trans(k,j,i,l) = D4trans(i,j,k,l)
        D4trans(k,i,l,j) = D4trans(i,j,k,l)
        D4trans(k,l,j,i) = D4trans(i,j,k,l)
        D4trans(k,l,i,j) = D4trans(i,j,k,l)
        D4trans(k,j,l,i) = D4trans(i,j,k,l)

        D4trans(l,i,j,k) = D4trans(i,j,k,l)
        D4trans(l,j,i,k) = D4trans(i,j,k,l)
        D4trans(l,k,j,i) = D4trans(i,j,k,l)
        D4trans(l,i,k,j) = D4trans(i,j,k,l)
        D4trans(l,k,i,j) = D4trans(i,j,k,l)
        D4trans(l,j,k,i) = D4trans(i,j,k,l)
      end do
    end do
  end do
end do

! Assign transformed values.
x(:) = q
grad(:) = tempgrad

Hess(:,:) = Temp

D3(:,:,:) = D3trans

D4(:,:,:,:) = D4trans

call mma_deallocate(s)
call mma_deallocate(t)
call mma_deallocate(u)
call mma_deallocate(v)
call mma_deallocate(q)
call mma_deallocate(tempgrad)
call mma_deallocate(Temp)
call mma_deallocate(D3trans)
call mma_deallocate(D4trans)

end subroutine qvar_to_var

subroutine funcval(x,coef,ipow,fval,nterm,nvar)
!  Purpose:
!    Return function value at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar)
real(kind=wp), intent(in) :: x(nvar), coef(nterm)
real(kind=wp), intent(out) :: fval
integer(kind=iwp) :: iterm, ivar, nsum
real(kind=wp) :: prod, rsum

rsum = Zero
do iterm=1,nterm
  prod = One
  do ivar=1,nvar
    nsum = ipow(iterm,ivar)
    prod = prod*(x(ivar)**(nsum))
  end do
  rsum = rsum+coef(iterm)*prod
end do
fval = rsum

end subroutine funcval

subroutine gradient(x,coef,ipow,grad,nterm,nvar)
!  Purpose:
!    Return gradient at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar)
real(kind=wp), intent(in) :: x(nvar), coef(nvar)
real(kind=wp), intent(out) :: grad(nvar)
integer(kind=iwp) :: iterm, ivar, jvar, nder, nsum
real(kind=wp) :: prod, rfactor, rsum
logical(kind=iwp) :: ijEq

do ivar=1,nvar
  rsum = Zero
  do iterm=1,nterm
    prod = One
    do jvar=1,nvar
      ijEq = (ivar == jvar)
      if ((ipow(iterm,jvar) >= 1) .and. ijEq) then
        nder = 1
        nsum = ipow(iterm,jvar)-nder
      else if ((ipow(iterm,jvar) >= 0) .and. (.not. ijEq)) then
        nder = 0
        nsum = ipow(iterm,jvar)
      else
        nder = -1
        nsum = 0
      end if
      call factor(nsum,nder,rfactor)
      prod = prod*rfactor*(x(jvar)**(nsum))
    end do
    rsum = rsum+coef(iterm)*prod
  end do
  grad(ivar) = rsum
end do

end subroutine gradient

subroutine Hessian(x,coef,ipow,Hess,nterm,nvar)
!  Purpose:
!    Return Hessian at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar)
real(kind=wp), intent(in) :: x(nvar), coef(nvar)
real(kind=wp), intent(out) :: Hess(nvar,nvar)
integer(kind=iwp) :: iterm, ivar, jvar, kvar, nder, nsum
real(kind=wp) :: prod, rfactor, rsum
logical(kind=iwp) :: ijEq, ikEq, jkEq

do ivar=1,nvar
  do jvar=ivar,nvar
    ijEq = (ivar == jvar)
    rsum = Zero
    do iterm=1,nterm
      prod = One
      do kvar=1,nvar
        ikEq = (ivar == kvar)
        jkEq = (jvar == kvar)
        if ((ipow(iterm,kvar) >= 2) .and. ikEq .and. jkEq) then
          nder = 2
          nsum = ipow(iterm,kvar)-nder
        else if ((ipow(iterm,kvar) >= 1) .and. (ikEq .or. jkEq) .and. (.not. ijEq)) then
          nder = 1
          nsum = ipow(iterm,kvar)-nder
        else if ((ipow(iterm,kvar) >= 0) .and. (.not. ikEq) .and. (.not. jkEq)) then
          nder = 0
          nsum = ipow(iterm,kvar)
        else
          nder = -1
          nsum = 0
        end if
        call factor(nsum,nder,rfactor)
        prod = prod*rfactor*(x(kvar)**(nsum))
      end do
      rsum = rsum+coef(iterm)*prod
    end do
    Hess(ivar,jvar) = rsum
    Hess(jvar,ivar) = rsum
  end do
end do

end subroutine Hessian

subroutine thirdDer(x,coef,ipow,D3,nterm,nvar)
!  Purpose:
!    Return third derivatives at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar)
real(kind=wp), intent(in) :: x(nvar), coef(nvar)
real(kind=wp), intent(out) :: D3(nvar,nvar,nvar)
integer(kind=iwp) :: iterm, ivar, jvar, kvar, lvar, nder, nsum
real(kind=wp) :: prod, rfactor, rsum
logical(kind=iwp) :: ilEq, jlEq, klEq

D3(:,:,:) = Zero
do ivar=1,nvar
  do jvar=ivar,nvar
    do kvar=jvar,nvar
      rsum = Zero
      do iterm=1,nterm
        prod = One
        do lvar=1,nvar
          ilEq = (ivar == lvar)
          jlEq = (jvar == lvar)
          klEq = (kvar == lvar)
          if ((ipow(iterm,lvar) >= 3) .and. ilEq .and. jlEq .and. klEq) then
            nder = 3
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 2) .and. ((ilEq .and. jlEq .and. (.not. klEq)) .or. &
                                                  (jlEq .and. klEq .and. (.not. ilEq)))) then
            nder = 2
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 1) .and. ((ilEq .and. (.not. jlEq) .and. (.not. klEq)) .or. &
                                                  (jlEq .and. (.not. ilEq) .and. (.not. klEq)) .or. &
                                                  (klEq .and. (.not. ilEq) .and. (.not. jlEq)))) then
            nder = 1
            nsum = ipow(iterm,lvar)-nder
          else if ((ipow(iterm,lvar) >= 0) .and. (.not. ilEq) .and. (.not. jlEq) .and. (.not. klEq)) then
            nder = 0
            nsum = ipow(iterm,lvar)
          else
            nder = -1
            nsum = 0
          end if
          call factor(nsum,nder,rfactor)
          prod = prod*rfactor*(x(lvar)**(nsum))
        end do
        rsum = rsum+coef(iterm)*prod
      end do
      D3(ivar,jvar,kvar) = rsum
      D3(kvar,ivar,jvar) = rsum
      D3(jvar,kvar,ivar) = rsum
      D3(jvar,ivar,kvar) = rsum
      D3(kvar,jvar,ivar) = rsum
      D3(ivar,kvar,jvar) = rsum
    end do
  end do
end do

end subroutine thirdDer

subroutine fourthDer(x,coef,ipow,D4,nterm,nvar)
!  Purpose:
!    Return fourth derivatives at x.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nterm, nvar, ipow(nterm,nvar)
real(kind=wp), intent(in) :: x(nvar), coef(nvar)
real(kind=wp), intent(out) :: D4(nvar,nvar,nvar,nvar)
integer(kind=iwp) :: iterm, ivar, jvar, kvar, lvar, mvar, nder, nsum
real(kind=wp) :: prod, rfactor, rsum
logical(kind=iwp) :: imEq, jmEq, kmEq, lmEq

D4(:,:,:,:) = Zero
do ivar=1,nvar
  do jvar=ivar,nvar
    do kvar=jvar,nvar
      do lvar=kvar,nvar
        rsum = Zero
        do iterm=1,nterm
          prod = One
          do mvar=1,nvar
            imEq = (ivar == mvar)
            jmEq = (jvar == mvar)
            kmEq = (kvar == mvar)
            lmEq = (lvar == mvar)
            if ((ipow(iterm,mvar) >= 4) .and. imEq .and. jmEq .and. kmEq .and. lmEq) then
              nder = 4
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 3) .and. ((imEq .and. jmEq .and. kmEq .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. kmEq .and. lmEq .and. (.not. imEq)))) then
              nder = 3
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 2) .and. ((imEq .and. jmEq .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. kmEq .and. (.not. imEq) .and. (.not. lmEq)) .or. &
                                                    (kmEq .and. lmEq .and. (.not. imEq) .and. (.not. jmEq)))) then
              nder = 2
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 1) .and. ((imEq .and. (.not. jmEq) .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (jmEq .and. (.not. imEq) .and. (.not. kmEq) .and. (.not. lmEq)) .or. &
                                                    (kmEq .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. lmEq)) .or. &
                                                    (lmEq .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. kmEq)))) then
              nder = 1
              nsum = ipow(iterm,mvar)-nder
            else if ((ipow(iterm,mvar) >= 0) .and. (.not. imEq) .and. (.not. jmEq) .and. (.not. kmEq) .and. (.not. lmEq)) then
              nder = 0
              nsum = ipow(iterm,mvar)
            else
              nder = -1
              nsum = 0
            end if
            call factor(nsum,nder,rfactor)
            prod = prod*rfactor*(x(mvar)**(nsum))
          end do
          rsum = rsum+coef(iterm)*prod
        end do
        D4(ivar,jvar,kvar,lvar) = rsum
        D4(ivar,kvar,jvar,lvar) = rsum
        D4(ivar,jvar,lvar,kvar) = rsum
        D4(ivar,lvar,kvar,jvar) = rsum
        D4(ivar,lvar,jvar,kvar) = rsum
        D4(ivar,kvar,lvar,jvar) = rsum

        D4(jvar,ivar,kvar,lvar) = rsum
        D4(jvar,kvar,ivar,lvar) = rsum
        D4(jvar,ivar,lvar,kvar) = rsum
        D4(jvar,lvar,kvar,ivar) = rsum
        D4(jvar,lvar,ivar,kvar) = rsum
        D4(jvar,kvar,lvar,ivar) = rsum

        D4(kvar,ivar,jvar,lvar) = rsum
        D4(kvar,jvar,ivar,lvar) = rsum
        D4(kvar,ivar,lvar,jvar) = rsum
        D4(kvar,lvar,jvar,ivar) = rsum
        D4(kvar,lvar,ivar,jvar) = rsum
        D4(kvar,jvar,lvar,ivar) = rsum

        D4(lvar,ivar,jvar,kvar) = rsum
        D4(lvar,jvar,ivar,kvar) = rsum
        D4(lvar,kvar,jvar,ivar) = rsum
        D4(lvar,ivar,kvar,jvar) = rsum
        D4(lvar,kvar,ivar,jvar) = rsum
        D4(lvar,jvar,kvar,ivar) = rsum
      end do
    end do
  end do
end do

end subroutine fourthDer

subroutine ShiftHess(Hess,shift,nDim,nDim2)
!  Purpose:
!    Shifts Hessian to make it positive definite.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nDim2
real(kind=wp), intent(inout) :: Hess(nDim,nDim2)
logical(kind=iwp), intent(out) :: shift
integer(kind=iwp) :: i, j, k
real(kind=wp) :: eigen_min, eps
real(kind=wp), allocatable :: Hess_lowT(:), U(:,:)

call mma_allocate(U,nDim,nDim,label='U')
call mma_allocate(Hess_lowT,nDim*(nDim+1)/2,label='Hess_lowT')

! Initialize.

k = 0
do i=1,nDim
  do j=1,i
    k = k+1
    Hess_lowT(k) = Hess(i,j)
  end do
end do
call unitmat(U,nDim)
call Jacob(Hess_lowT,U,nDim,nDim)
call Jacord(Hess_lowT,U,nDim,nDim)
eigen_min = Hess_lowT(1)
shift = eigen_min < Zero
if (shift) then
  eps = 2*eigen_min
  do i=1,nDim
    Hess(i,i) = Hess(i,i)-eps
  end do
end if
call mma_deallocate(U)
call mma_deallocate(Hess_lowT)

end subroutine ShiftHess

!end module OptMod
