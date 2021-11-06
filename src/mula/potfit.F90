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
!    var          : Real*8 two dimensional array - coordinates
!                   for which we know the energy.
!    yin          : Real*8 array - energy in points given in var.
!    trfName      : Character array - transformation associated with each
!                   internal coordinate.
!
!  Output:
!    coef         : Real*8 two dimensional array - coefficients
!                   for each term specified by ipow.
!    x            : Real*8 array - geometry in internal coordinates.
!    energy       : Real*8 - energy in minimum.
!    Hess         : Real*8 two dimensional array - Hessian in
!                   minimum.
!    D3           : Real*8 three dimensional array - cubic force
!                   constants.
!    D4           : Real*8 four dimensional array - quartic force
!                   constants.
!    stand_dev    : Real*8 variable - standard deviation of fitted
!                   values in Subroutine PolFit.
!    max_err      : Real*8 variable - maximum error of fitted
!                   values in Subroutine PolFit.
!    find_minimum : Logical variable - whether or not to optimize structure.
!    max_term     : Integer - highest power if term in polynomial fit.
!    use_weight   : Logical variable - whether or not to use weights to
!                   make the minimum more important in the fit.
!
!  Uses:
!    LinAlg
!    VibMod
!    RandomMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!use LinAlg
!use VibMod
!use RandomMod
implicit real*8(a-h,o-z)
integer iPow(nterm,nvar)
real*8 var(ndata,nvar), yin(ndata)
real*8 coef(nterm,1)
real*8 x(nosc)
real*8 grad(nosc)
real*8 Hess(l_Hess_1,l_Hess_2)
real*8 D3(l_D,l_D,l_D)
real*8 D4(l_D,l_D,l_D,l_D)
character*80 trfName(nvar)
character*32 trfcode
real*8 stand_dev, max_err
logical find_minimum
logical use_weight
logical fit_alpha
real*8 diff_vec(ndata)
character*32 Inline
#include "WrkSpc.fh"


! Initialize.
do iv=1,nterm
  coef(iv,1) = 0.0d0
end do
call GetMem('qvar','Allo','Real',ipqvar,ndata*nvar)
call GetMem('ref','Allo','Real',ipref,nvar)
call GetMem('qref','Allo','Real',ipqref,nvar)

! Check if a fit of alpha is wanted.
fit_alpha = .false.
call GetMem('alphaStart','Allo','Real',ipalphaStart,nvar)
call GetMem('alphaStop','Allo','Real',ipalphaStop,nvar)
call GetMem('alpha','Allo','Real',ipalpha,nvar)
call GetMem('alpha0','Allo','Real',ipalpha0,nvar)
call GetMem('alphaTemp','Allo','Real',ipalphaTemp,nvar)

call GetMem('alphaRad','Allo','Real',ipalphaRad,nvar)
call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStart),1)
call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStop),1)
!alphaStart =-1.0d0
!alphaStop =-1.0d0
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
    read(Inline(1:istop),*) Work(ipalphaStart+ivar-1),Work(ipalphaStop+ivar-1)
    Work(ipalphaRad+ivar-1) = (Work(ipalphaStop+ivar-1)-Work(ipalphaStart+ivar-1))/2.0d0
    Work(ipalpha0+ivar-1) = Work(ipalphaStart+ivar-1)+Work(ipalphaRad+ivar-1)
  end if
end do

if (fit_alpha) then

  ! Determine which alpha radius is the largest.
  alphaMax = Work(ipalphaRad)
  ivarMax = 1
  do ivar=2,nvar
    if (Work(ipalphaRad+ivar-1) > alphaMax) then
      alphaMax = Work(ipalphaRad+ivar-1)
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
  call GetMem('RandVec','Allo','Real',ipRandVec,length)

  call Ranmar(Work(ipRandVec),length)

  ! Simulated annealing.
  stand_dev_best = 100.0d0
  iRandNum = 1
  do while (Work(ipalphaRad+ivarMax-1) > 1.0d-1)
    do iv=1,nvar
      Work(ipalphaTemp+iv-1) = Work(ipalpha0+iv-1)
    end do
    do iPoints=1,1000
      do ivar=1,nvar
        if (Work(ipalphaStart+ivar-1) > 0.0d0) then
          Work(ipalpha+ivar-1) = Work(ipalpha0+ivar-1)+Work(ipalphaRad+ivar-1)*(2.0d0*Work(ipRandVec+iRandNum+ivar-1-1)-1.0d0)
        end if
      end do
      iRandNum = iRandNum+nvar
      call var_to_qvar(var,Work(ipqvar),Work(ipref),Work(ipqref),Work(ipalpha),trfName,ndata,nvar)
      call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,coef,nterm,stand_dev,max_err,diff_vec,use_weight)
      if (stand_dev < stand_dev_best) then
        stand_dev_best = stand_dev
        do iv=1,nvar
          Work(ipalphaTemp+iv-1) = Work(ipalpha+iv-1)
        end do
      end if
    end do
    do iv=1,nvar
      Work(ipalpha0+iv-1) = Work(ipalphaTemp+iv-1)
    end do
    write(6,*) (Work(ipalpha0+iv-1),iv=1,nvar),stand_dev_best
    do iv=1,nvar
      Work(ipalphaRad+iv-1) = Work(ipalphaRad+iv-1)/2.0d0
    end do
  end do
  do iv=1,nvar
    Work(ipalpha+iv-1) = Work(ipalpha0+iv-1)
  end do
  write(6,*) Work(ipalpha),Work(ipalpha+1)
  call GetMem('RandVec','Free','Real',ipRandVec,length)

end if

call GetMem('alpha0','Free','Real',ipalpha0,nvar)
call GetMem('alphaTemp','Free','Real',ipalphaTemp,nvar)

call GetMem('alphaRad','Free','Real',ipalphaRad,nvar)
call GetMem('alphaStart','Free','Real',ipalphaStart,nvar)
call GetMem('alphaStop','Free','Real',ipalphaStop,nvar)

! Transform coordinates to get better numerical fit.
call var_to_qvar(var,Work(ipqvar),Work(ipref),Work(ipqref),Work(ipalpha),trfName,ndata,nvar)

! Fit polynomial to energies.
call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,coef,nterm,stand_dev,max_err,diff_vec,use_weight)

!write(6,*) ndata

!do i=1,ndata
!  write(6,'(4f15.8,es15.6)') (var(i,j),j=1,nvar),yin(i),diff_vec(i)
!end do

if (find_minimum) then
  call Optimize(ipow,Work(ipqvar),coef,x,energy,Hess,nterm,nvar,ndata)
else
  call x_to_qvar(x,Work(ipref),Work(ipqref),Work(ipalpha),trfName,nvar)
end if

! Calculate gradient and quadratic, cubic and quartic force constants.
call funcval(x,coef,ipow,energy,nterm,nvar)
call gradient(x,coef,ipow,grad,nterm,nvar)
call Hessian(x,coef,ipow,Hess,nterm,nvar)
if (max_term > 2) then
  call thirdDer(x,coef,ipow,D3,nterm,nvar)
else
  call dcopy_(l_D*l_D*l_D,[0.0d0],0,D3,1)
  !D3 = 0.0d0
end if
if (max_term > 3) then
  call fourthDer(x,coef,ipow,D4,nterm,nvar)
else
  !D4 = 0.0d0
  call dcopy_(l_D*l_D*l_D*l_D,[0.0d0],0,D4,1)
end if

! Transform back to original coordinates.
call qvar_to_var(var,x,grad,Hess,D3,D4,Work(ipref),Work(ipqref),trfName,Work(ipalpha),max_term,ndata,nvar)

call GetMem('alpha','Free','Real',ipalpha,nvar)
call GetMem('qvar','Free','Real',ipqvar,ndata*nvar)
call GetMem('ref','Free','Real',ipref,nvar)
call GetMem('qref','Free','Real',ipqref,nvar)

end subroutine PotFit
