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
!!-----------------------------------------------------------------------!
!!
!       Module OptMod
!!
!!  Contains:
!!    PotFit       (ipow,var,yin,coef,x,energy,Hess,D3,D4,trfName,
!!                  stand_dev,max_err,find_minimum)
!!    Optimize     (ipow,var,coef,x,energy,Hess,D3,D4)
!!    var_to_qvar  (var,qvar,ref,qref,alpha,trfName)
!!    qvar_to_var  (var,x,grad,Hess,D3,D4,ref,qref,trfName,alpha)
!!    funcval      (x,coef,ipow) Result(sum)
!!    gradient     (x,coef,ipow) Result(grad)
!!    Hessian      (x,coef,ipow) Result(Hess)
!!    thirdDer     (x,coef,ipow) Result(D3)
!!    fourthDer    (x,coef,ipow) Result(D4)
!!    ShiftHess    (Hess,shift)
!!
!!  Uses:
!!    Constants
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
!!-----------------------------------------------------------------------!
!!
!       Contains

!!-----------------------------------------------------------------------!
!!
      Subroutine PotFit(nterm,nvar,ndata,ipow,var,yin,                  &
     &       coef,x,nOsc,energy,                                        &
     &       grad,Hess,D3,D4,trfName,stand_dev,max_err,                 &
     &       find_minimum,max_term,use_weight,                          &
     &       l_Hess_1,l_Hess_2,l_D)
!!
!!  Purpose:
!!    Fit polynomial to data and then depending upun the variable
!!    find_minimum, we either calculate the energy, gradient, Hessian,
!!    and qubic and quartic force constants in x or we optimize the
!!    structure before we perform the calculations.
!!
!!  Input:
!!    ipow         : Two dimensional integer array - terms in polynomial.
!!    var          : Real*8 two dimensional array - coordinates
!!                   for which we know the energy.
!!    yin          : Real*8 array - energy in points given in var.
!!    trfName      : Character array - transformation associated with each
!!                   internal coordinate.
!!
!!  Output:
!!    coef         : Real*8 two dimensional array - coefficients
!!                   for each term specified by ipow.
!!    x            : Real*8 array - geometry in internal coordinates.
!!    energy       : Real*8 - energy in minimum.
!!    Hess         : Real*8 two dimensional array - Hessian in
!!                   minimum.
!!    D3           : Real*8 three dimensional array - cubic force
!!                   constants.
!!    D4           : Real*8 four dimensional array - quartic force
!!                   constants.
!!    stand_dev    : Real*8 variable - standard deviation of fitted
!!                   values in Subroutine PolFit.
!!    max_err      : Real*8 variable - maximum error of fitted
!!                   values in Subroutine PolFit.
!!    find_minimum : Logical variable - whether or not to optimize structure.
!!    max_term     : Integer - highest power if term in polynomial fit.
!!    use_weight   : Logical variable - whether or not to use weights to
!!                   make the minimum more important in the fit.
!!
!!  Uses:
!!    LinAlg
!!    VibMod
!!    RandomMod
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
!       Use LinAlg
!       Use VibMod
!! Use RandomMod
      Implicit Real*8 ( a-h,o-z )
      Integer iPow(nterm,nvar)
      Real*8   var(ndata,nvar),yin(ndata)
      Real*8   coef(nterm,1)
      Real*8  x(nosc)
      Real*8  grad(nosc)
      Real*8  Hess(l_Hess_1,l_Hess_2)
      Real*8  D3(l_D,l_D,l_D)
      Real*8  D4(l_D,l_D,l_D,l_D)
      Character*80 trfName(nvar)
      Character*32  trfcode
      Real*8    stand_dev,max_err
      Logical   find_minimum
      Logical   use_weight
      Logical   fit_alpha
      Real*8 diff_vec(ndata)

      Character*32   Inline
#include "WrkSpc.fh"

!!
!!---- Initialize.
      do iv=1,nterm
      coef(iv,1) = 0.0d0
      enddo
      Call GetMem('qvar','Allo','Real',ipqvar,ndata*nvar)
      Call GetMem('ref','Allo','Real',ipref,nvar)
      Call GetMem('qref','Allo','Real',ipqref,nvar)

!!
!!---- Check if a fit of alpha is wanted.
      fit_alpha = .false.
      Call GetMem('alphaStart','Allo','Real',ipalphaStart,nvar)
      Call GetMem('alphaStop','Allo','Real',ipalphaStop,nvar)
      Call GetMem('alpha','Allo','Real',ipalpha,nvar)
      Call GetMem('alpha0','Allo','Real',ipalpha0,nvar)
      Call GetMem('alphaTemp','Allo','Real',ipalphaTemp,nvar)

      Call GetMem('alphaRad','Allo','Real',ipalphaRad,nvar)
      call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStart),1)
      call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStop),1)
!       alphaStart =-1.0d0
!       alphaStop  =-1.0d0
      Do ivar = 1,nvar
      trfcode = trfName(ivar)(1:32)
      ifit = index(trfcode,'FIT')
      If ( ifit.gt.0 ) Then
      fit_alpha = .true.
      Inline = trfName(ivar)(1:32)
      istart = index(Inline,'ALPHA=')
      istart = istart+6
      Inline = trfCode(istart:len(trfCode))
      istop = len(Inline)
      Read(Inline(1:istop),*) Work(ipalphaStart+ivar-1),                &
     &                               Work(ipalphaStop+ivar-1)
      Work(ipalphaRad+ivar-1) = (Work(ipalphaStop+ivar-1)-              &
     &                 Work(ipalphaStart+ivar-1))/2.0d0
      Work(ipalpha0+ivar-1) = Work(ipalphaStart+ivar-1)+                &
     &                 Work(ipalphaRad+ivar-1)
      End If
      End Do
!!
      If ( fit_alpha ) Then
!!
!!---- Determine which alpha radius is the largest.
      alphaMax = Work(ipalphaRad)
      ivarMax = 1
      Do ivar = 2,nvar
      If ( Work(ipalphaRad+ivar-1).gt.alphaMax ) Then
      alphaMax = Work(ipalphaRad+ivar-1)
      ivarMax = ivar
      End If
      End Do
!!
!!---- Initialize random number calculation.
      nPoints = 20000
      Call RmarIn(1802,9373)
!!
!!---- Start with the (nRandom+1)'th set of random numbers.
!!    nRandom = 4
!!    Do i = 1,nRandom ; Call Ranmar(temp,2) ; End Do
!!
!!---- Store nPoints random numbers in vector RandVec.
      length = nvar*nPoints
      Call GetMem('RandVec','Allo','Real',ipRandVec,length)

      Call Ranmar(Work(ipRandVec),length)
!!
!!---- Simulated annealing.
      stand_dev_best = 100.0d0
      iRandNum = 1
      Do While ( Work(ipalphaRad+ivarMax-1).gt.1.0d-1 )
      do iv=1,nvar
      Work(ipalphaTemp+iv-1) = Work(ipalpha0+iv-1)
      enddo
      Do iPoints = 1,1000
      Do ivar = 1,nvar
      If ( Work(ipalphaStart+ivar-1).gt.0.0d0 ) Then
      Work(ipalpha+ivar-1) = Work(ipalpha0+ivar-1)+                     &
     &                           Work(ipalphaRad+ivar-1)*(2.0d0*        &
     &     Work(ipRandVec+iRandNum+ivar-1-1)-1.0d0)
      End If
      End Do
      iRandNum = iRandNum+nvar
      Call var_to_qvar(var,Work(ipqvar),Work(ipref),                    &
     &                      Work(ipqref),Work(ipalpha),trfName,         &
     &                            ndata,nvar)
      Call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,                     &
     &          coef,nterm,                                             &
     &          stand_dev,                                              &
     &          max_err,diff_vec,use_weight)
      If ( stand_dev.lt.stand_dev_best ) Then
      stand_dev_best = stand_dev
      do iv=1,nvar
      Work(ipalphaTemp+iv-1) = Work(ipalpha+iv-1)
      enddo
      End If
      End Do
      do iv=1,nvar
      Work(ipalpha0+iv-1) = Work(ipalphaTemp+iv-1)
      enddo
      Write(6,*) (Work(ipalpha0+iv-1),iv=1,nvar),                       &
     &       stand_dev_best
      do iv=1,nvar
      Work(ipalphaRad+iv-1) = Work(ipalphaRad+iv-1)/2.0d0
      enddo
      End Do
      do iv=1,nvar
      Work(ipalpha+iv-1) = Work(ipalpha0+iv-1)
      enddo
      Write(6,*) Work(ipalpha),Work(ipalpha+1)
      Call GetMem('RandVec','Free','Real',ipRandVec,length)
!!
      End If
!!
      Call GetMem('alpha0','Free','Real',ipalpha0,nvar)
      Call GetMem('alphaTemp','Free','Real',ipalphaTemp,nvar)

      Call GetMem('alphaRad','Free','Real',ipalphaRad,nvar)
      Call GetMem('alphaStart','Free','Real',ipalphaStart,nvar)
      Call GetMem('alphaStop','Free','Real',ipalphaStop,nvar)
!!
!!---- Transform coordinates to get better numerical fit.
      Call var_to_qvar(var,Work(ipqvar),Work(ipref),Work(ipqref),       &
     &       Work(ipalpha),trfName,ndata,nvar)
!!
!!---- Fit polynomial to energies.
      Call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,coef,nterm,          &
     &  stand_dev,max_err,                                              &
     &       diff_vec,use_weight)
!!
!! Write(6,*) ndata

!! Do i = 1,ndata
!!   Write(6,'(4f15.8,es15.6)') (var(i,j),j=1,nvar),yin(i),diff_vec(i)
!! End Do
!!
      If ( find_minimum ) Then
      Call Optimize(ipow,Work(ipqvar),coef,x,energy,Hess,               &
     &    nterm,nvar, ndata)
      Else
      Call x_to_qvar(x,Work(ipref),Work(ipqref),Work(ipalpha),          &
     &    trfName,nvar)
      End If
!!
!!---- Calculate gradient and quadratic, cubic and quartic force constants.
      call funcval(x,coef,ipow,energy,nterm,nvar)
      call gradient(x,coef,ipow,grad,nterm,nvar)
      call Hessian(x,coef,ipow,Hess,nterm,nvar)
      If ( max_term.gt.2 ) Then
      call thirdDer(x,coef,ipow,D3,nterm,nvar)
      Else
      call dcopy_(l_D*l_D*l_D,[0.0d0],0,D3,1)
!              D3 = 0.0d0
      End If
      If ( max_term.gt.3 ) Then
      call fourthDer(x,coef,ipow,D4,nterm,nvar)
      Else
!              D4 = 0.0d0
      call dcopy_(l_D*l_D*l_D*l_D,[0.0d0],0,D4,1)
      End If
!!
!!---- Transform back to original coordinates.
      Call qvar_to_var(var,x,grad,Hess,D3,D4,Work(ipref),               &
     &       Work(ipqref),trfName,                                      &
     &       Work(ipalpha),max_term,ndata,nvar)
!!
      Call GetMem('alpha','Free','Real',ipalpha,nvar)
      Call GetMem('qvar','Free','Real',ipqvar,ndata*nvar)
      Call GetMem('ref','Free','Real',ipref,nvar)
      Call GetMem('qref','Free','Real',ipqref,nvar)
!!
      End
