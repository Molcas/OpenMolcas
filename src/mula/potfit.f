************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module OptMod
C!
C!  Contains:
C!    PotFit       (ipow,var,yin,coef,x,energy,Hess,D3,D4,trfName,
C!                  stand_dev,max_err,find_minimum)
C!    Optimize     (ipow,var,coef,x,energy,Hess,D3,D4)
C!    var_to_qvar  (var,qvar,ref,qref,alpha,trfName)
C!    qvar_to_var  (var,x,grad,Hess,D3,D4,ref,qref,trfName,alpha)
C!    funcval      (x,coef,ipow) Result(sum)
C!    gradient     (x,coef,ipow) Result(grad)
C!    Hessian      (x,coef,ipow) Result(Hess)
C!    thirdDer     (x,coef,ipow) Result(D3)
C!    fourthDer    (x,coef,ipow) Result(D4)
C!    ShiftHess    (Hess,shift)
C!
C!  Uses:
C!    Constants
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
C!-----------------------------------------------------------------------!
C!
c       Contains

C!-----------------------------------------------------------------------!
C!
      Subroutine PotFit(nterm,nvar,ndata,ipow,var,yin,
     &       coef,x,nOsc,energy,
     &       grad,Hess,D3,D4,trfName,stand_dev,max_err,
     &       find_minimum,max_term,use_weight,
     &       l_Hess_1,l_Hess_2,l_D)
C!
C!  Purpose:
C!    Fit polynomial to data and then depending upun the variable
C!    find_minimum, we either calculate the energy, gradient, Hessian,
C!    and qubic and quartic force constants in x or we optimize the
C!    structure before we perform the calculations.
C!
C!  Input:
C!    ipow         : Two dimensional integer array - terms in polynomial.
C!    var          : Real*8 two dimensional array - coordinates
C!                   for which we know the energy.
C!    yin          : Real*8 array - energy in points given in var.
C!    trfName      : Character array - transformation associated with each
C!                   internal coordinate.
C!
C!  Output:
C!    coef         : Real*8 two dimensional array - coefficients
C!                   for each term specified by ipow.
C!    x            : Real*8 array - geometry in internal coordinates.
C!    energy       : Real*8 - energy in minimum.
C!    Hess         : Real*8 two dimensional array - Hessian in
C!                   minimum.
C!    D3           : Real*8 three dimensional array - cubic force
C!                   constants.
C!    D4           : Real*8 four dimensional array - quartic force
C!                   constants.
C!    stand_dev    : Real*8 variable - standard deviation of fitted
C!                   values in Subroutine PolFit.
C!    max_err      : Real*8 variable - maximum error of fitted
C!                   values in Subroutine PolFit.
C!    find_minimum : Logical variable - whether or not to optimize structure.
C!    max_term     : Integer - highest power if term in polynomial fit.
C!    use_weight   : Logical variable - whether or not to use weights to
C!                   make the minimum more important in the fit.
C!
C!  Uses:
C!    LinAlg
C!    VibMod
C!    RandomMod
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
c       Use LinAlg
c       Use VibMod
C! Use RandomMod
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

C!
C!---- Initialize.
      do iv=1,nterm
      coef(iv,1) = 0.0d0
      enddo
      Call GetMem('qvar','Allo','Real',ipqvar,ndata*nvar)
      Call GetMem('ref','Allo','Real',ipref,nvar)
      Call GetMem('qref','Allo','Real',ipqref,nvar)

C!
C!---- Check if a fit of alpha is wanted.
      fit_alpha = .false.
      Call GetMem('alphaStart','Allo','Real',ipalphaStart,nvar)
      Call GetMem('alphaStop','Allo','Real',ipalphaStop,nvar)
      Call GetMem('alpha','Allo','Real',ipalpha,nvar)
      Call GetMem('alpha0','Allo','Real',ipalpha0,nvar)
      Call GetMem('alphaTemp','Allo','Real',ipalphaTemp,nvar)

      Call GetMem('alphaRad','Allo','Real',ipalphaRad,nvar)
      call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStart),1)
      call dcopy_(nvar,[-1.0d0],0,Work(ipalphaStop),1)
c       alphaStart =-1.0d0
c       alphaStop  =-1.0d0
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
      Read(Inline(1:istop),*) Work(ipalphaStart+ivar-1),
     &                               Work(ipalphaStop+ivar-1)
      Work(ipalphaRad+ivar-1) = (Work(ipalphaStop+ivar-1)-
     &                 Work(ipalphaStart+ivar-1))/2.0d0
      Work(ipalpha0+ivar-1) = Work(ipalphaStart+ivar-1)+
     &                 Work(ipalphaRad+ivar-1)
      End If
      End Do
C!
      If ( fit_alpha ) Then
C!
C!---- Determine which alpha radius is the largest.
      alphaMax = Work(ipalphaRad)
      ivarMax = 1
      Do ivar = 2,nvar
      If ( Work(ipalphaRad+ivar-1).gt.alphaMax ) Then
      alphaMax = Work(ipalphaRad+ivar-1)
      ivarMax = ivar
      End If
      End Do
C!
C!---- Initialize random number calculation.
      nPoints = 20000
      Call RmarIn(1802,9373)
C!
C!---- Start with the (nRandom+1)'th set of random numbers.
C!    nRandom = 4
C!    Do i = 1,nRandom ; Call Ranmar(temp,2) ; End Do
C!
C!---- Store nPoints random numbers in vector RandVec.
      length = nvar*nPoints
      Call GetMem('RandVec','Allo','Real',ipRandVec,length)

      Call Ranmar(Work(ipRandVec),length)
C!
C!---- Simulated annealing.
      stand_dev_best = 100.0d0
      iRandNum = 1
      Do While ( Work(ipalphaRad+ivarMax-1).gt.1.0d-1 )
      do iv=1,nvar
      Work(ipalphaTemp+iv-1) = Work(ipalpha0+iv-1)
      enddo
      Do iPoints = 1,1000
      Do ivar = 1,nvar
      If ( Work(ipalphaStart+ivar-1).gt.0.0d0 ) Then
      Work(ipalpha+ivar-1) = Work(ipalpha0+ivar-1)+
     &                           Work(ipalphaRad+ivar-1)*(2.0d0*
     &     Work(ipRandVec+iRandNum+ivar-1-1)-1.0d0)
      End If
      End Do
      iRandNum = iRandNum+nvar
      Call var_to_qvar(var,Work(ipqvar),Work(ipref),
     &                      Work(ipqref),Work(ipalpha),trfName,
     &                            ndata,nvar)
      Call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,
     &          coef,nterm,
     &          stand_dev,
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
      Write(6,*) (Work(ipalpha0+iv-1),iv=1,nvar),
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
C!
      End If
C!
      Call GetMem('alpha0','Free','Real',ipalpha0,nvar)
      Call GetMem('alphaTemp','Free','Real',ipalphaTemp,nvar)

      Call GetMem('alphaRad','Free','Real',ipalphaRad,nvar)
      Call GetMem('alphaStart','Free','Real',ipalphaStart,nvar)
      Call GetMem('alphaStop','Free','Real',ipalphaStop,nvar)
C!
C!---- Transform coordinates to get better numerical fit.
      Call var_to_qvar(var,Work(ipqvar),Work(ipref),Work(ipqref),
     &       Work(ipalpha),trfName,ndata,nvar)
C!
C!---- Fit polynomial to energies.
      Call PolFit(ipow,nvar,Work(ipqvar),yin,ndata,coef,nterm,
     &  stand_dev,max_err,
     &       diff_vec,use_weight)
C!
C! Write(6,*) ndata

C! Do i = 1,ndata
C!   Write(6,'(4f15.8,es15.6)') (var(i,j),j=1,nvar),yin(i),diff_vec(i)
C! End Do
C!
      If ( find_minimum ) Then
      Call Optimize(ipow,Work(ipqvar),coef,x,energy,Hess,
     &    nterm,nvar, ndata)
      Else
      Call x_to_qvar(x,Work(ipref),Work(ipqref),Work(ipalpha),
     &    trfName,nvar)
      End If
C!
C!---- Calculate gradient and quadratic, cubic and quartic force constants.
      call funcval(x,coef,ipow,energy,nterm,nvar)
      call gradient(x,coef,ipow,grad,nterm,nvar)
      call Hessian(x,coef,ipow,Hess,nterm,nvar)
      If ( max_term.gt.2 ) Then
      call thirdDer(x,coef,ipow,D3,nterm,nvar)
      Else
      call dcopy_(l_D*l_D*l_D,[0.0d0],0,D3,1)
c              D3 = 0.0d0
      End If
      If ( max_term.gt.3 ) Then
      call fourthDer(x,coef,ipow,D4,nterm,nvar)
      Else
c              D4 = 0.0d0
      call dcopy_(l_D*l_D*l_D*l_D,[0.0d0],0,D4,1)
      End If
C!
C!---- Transform back to original coordinates.
      Call qvar_to_var(var,x,grad,Hess,D3,D4,Work(ipref),
     &       Work(ipqref),trfName,
     &       Work(ipalpha),max_term,ndata,nvar)
C!
      Call GetMem('alpha','Free','Real',ipalpha,nvar)
      Call GetMem('qvar','Free','Real',ipqvar,ndata*nvar)
      Call GetMem('ref','Free','Real',ipref,nvar)
      Call GetMem('qref','Free','Real',ipqref,nvar)
C!
      End
