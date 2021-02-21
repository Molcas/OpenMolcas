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
      Subroutine Optimize(ipow,var,coef,x,energy,Hess,
     &  nterm, nvar, ndata  )
C!
C!  Purpose:
C!    Optimize structure given the polynomial fit from PolFit.
C!
C!  Input:
C!    ipow      : Two dimensional integer array - terms in polynomial.
C!    var       : Real*8 two dimensional array - coordinates
C!                for which we know the energy.
C!    coef      : Real*8 two dimensional array - coefficients
C!                for each term specified by ipow.
C!
C!  Output:
C!    x         : Real*8 array - geometry in minimum.
C!    energy    : Real*8 - energy in minimum.
C!    Hess      : Real*8 two dimensional array - Hessian in
C!                minimum.
C!
C!  Calls:
C!    Dool_MULA
C!    funcval
C!    gradient
C!    Hessian
C!    thirdDer
C!    fourthDer
C!    ShiftHess
C!
C!  Uses:
C!    LinAlg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
c        Use Linalg
      Implicit Real*8 ( a-h,o-z )
      Parameter ( maxiter = 100 )
      Parameter (delta_max = 1.0d0)
      Integer ipow (nterm,nvar)
      Real*8 var (ndata,nvar)
      Real*8 coef (nterm,1)
      Real*8 x (nvar)
      Real*8 xmin (nvar)
      Real*8 grad (nvar)
      Real*8 delta (nvar,1)
      Real*8 Hess (nterm,nvar)
      Real*8 var_intervals (nvar,2)
      Logical  shift
C!
C!---- Initialize.
C!
C!---- Determine intervals for random starting point in optimization.
      do iv=1,nvar
      var_intervals(iv,1) = var(1,iv)
      var_intervals(iv,2) = var(1,iv)
      enddo
      Do ivar = 1,nvar
      Do iterm = 2,nterm
      If ( var(iterm,ivar).lt.var_intervals(ivar,1) ) Then
      var_intervals(ivar,1) = var(iterm,ivar)
      Else If ( var(iterm,ivar).gt.var_intervals(ivar,2) ) Then
      var_intervals(ivar,2) = var(iterm,ivar)
      End If
      End Do
      End Do
C!
C!---- Optimize using Newton-Raphson.
      iseed=12345
      energy = 1000.0d0
      !Call Srand(0.5)
      Do i = 1,1
      Do j = 1,nvar
*             rand_number = rand()
* Random_molcas takes an integer seed 0<iseed<2.0d0**31
* and produces a pseudorandom number 0<z<1.0d0
      rand_number = Random_molcas(iseed)
      iseed=nint(rand_number*2.0d0**31)
      x(j) = var_intervals(j,1)+
     &       (var_intervals(j,2)-var_intervals(j,1))*rand_number
      End Do
      call gradient(x,coef,ipow,grad,nterm,nvar)
      Call Hessian(x,coef,ipow,Hess,nterm,nvar)
      Call ShiftHess(Hess,shift,nterm,nvar)
      do iv=1,nvar
      delta(iv,1) = -grad(iv)
      enddo
      Call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
c          Call calcNorm(delta(:,1),delta_norm)
      sum=0.0d0
      do iv=1,nvar
      sum=sum+delta(iv,1)**2
      enddo
      delta_norm=sqrt(sum)

      scale = 1.0d0
      If ( delta_norm.gt.delta_max ) scale =
     &        delta_max/delta_norm
      do iv=1,nvar
      x(iv) = x(iv)+delta(iv,1)*scale
      enddo
      iter = 0
      Do While (( delta_norm.gt.1.0d-12 ).and.( iter.le.maxiter ))
      iter = iter+1
      call gradient(x,coef,ipow,grad,nterm,nvar)
      Call Hessian(x,coef,ipow,Hess,nterm,nvar)
      Call ShiftHess(Hess,shift,nterm,nvar)
      do iv=1,nvar
      delta(iv,1) = -grad(iv)
      enddo
      Call Dool_MULA(Hess,nterm,nvar,delta,nvar,1,det)
c             Call calcNorm(delta(:,1), delta_norm)
      sum=0.0d0
      do iv=1,nvar
      sum=sum+delta(iv,1)**2
      enddo
      delta_norm=sqrt(sum)

      scale = 1.0d0
      If ( delta_norm.gt.delta_max ) scale =
     &        delta_max/delta_norm
      do iv=1,nvar
      x(iv) = x(iv)+delta(iv,1)*scale
      enddo
      End Do
      If ( iter.ge.maxiter ) Write(6,*)
     &    'WARNING!! No convergence in Optimize'
      call funcval(x,coef,ipow,fval,nterm,nvar)
      If ( fval.lt.energy ) Then
      energy = fval
      do iv=1,nvar
      xmin(iv) = x(iv)
      enddo
      End If
      End Do
      do iv=1,nvar
      x(iv) = xmin(iv)
      enddo
C!
      End
