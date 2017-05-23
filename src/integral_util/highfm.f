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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine HighFm(F,T,m,n)
************************************************************************
*  Object: to compute the auxiliary function for orders which we do    *
*          not use Shavitt's method of tables.                         *
*                                                                      *
* Called from: Auxil                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 F(n), T(n)
#include "real.fh"
*
*     Find T for which the asympotic formula can be used
*
      Tmax=50.0D0
 88   gTmp=Gamma2(m,Tmax)
         i = 1
         ii = 2*m-1
         sum1 = One
         sum0 =  One
 77      Sum1 = Sum1 * DBLE(ii) / (Two*Tmax)
         Sum0 = Sum0 + Sum1
         ii = ii - 2
         i = i + 1
         If (i.ge.m .and. Sum1/sum0.gt.1.0d-11) Go to 77
      Tnew = Log( sum0 / (2.0D-16*Tmax*gTmp)  )
      If (Abs(Tnew-Tmax).lt.1.0d-9) Go To 97
      Tmax = Tnew
      Go To 88
*
*     Compute the auxiliary functions
*
 97   Tmax = Tnew
      Do 100 k = 1, n
         If (T(k).lt.Tmax) Then
            Fvalue = Zero
            i = 0
            Term = One
 99         Continue
               Term = Term  / DBLE(2*(m+i)+1)
               Fvalue = Fvalue + Term
               i = i + 1
               Term = Term * (Two*T(k))
               If (Abs(Term/Fvalue).gt.1.0d-18) Go To 99
            F(k) = Exp(-T(k))*Fvalue
         Else
            F(k) = Gamma2(m,T(k))
         End If
 100  Continue
*
      Return
      End
