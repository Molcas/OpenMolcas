************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Function Eval_RMax(alpha,m,R_L)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Eval_RMax
*     Write (6,*) 'alpha,m,R_L=',
*    &            alpha,m,R_L
*
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute r_k_H as a function of m and R_L
*
*     Eq(19) R. Lindh, P.-A. Malmqvist, L. Gagliardi,
*     TCA, 106:178-187 (2001)

      If (MOD(m+3,2).eq.0) Then
         Gamma=One
         Do i = 2, (m+3)/2, 1
            Gamma = Gamma * DBLE(i-1)
         End Do
      Else
         Gamma=Sqrt(Pi)
         Do i = 5, m+3, 2
            Gamma = Gamma * DBLE(i-1)/Two
         End Do
      End If
*
*     x = Alpha * (r_k_H)**2
*
      x = 10.0D0 ! Start value
 123  Continue
      x_new = LOG((Gamma/R_L)*x**(Half*(DBLE(m)+One)))
C     Write (6,*) 'x,x_new=',x,x_new
      If (ABS(x-x_new).gt.1.0D-8) Then
         x=x_new
         Go To 123
      End If
*
      Eval_RMax=Sqrt(x/alpha)
C     Write (*,*) 'Eval_RMax=',Eval_RMax
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
