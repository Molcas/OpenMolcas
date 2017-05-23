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
      SubRoutine Rys11(Arg,nArg,Root,Weight,iPntr,nPntr,
     &                 x0,nMax,R6,R5,R4,R3,R2,R1,R0,
     &                         W6,W5,W4,W3,W2,W1,W0,ddx,
     &                 HerW,HerR2,TMax)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Arg(nArg), Root(nArg), Weight(nArg), x0(nMax),
     &       R6(nMax), R5(nMax),
     &       R4(nMax), R3(nMax), R2(nMax), R1(nMax), R0(nMax),
     &       W6(nMax), W5(nMax),
     &       W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax)
      Integer iPntr(nPntr)
*
      xdInv=One/ddx
      dddx=ddx/10d0 + ddx
      Do iArg = 1, nArg
         If (Arg(iArg).lt.TMax) Then
            n = iPntr(Int((Arg(iArg)+dddx)*xdInv))
            z = Arg(iArg) - x0(n)
            Root(iArg)=(((((r6(n)*z+r5(n))*z+r4(n))*z+r3(n))*z
     &       +r2(n))*z+r1(n))*z+r0(n)
            Weight(iArg)=(((((w6(n)*z+w5(n))*z+w4(n))*z+w3(n))*z
     &       +w2(n))*z+w1(n))*z+w0(n)
         Else
            ai=1.0d0/Arg(iArg)
            Root(iArg) = HerR2*ai
            Weight(iArg) = HerW*Sqrt(ai)
         End If
      End Do
      Return
      End
