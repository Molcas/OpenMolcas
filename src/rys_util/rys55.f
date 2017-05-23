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
      SubRoutine Rys55(Arg,nArg,Root,Weight,iPntr,nPntr,
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
      Real*8 Arg(nArg), Root(5,nArg), Weight(5,nArg), x0(nMax),
     &       R6(nMax,5), R5(nMax,5),
     &       R4(nMax,5), R3(nMax,5), R2(nMax,5), R1(nMax,5), R0(nMax,5),
     &       W6(nMax,5), W5(nMax,5),
     &       W4(nMax,5), W3(nMax,5), W2(nMax,5), W1(nMax,5), W0(nMax,5),
     &       HerW(5), HerR2(5)
      Integer iPntr(nPntr)
*
      xdInv=One/ddx
      dddx=ddx/10d0 + ddx
      Do iArg = 1, nArg
         If (ARg(iArg).lt.TMax) Then
            n = iPntr(Int((Arg(iArg)+dddx)*xdInv))
            z = Arg(iArg) - x0(n)
            r = (((((R6(n,1)*z+R5(n,1))*z+R4(n,1))*z+R3(n,1))*z+R2(n,1))
     &          *z+R1(n,1))*z+R0(n,1)
            Root(1,iArg)= r
            r = (((((R6(n,2)*z+R5(n,2))*z+R4(n,2))*z+R3(n,2))*z+R2(n,2))
     &          *z+R1(n,2))*z+R0(n,2)
            Root(2,iArg)= r
            r = (((((R6(n,3)*z+R5(n,3))*z+R4(n,3))*z+R3(n,3))*z+R2(n,3))
     &          *z+R1(n,3))*z+R0(n,3)
            Root(3,iArg)= r
            r = (((((R6(n,4)*z+R5(n,4))*z+R4(n,4))*z+R3(n,4))*z+R2(n,4))
     &          *z+R1(n,4))*z+R0(n,4)
            Root(4,iArg)= r
            r = (((((R6(n,5)*z+R5(n,5))*z+R4(n,5))*z+R3(n,5))*z+R2(n,5))
     &          *z+R1(n,5))*z+R0(n,5)
            Root(5,iArg)= r
            r = (((((W6(n,1)*z+W5(n,1))*z+W4(n,1))*z+W3(n,1))*z+W2(n,1))
     &          *z+W1(n,1))*z+W0(n,1)
            Weight(1,iArg) = r
            r = (((((W6(n,2)*z+W5(n,2))*z+W4(n,2))*z+W3(n,2))*z+W2(n,2))
     &          *z+W1(n,2))*z+W0(n,2)
            Weight(2,iArg) = r
            r = (((((W6(n,3)*z+W5(n,3))*z+W4(n,3))*z+W3(n,3))*z+W2(n,3))
     &          *z+W1(n,3))*z+W0(n,3)
            Weight(3,iArg) = r
            r = (((((W6(n,4)*z+W5(n,4))*z+W4(n,4))*z+W3(n,4))*z+W2(n,4))
     &          *z+W1(n,4))*z+W0(n,4)
            Weight(4,iArg) = r
            r = (((((W6(n,5)*z+W5(n,5))*z+W4(n,5))*z+W3(n,5))*z+W2(n,5))
     &          *z+W1(n,5))*z+W0(n,5)
            Weight(5,iArg) = r
         Else
            ai = 1.0D0/Arg(iArg)
            si = Sqrt(ai)
            Root(1,iArg) = HerR2(1)*ai
            Root(2,iArg) = HerR2(2)*ai
            Root(3,iArg) = HerR2(3)*ai
            Root(4,iArg) = HerR2(4)*ai
            Root(5,iArg) = HerR2(5)*ai
            Weight(1,iArg) = HerW(1)*si
            Weight(2,iArg) = HerW(2)*si
            Weight(3,iArg) = HerW(3)*si
            Weight(4,iArg) = HerW(4)*si
            Weight(5,iArg) = HerW(5)*si
         End If
      End Do
      Return
      End
