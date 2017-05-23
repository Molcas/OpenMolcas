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
      SubRoutine Rys01(Arg,nArg,Weight,iPntr,nPntr,
     &                 x0,nMax,W6,W5,W4,W3,W2,W1,W0,ddx,HerW,Tmax)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Arg(nArg), Weight(nArg), x0(nMax), W6(nMax), W5(nMax),
     &       W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax)
      Integer iPntr(nPntr)
*
      xdInv=One/ddx
      dddx = ddx/10d0 + ddx
      Do 10 iArg = 1, nArg
         If (Arg(iArg).lt.TMax) Then
            n = iPntr(Int((Arg(iArg)+dddx)*xdInv))
            z = Arg(iArg) - x0(n)
            w =(((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))*z+W1(n))
     &          *z+w0(n)
            Weight(iArg) = w
         Else
            ai = 1.0D0/Arg(iArg)
            Weight(iArg) = HerW*Sqrt(ai)
         End If
 10   Continue
      Return
      End
