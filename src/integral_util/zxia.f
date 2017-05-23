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
      SubRoutine ZXia(Zeta,ZInv,N,M,Alpha,Beta)
************************************************************************
*                                                                      *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : GetMem                                                  *
*              RecPrt                                                  *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Zeta(N,M), ZInv(N,M), Alpha(N), Beta(M)
*
      iRout = 113
      iPrint = nPrint(iRout)
*     Call qEnter('ZXia')
*     Call GetMem(' Enter ZXia','CHECK','REAL',iDum,iDum)
      If (iPrint.ge.99) Then
         Call RecPrt(' Alpha',' ',Alpha,N,1)
         Call RecPrt(' Beta ',' ',Beta ,M,1)
      End If
*
      Do 10 i = 1, N
         Do 20 j = 1, M
            Zeta(i,j) = (Alpha(i)+Beta(j))
            ZInv(i,j) = One/Zeta(i,j)
 20      Continue
 10   Continue
      If (iPrint.ge.99) Then
         Call RecPrt( ' In ZXia: Zeta',' ',Zeta,N,M)
      End If
*     Call GetMem('Exit ZXia','CHECK','REAL',iDum,iDum)
*
*     Call qExit('ZXia')
      Return
      End
