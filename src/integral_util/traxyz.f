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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine Traxyz(nZeta,la,WInt,Scr,A)
************************************************************************
*                                                                      *
* Object: to transform the well-integrals from the local coordinate    *
*         system to the global one. The transformation matrix is       *
*         stored in A.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*              October '92.                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 WInt(nZeta*3**(la-1),3), Scr(nZeta*3**la),
     &       A(nZeta,3,3)
*
      iRout = 233
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' Enter Traxyz: WInt',' ',Wint,nZeta,3**la)
         Call RecPrt(' The transformation matrix',' ',A,nZeta,9)
      End If
*
*-----Transform
*
      nLen = 3
      mLen = 3**(la-1)
      kLen = 3**la
*     Write (*,*) ' nLen, mLen=',nLen, mLen
*-----Loop over all angular indices to be transformed.
      Do 210 i = 1, la
         iVec = 0
         Do 200 iLen = 1, mLen
            iOff = (iLen-1)*nLen
            Do 220 iZeta = 1, nZeta
               iVec = iVec + 1
*
               jVec = iOff*nZeta + iZeta
               Scr(jVec) = A(iZeta,1,1) * WInt(iVec,1) +
     &                     A(iZeta,1,2) * WInt(iVec,2) +
     &                     A(iZeta,1,3) * WInt(iVec,3)
               jVec = (iOff+1)*nZeta + iZeta
               Scr(jVec) = A(iZeta,2,1) * WInt(iVec,1) +
     &                     A(iZeta,2,2) * WInt(iVec,2) +
     &                     A(iZeta,2,3) * WInt(iVec,3)
               jVec = (iOff+2)*nZeta + iZeta
               Scr(jVec) = A(iZeta,3,1) * WInt(iVec,1) +
     &                     A(iZeta,3,2) * WInt(iVec,2) +
     &                     A(iZeta,3,3) * WInt(iVec,3)
*
 220        Continue
 200     Continue
*        Call RecPrt(' Partially transformed WInt',' ',Scr,nZeta,kLen)
         call dcopy_(nZeta*kLen,Scr,1,WInt,1)
 210  Continue
*
      If (iPrint.ge.99)
     &   Call RecPrt('Exit Traxyz :Global well integrals',' ',
     &                WInt,nZeta,kLen)
*     Call GetMem('Traxyz','Check','Real',iDum,iDum)
      Return
      End
