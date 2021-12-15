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
      SubRoutine Change(iBas, iBsInc,QiBas, kBas, kBsInc,QkBas,
     &                  jBas, jBsInc,QjBas, lBas, lBsInc,QlBas,
     &                  jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim,MaxReq,
     &                  Fail)
************************************************************************
*                                                                      *
* Object: to change the length of the primitive and basis functions    *
*         vectors of center D and B, and center  D, B, C, and A,       *
*         respectively.                                                *
*                                                                      *
* Called from: PSOAOx                                                  *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
*
      Logical QiBas, QjBas, QkBas, QlBas, QjPrim, QlPrim, Fail
*
c     iRout = 11
c     iPrint = nPrint(iRout)
      Fail = .False.
      If (QlPrim) Then
         If (lPrInc.eq.1) Then
            QlPrim=.False.
            QjPrim=.True.
            Go To 100
         End If
         Do 10 i = 2, lPrim
            If ((lPrim+1)/i.lt.lPrInc) Then
               lPrInc=Max(1,(lPrim+1)/i)
               Return
            End If
 10      Continue
      End If
 100  Continue
      If (QjPrim) Then
         lPrInc = lPrim
         If (jPrInc.eq.1) Then
            QjPrim=.False.
            QlBas =.True.
            Go To 200
         End If
         Do 20 i = 2, jPrim
            If ((jPrim+1)/i.lt.jPrInc) Then
               jPrInc=Max(1,(jPrim+1)/i)
               Return
            End If
 20      Continue
      End If
 200  Continue
      lPrInc = lPrim
      jPrInc = jPrim
      If (QlBas) Then
         If (lBsInc.eq.1) Then
            QlBas =.False.
            QjBas =.True.
            Go To 300
         End If
         Do 30 i = 2, lBas
            If ((lBas+1)/i.lt.lBsInc) Then
               lBsInc=Max(1,(lBas+1)/i)
               QlPrim=.True.
               Return
            End If
 30      Continue
      End If
 300  Continue
      If (QjBas) Then
         lBsInc = lBas
         If (jBsInc.eq.1) Then
            QjBas =.False.
            QkBas =.True.
            Go To 400
         End If
         Do 40 i = 2, jBas
            If ((jBas+1)/i.lt.jBsInc) Then
               jBsInc=Max((jBas+1)/i,1)
               QlPrim=.True.
               Return
            End If
 40      Continue
      End If
 400  Continue
      If (QkBas) Then
         lBsInc = lBas
         jBsInc = jBas
         If (kBsInc.eq.1) Then
            QkBas =.False.
            QiBas =.True.
            Go To 500
         End If
         Do 50 i = 2, kBas
            If ((kBas+1)/i.lt.kBsInc) Then
               kBsInc=Max((kBas+1)/i,1)
               QlPrim=.True.
               Return
            End If
 50      Continue
      End If
 500  Continue
      If (QiBas) Then
         lBsInc = lBas
         jBsInc = jBas
         kBsInc = kBas
         If (iBsInc.eq.1) Then
            Fail = .True.
            Return
         End If
         Do 60 i = 2, iBas
            If ((iBas+1)/i.lt.iBsInc) Then
               iBsInc=Max(1,(iBas+1)/i)
               QlPrim=.True.
               Return
            End If
 60      Continue
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(MaxReq)
      End
