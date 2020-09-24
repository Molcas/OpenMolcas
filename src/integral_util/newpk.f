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
* Copyright (C) 1990-1992,1999, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine NewPK(A,B,P,mZeta,nZeta,rKappa,Alpha,Beta)
************************************************************************
* Object : to compute P and kappa.                                     *
*                                                                      *
* Called from: k2Loop                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*             May '90, modified for integral cutoff.                   *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN.                              *
*             June '91, modified for k2 loop.                          *
*             January '92, modified for gradient calculations.         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Alpha(nZeta), Beta(nZeta), A(3), B(3), P(nZeta,3),
     &       rKappa(nZeta)
*
      iRout = 243
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In NewPK:Alpha',' ',Alpha,mZeta,1)
         Call RecPrt(' In NewPK:Beta',' ',Beta,mZeta,1)
      End If
*
      AB2=(A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
      Do iZeta = 1, mZeta
         Tmp0=One/(Alpha(iZeta)+Beta(iZeta))
         Tmp1=TwoP54*Exp(-Alpha(iZeta)*Beta(iZeta)*AB2*Tmp0)*Tmp0
         If (Tmp1.lt.1.0D-99) Tmp1=1.0D-99
         rKappa(iZeta) = Tmp1
         P(iZeta,1) = (Alpha(iZeta)*A(1)+Beta(iZeta)*B(1))*Tmp0
         P(iZeta,2) = (Alpha(iZeta)*A(2)+Beta(iZeta)*B(2))*Tmp0
         P(iZeta,3) = (Alpha(iZeta)*A(3)+Beta(iZeta)*B(3))*Tmp0
      End Do
      Do iZeta = mZeta+1, nZeta
         rKappa(iZeta) = Zero
         P(iZeta,1) = Zero
         P(iZeta,2) = Zero
         P(iZeta,3) = Zero
      End Do
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In NewPK: Kappa',' ',rKappa,mZeta,1)
         Call RecPrt(' In NewPK: Px',' ',P(1,1),mZeta,1)
         Call RecPrt(' In NewPK: Py',' ',P(1,2),mZeta,1)
         Call RecPrt(' In NewPK: Px',' ',P(1,3),mZeta,1)
      End If
*
      Return
      End
