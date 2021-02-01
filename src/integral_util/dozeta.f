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
* Copyright (C) 1990-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine DoZeta(Alpha,nAlpha,Beta,nBeta,A,B,P,Zeta,rKappa,
     &                  ZInv,Alpha_,Beta_,Ind_Pair)
************************************************************************
*                                                                      *
* Object : to compute P and kappa.                                     *
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
      Real*8 Alpha(nAlpha), Beta(nBeta), Zeta(nAlpha*nBeta),
     &       Alpha_(nAlpha*nBeta), Beta_(nAlpha*nBeta),
     &       ZInv(nAlpha*nBeta), A(3), B(3),
     &       P(nAlpha*nBeta,3), rKappa(nAlpha*nBeta)
      Integer Ind_Pair(nAlpha*nBeta+1)
*
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt(' In DoZeta:Alpha',' ',Alpha,nAlpha,1)
      Call RecPrt(' In DoZeta:Beta',' ',Beta,nBeta,1)
#endif
*
      AB2=(A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
      Do iBeta = 1, nBeta
         Do iAlpha = 1, nAlpha
            iZeta = nAlpha*(iBeta-1) + iAlpha
            Zeta  (iZeta)=Alpha(iAlpha)+Beta(iBeta)
            Alpha_(iZeta)=Alpha(iAlpha)
            Beta_ (iZeta)=Beta(iBeta)
            ZInv  (iZeta)=One/Zeta(iZeta)
            Tmp0=ZInv(iZeta)
            Tmp1=TwoP54*Exp(-Alpha(iAlpha)*Beta(iBeta)*AB2*Tmp0)*Tmp0
            If (Tmp1.lt.1.0D-99) Tmp1=1.0D-99
            rKappa(iZeta) = Tmp1
            P(iZeta,1) = (Alpha(iAlpha)*A(1)+Beta(iBeta)*B(1))*Tmp0
            P(iZeta,2) = (Alpha(iAlpha)*A(2)+Beta(iBeta)*B(2))*Tmp0
            P(iZeta,3) = (Alpha(iAlpha)*A(3)+Beta(iBeta)*B(3))*Tmp0
            Ind_Pair(iZeta)=iZeta
         End Do
      End Do
      Ind_Pair(nAlpha*nBeta+1)=nAlpha*nBeta
*
*     Sort from Large to Small
*
*define _New_Code_
#if defined(_New_Code_) || defined(_DEBUGPRINT_)
      nZeta=nAlpha*nBeta
#endif
#ifdef _New_Code_
      Do iZeta = 1, nZeta-1
         Tmp1 = rKappa(iZeta)
         Do jZeta = iZeta+1, nZeta
            If (Tmp1.lt.rKappa(jZeta)) Then
               Tmp1 = rKappa(jZeta)
               Tmp2 = Zeta(iZeta)
               Zeta(iZeta)=Zeta(jZeta)
               Zeta(jZeta)=Tmp2
               Tmp2 = Alpha_(iZeta)
               Alpha_(iZeta)=Alpha_(jZeta)
               Alpha_(jZeta)=Tmp2
               Tmp2 = Beta_ (iZeta)
               Beta_ (iZeta)=Beta_ (jZeta)
               Beta_ (jZeta)=Tmp2
               Tmp2 = ZInv  (iZeta)
               ZInv  (iZeta)=ZInv  (jZeta)
               ZInv  (jZeta)=Tmp2
               Tmp2 = rKappa(iZeta)
               rKappa(iZeta)=rKappa(jZeta)
               rKappa(jZeta)=Tmp2
               Tmp2 = P(iZeta,1)
               P(iZeta,1)=P(jZeta,1)
               P(jZeta,1)=Tmp2
               Tmp2 = P(iZeta,2)
               P(iZeta,2)=P(jZeta,2)
               P(jZeta,2)=Tmp2
               Tmp2 = P(iZeta,3)
               P(iZeta,3)=P(jZeta,3)
               P(jZeta,3)=Tmp2
               iTmp = Ind_Pair(iZeta)
               Ind_Pair(iZeta)=Ind_Pair(jZeta)
               Ind_Pair(jZeta)=iTmp
            End If
         End Do
      End Do
#endif
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' In DoZeta: Kappa',' ',rKappa,nZeta,1)
      Call RecPrt(' In DoZeta: P',' ',P,nZeta,3)
#endif
*
      Return
      End
