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
* Copyright (C) 1992,1993, Roland Lindh                                *
************************************************************************
      SubRoutine Screen(nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                  Zeta,ZInv,P,KappAB,IndZet,Data1,nAlpha,nBeta,
     &                  IndZ,ZtMax,abMax,ZtMaxD,abMaxD,
     &                  Eta,EInv,Q,KappCD,IndEta,Data2,nGamma,nDelta,
     &                  IndE,EtMax,cdMax,EtMaxD,cdMaxD,
     &                  Dij,Dkl,
     &                  iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutDInt,
     &                  CutInt,vij,vkl,vik,vil,vjk,vjl,
     &                  Prescreen_On_Int_Only)
************************************************************************
*                                                                      *
* Object: to prescreen the integrals for Direct SCF                    *
*                                                                      *
*   nZeta, nEta : unpartioned length of primitives.                    *
*                                                                      *
*   mZeta, mEta : section length due to partioning. These are usually  *
*                 equal to nZeta and nEta.                             *
*                                                                      *
*   lZeta, lEta : section length after prescreening.                   *
*                                                                      *
*   Observe that the integrals are ordered for optimal performance in  *
*   the integral generation step whereas the 1st order density is      *
*   canonically ordered.                                               *
*                                                                      *
* Called from: Twoel                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             March '92                                                *
*             April '92 modified for gradient estimate                 *
*             January '93 modified to Direct SCF                       *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
      Real*8 Zeta(mZeta), ZInv(mZeta), KappAB(mZeta), P(nZeta,3),
     &        Eta(mEta),  EInv(mEta),  KappCD(mEta),  Q(nEta, 3),
     &       Data1(nZeta*(nDArray-1)), Data2(nEta*(nDArray-1)),
     &       Dij(nZeta), Dkl(nEta)
      Real*8 ZtMax,EtMax,abMax,cdMax,ZtMaxD,EtMaxD,abMaxD,cdMaxD
      Integer   IndZet(nZeta), IndEta(nEta),
     &          IndZ(nZeta), IndE(nEta)
      Logical Prescreen_On_Int_Only
*
#include "print.fh"
#include "real.fh"
*
*     decalaration of local variables...
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt(' In Screen: Data1',' ',Data1,nZeta,nDArray-1)
      Call RecPrt(' In Screen: Data2',' ',Data2,nEta ,nDArray-1)
      Call FZero(P,nZeta*3)
      Call FZero(Q,nEta *3)
#endif
*
      If (.Not.Prescreen_On_Int_Only) Then
*
*        The former technique will be used for integrals which will
*        be partially stored.
*
*-----   Set up prescreening parameters
*
         DMax = Max(vik,vil,vjk,vjl)/Four            ! AO basis
*
*-----   Combine threshold with maximum of density contribution.
*        Note that we are using the square to avoid doing the
*        square root in the denominator!
         Cut  = Max(CutDInt/Max(DMax,vij,vkl),CutDInt)      ! AO basis
*
*-----   Test on the largest integral which will be produced by this
*        batch.
*
         vMax = abMax*cdMax                     ! integral in AO basis
*-----   Skip prescreening if too low.
         If (vMax.lt.Cut) Then
            lZeta = 0
            lEta  = 0
            Go To 999
         End If
      Else
         DMax=Zero
      End If
*
      lZeta=0
      lEta =0
*
*---- prescreen zeta pairs
      If (.Not.Prescreen_On_Int_Only) Then
         lZeta=0
         Do iZeta = 1, mZeta
            ppaa= Data1(ip_ab   (iZeta,nZeta)) * cdMax
            aaaa= Data1(ip_abCon(iZeta,nZeta)) * cdMax
            jZeta = IndZ(iZeta)
            Test=ppaa*Dij(jZeta)+aaaa*(DMax+vkl)
            If (Test.ge.CutDInt) Then
               lZeta=lZeta+1
               Zeta(lZeta)  = Data1(ip_Z    (iZeta,nZeta))
               KappAB(lZeta)= Data1(ip_Kappa(iZeta,nZeta))
               P(lZeta,1)   = Data1(ip_Pcoor(iZeta,        nZeta))
               P(lZeta,2)   = Data1(ip_Pcoor(iZeta+  nZeta,nZeta))
               P(lZeta,3)   = Data1(ip_Pcoor(iZeta+2*nZeta,nZeta))
               IndZet(lZeta)= IndZ(iZeta)
               ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
            End If
         End Do
      Else
         lZeta=0
         Do iZeta = 1, mZeta
            aaaa= Data1(ip_abCon(iZeta,nZeta)) * cdMax
            If (aaaa.ge.CutInt) Then
               lZeta=lZeta+1
               Zeta(lZeta)  = Data1(ip_Z    (iZeta,nZeta))
               KappAB(lZeta)= Data1(ip_Kappa(iZeta,nZeta))
               P(lZeta,1)   = Data1(ip_Pcoor(iZeta,        nZeta))
               P(lZeta,2)   = Data1(ip_Pcoor(iZeta+  nZeta,nZeta))
               P(lZeta,3)   = Data1(ip_Pcoor(iZeta+2*nZeta,nZeta))
               IndZet(lZeta)= IndZ(iZeta)
               ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
            End If
         End Do
      End If
      If (lZeta.eq.0) Go To 999
*
      If (iphX1.ne.1) Call DScal_(lZeta,-One,P(1,1),1)
      If (iphY1.ne.1) Call DScal_(lZeta,-One,P(1,2),1)
      If (iphZ1.ne.1) Call DScal_(lZeta,-One,P(1,3),1)
*
*---- prescreen eta pairs
      If (.Not.Prescreen_On_Int_Only) Then
         lEta=0
         Do iEta = 1, mEta
            ppaa= Data2(ip_ab   (iEta,nEta)) * abMax
            aaaa= Data2(ip_abCon(iEta,nEta)) * abMax
            jEta = IndE(iEta)
            Test=ppaa*Dkl(jEta)+aaaa*(DMax+vij)
            If (Test.ge.CutDInt) Then
               lEta=lEta+1
               IndEta(lEta)= IndE(iEta)
               Eta(lEta)   = Data2(ip_Z    (iEta,nEta))
               KappCD(lEta)= Data2(ip_Kappa(iEta,nEta))
               Q(lEta,1)   = Data2(ip_Pcoor(iEta,       nEta))
               Q(lEta,2)   = Data2(ip_Pcoor(iEta+  nEta,nEta))
               Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
               EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
            End If
         End Do
      Else
         lEta=0
         Do iEta = 1, mEta
            aaaa= Data2(ip_abCon(iEta,nEta)) * abMax
            If (aaaa.ge.CutInt) Then
               lEta=lEta+1
               IndEta(lEta)= IndE(iEta)
               Eta(lEta)   = Data2(ip_Z    (iEta,nEta))
               KappCD(lEta)= Data2(ip_Kappa(iEta,nEta))
               Q(lEta,1)   = Data2(ip_Pcoor(iEta,       nEta))
               Q(lEta,2)   = Data2(ip_Pcoor(iEta+  nEta,nEta))
               Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
               EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
            End If
         End Do
      End If
      If (lEta.eq.0) Go To 999
*
      If (iphX2.ne.1) Call DScal_(lEta,-One,Q(1,1),1)
      If (iphY2.ne.1) Call DScal_(lEta,-One,Q(1,2),1)
      If (iphZ2.ne.1) Call DScal_(lEta,-One,Q(1,3),1)
*
 999  Continue
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) ' In Screen'
      Call RecPrt(' Zeta  ',' ',Zeta,lZeta,1)
      Call RecPrt(' Eta   ',' ', Eta, lEta,1)
      Call RecPrt(' ZInv  ',' ',ZInv,lZeta,1)
      Call RecPrt(' EInv  ',' ',EInv,lEta, 1)
      Call RecPrt(' Px    ',' ',P(1,1),lZeta,1)
      Call RecPrt(' Py    ',' ',P(1,2),lZeta,1)
      Call RecPrt(' Pz    ',' ',P(1,3),lZeta,1)
      Call RecPrt(' Qx    ',' ',Q(1,1),lEta, 1)
      Call RecPrt(' Qy    ',' ',Q(1,2),lEta, 1)
      Call RecPrt(' Qz    ',' ',Q(1,3),lEta, 1)
      Call RecPrt(' KappAB',' ',KappAB,lZeta,1)
      Call RecPrt(' KappCD',' ',KappCD,lEta, 1)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nAlpha)
         Call Unused_integer(nBeta)
         Call Unused_real(ZtMax)
         Call Unused_real(ZtMaxD)
         Call Unused_real(abMaxD)
         Call Unused_integer(nGamma)
         Call Unused_integer(nDelta)
         Call Unused_real(EtMax)
         Call Unused_real(EtMaxD)
         Call Unused_real(cdMaxD)
      End If
      End
