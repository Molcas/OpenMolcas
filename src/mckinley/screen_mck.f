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
*               1995, Anders Bernhardsson                              *
************************************************************************
#define _Old_Code_
#ifdef _Old_Code_
      SubRoutine Screen_mck(PAO,Scrtch,mPAO,
     &                  nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                  Zeta,ZInv,P,xA,xB,rKA,
     &                  Data1,IndZ,ztmx,abmax,zexpmax,
     &                  nAlpha,nBeta,
     &                  Eta, EInv,Q,xG,xD,rKC,Data2,IndE,
     &                  etmx,cdmax,eexpmax,nGamma,nDelta,
     &                  xpre,
     &                  iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutInt,
     &                  PreScr,IndZet,IndEta,ldot)
************************************************************************
*                                                                      *
* Object: to prescreen the integral derivatives.                       *
*                                                                      *
*   nZeta, nEta : unpartioned length of primitives.                    *
*                                                                      *
*   mZeta, mEta : section length due to partioning. These are usually  *
*                 equal to nZeta and nEta.                             *
*                                                                      *
*   lZeta, lEta : section length after prescreening.                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             March '92                                                *
*                                                                      *
*             April '92 modified for gradient estimate                 *
*                                                                      *
*             Anders Bernhardsson 1995                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
      Real*8 PAO(mZeta*mEta*mPAO),Scrtch(mZeta*mEta*(1+mPAO*2)),
     &       Zeta(nZeta), ZInv(nZeta),  P(nZeta,3),
     &        Eta(nEta),  EInv(nEta),   Q(nEta, 3),
     &       xA(nZeta), xB(nZeta), xG(nEta), xD(nEta),
     &       Data1(nZeta*nDArray),
     &       Data2(nEta *nDArray),
     &       rKA(nZeta),rKC(nEta),
     &       xpre(mZeta*mEta)
      Logical PreScr,ldot
      Integer IndEta(nEta),IndZet(nZeta), IndZ(mZeta), IndE(mEta)
#include "real.fh"
#ifdef _DEBUG_
#include "print.fh"
#endif
*
#ifdef _DEBUG_
      iRout = 180
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt(' Data1',' ',Data1,nZeta,nDArray)
         Call RecPrt(' Data2',' ',Data2,nEta ,nDArray)
         Call RecPrt('2nd order density matrix',' ',
     &               PAO,mZeta*mEta,mPAO)
      End If
#endif
*
      ip=1
      ipPAO = ip
      ip = ip + mZeta*mEta*mPAO
*
*-----Compress all indices except zeta
*
      ipOAP = ip
      ip = ip + mZeta*mEta*mPAO
      If (ldot) Call DGetMO(PAO,mZeta,mZeta,mEta*mPAO,
     &                      Scrtch(ipOAP),mEta*mPAO)
*
*-----Prescreen Zeta
*
      lZeta=0
      Call IZero(IndZet,nZeta)
      If (PreScr) Then
         Do iZeta = 1, mZeta
            jZeta = IndZ(iZeta)
            IndZet(jZeta) = -lZeta
            abcd= Data1(ip_ab(iZeta,nZeta)) * cdMax
            If (Abs(abcd).ge.CutInt) Then
               lZeta=lZeta+1
               IndZet(jZeta) = lZeta
               Zeta(lZeta)  = Data1(ip_Z(iZeta,nZeta))
               rKA(lZeta)   = Data1(ip_Kappa(iZeta,nZeta))
               P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))
               P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))
               P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
               xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
               xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
               ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
               ip1 = ipOAP + mEta*mPAO*(iZeta-1)
               ip2 = ipPAO + mEta*mPAO*(lZeta-1)
               If (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,
     &                                        Scrtch(ip2),1)
            End If
         End Do
      Else
         Do iZeta = 1, mZeta
            lZeta=lZeta+1
            jZeta = IndZ(iZeta)
            IndZet(jZeta) = lZeta
            Zeta(lZeta)  = Data1(ip_Z(iZeta,nZeta))
            rKA(lZeta)   = Data1(ip_Kappa(iZeta,nZeta))
            P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))
            P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))
            P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
            xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
            xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
            ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
            ip1 = ipOAP + mEta*mPAO*(iZeta-1)
            ip2 = ipPAO + mEta*mPAO*(lZeta-1)
            If (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,
     &                                     Scrtch(ip2),1)
         End Do
      End If
      If (lZeta.eq.0) Go To 999
*
      If (iphX1.ne.1) Call DScal_(lZeta,-One,P(1,1),1)
      If (iphY1.ne.1) Call DScal_(lZeta,-One,P(1,2),1)
      If (iphZ1.ne.1) Call DScal_(lZeta,-One,P(1,3),1)
*
*-----Transpose eta,mPAO,zeta to mPAO,zeta,eta
*
      If (lDot) Call DGetMO(Scrtch(ipPAO),mEta,mEta,mPAO*lZeta,
     &                      Scrtch(ipOAP),mPAO*lZeta)
*
*-----Prescreen Eta
*
      lEta=0
      Call IZero(IndEta,nEta)
      If (PreScr) Then
         Do iEta = 1, mEta
            jEta = IndE(iEta)
            IndEta(jEta) = - lEta ! To be removed
            abcd= Data2(ip_ab(iEta,nEta)) * abMax
            If (Abs(abcd).ge.CutInt) Then
               lEta=lEta+1
               IndEta(jEta) = lEta
               Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
               rKC(lEta)   = Data2(ip_Kappa(iEta,nEta))
               Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))
               Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))
               Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
               xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
               xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
               EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
               ip1 = ipOAP + mPAO*lZeta*(iEta-1)
               ip2 = ipPAO + mPAO*lZeta*(lEta-1)
               If (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,
     &                                         Scrtch(ip2),1)
            End If
         End Do
      Else
         Do iEta = 1, mEta
            lEta=lEta+1
            jEta = IndE(iEta)
            IndEta(jEta) = lEta
            Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
            rKC(lEta)   = Data2(ip_Kappa(iEta,nEta))
            Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))
            Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))
            Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
            xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
            xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
            EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
            ip1 = ipOAP + mPAO*lZeta*(iEta-1)
            ip2 = ipPAO + mPAO*lZeta*(lEta-1)
            if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,
     &                                      Scrtch(ip2),1)
         End Do
      End If
      If (lEta.eq.0) Go To 999
*
      If (iphX2.ne.1) Call DScal_(lEta,-One,Q(1,1),1)
      If (iphY2.ne.1) Call DScal_(lEta,-One,Q(1,2),1)
      If (iphZ2.ne.1) Call DScal_(lEta,-One,Q(1,3),1)
*
*-----Transpose mPAO,zeta,eta to zeta,eta,mPAO
*
      If (ldot) Call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,
     &                                          PAO,lZeta*lEta)
*
 999  Continue
*
      ij = 0
      Do iEta = 1,lEta
         Et   = Eta(iEta)
         rKCD = rkC(iEta)
         Do iZeta = 1,lZeta
            Zt   = Zeta(iZeta)
            rKAB = rkA(iZeta)
            ij = ij + 1
            xpre(ij) = rKAB*rKCD*Sqrt(1.0d0/(Zt+Et))
         End Do
      End Do
      If (ldot) Then
         jPAO = 0
         Do iPAO = 1, mPAO
            Do iZE = 0, lZeta*lEta-1
               jPAO=jPAO+1
               PAO(jPAO)=xpre(iZE+1)*PAO(jPAO)
            End Do
         End Do
      End If
#ifdef _DEBUG_
      If (iPrint.ge.39) Call RecPrt(' PAO',' ',
     &                              PAO,lZeta*lEta,mPAO)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(ztmx)
         Call Unused_real(zexpmax)
         Call Unused_integer(nAlpha)
         Call Unused_integer(nBeta)
         Call Unused_real(etmx)
         Call Unused_real(eexpmax)
         Call Unused_integer(nGamma)
         Call Unused_integer(nDelta)
      End If
      End
#else
      SubRoutine Screen_mck(PAO,Scrtch,mPAO,
     &                  nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                  Zeta,ZInv,P,xA,xB,rKA,
     &                  Data1,IndZ,ztmx,abmax,zexpmax,
     &                  nAlpha,nBeta,
     &                  Eta, EInv,Q,xG,xD,rKC,Data2,IndE,
     &                  etmx,cdmax,eexpmax,nGamma,nDelta,
     &                  xpre,
     &                  iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutInt,
     &                  PreScr,IndZet,IndEta,ldot)
************************************************************************
*                                                                      *
* Object: to prescreen the integral derivatives.                       *
*                                                                      *
*   nZeta, nEta : unpartioned length of primitives.                    *
*                                                                      *
*   mZeta, mEta : section length due to partioning. These are usually  *
*                 equal to nZeta and nEta.                             *
*                                                                      *
*   lZeta, lEta : section length after prescreening.                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             March '92                                                *
*                                                                      *
*             April '92 modified for gradient estimate                 *
*                                                                      *
*             Anders Bernhardsson 1995                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
      Real*8 PAO(mZeta*mEta*mPAO),Scrtch(mZeta*mEta*(1+mPAO*2)),
     &       Zeta(nZeta), ZInv(nZeta),  P(nZeta,3),
     &        Eta(nEta),  EInv(nEta),   Q(nEta, 3),
     &       xA(nZeta), xB(nZeta), xG(nEta), xD(nEta),
     &       Data1(nZeta*nDArray),
     &       Data2(nEta *nDArray),
     &       rKA(nZeta),rKC(nEta),
     &       xpre(mZeta*mEta)
      Logical PreScr,ldot
      Integer IndEta(nEta),IndZet(nZeta), IndZ(mZeta), IndE(mEta)
#include "real.fh"
#ifdef _DEBUG_
#include "print.fh"
#endif
*
#ifdef _DEBUG_
      iRout = 180
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt(' Data1',' ',Data1,nZeta,nDArray)
         Call RecPrt(' Data2',' ',Data2,nEta ,nDArray)
         Call RecPrt('2nd order density matrix',' ',
     &               PAO,mZeta*mEta,mPAO)
      End If
#endif
*
      ip=1
      ipPAO = ip
      ip = ip + mZeta*mEta*mPAO
*
*-----Compress all indices except zeta
*
      ipOAP = ip
      ip = ip + mZeta*mEta*mPAO
      If (ldot) Call DGetMO(PAO,mZeta,mZeta,mEta*mPAO,
     &                      Scrtch(ipOAP),mEta*mPAO)
*
*-----Prescreen Zeta
*
      lZeta=0
      If (PreScr) Then
         Do iZeta = 1, mZeta
            jZeta = IndZ(iZeta)
            abcd= Data1(ip_ab(jZeta,nZeta)) * cdMax
            If (Abs(abcd).ge.CutInt) Then
               lZeta=lZeta+1
               IndZet(lZeta) = IndZ(iZeta)
               Zeta(lZeta)  = Data1(ip_Z(iZeta,nZeta))
               rKA(lZeta)   = Data1(ip_Kappa(iZeta,nZeta))
               P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))
               P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))
               P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
               xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
               xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
               ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
               ip1 = ipOAP + mEta*mPAO*(iZeta-1)
               ip2 = ipOAP + mEta*mPAO*(lZeta-1)
               If (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,
     &                                        Scrtch(ip2),1)
            End If
         End Do
      Else
         Do iZeta = 1, mZeta
            lZeta=lZeta+1
            IndZet(lZeta) = IndZ(iZeta)
            Zeta(lZeta)  = Data1(ip_Z(iZeta,nZeta))
            rKA(lZeta)   = Data1(ip_Kappa(iZeta,nZeta))
            P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))
            P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))
            P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))
            xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
            xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
            ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
            ip1 = ipOAP + mEta*mPAO*(iZeta-1)
            ip2 = ipOAP + mEta*mPAO*(lZeta-1)
            If (lDot) call dcopy_(mEta*mPAO,Scrtch(ip1),1,
     &                                     Scrtch(ip2),1)
         End Do
      End If
      If (lZeta.eq.0) Go To 999
*
      If (iphX1.ne.1) Call DScal_(lZeta,-One,P(1,1),1)
      If (iphY1.ne.1) Call DScal_(lZeta,-One,P(1,2),1)
      If (iphZ1.ne.1) Call DScal_(lZeta,-One,P(1,3),1)
*
*-----Transpose eta,mPAO,zeta to mPAO,zeta,eta
*
      If (lDot) Call DGetMO(Scrtch(ipOAP),mEta,mEta,mPAO*lZeta,
     &                      Scrtch(ipPAO),mPAO*lZeta)
*
*-----Prescreen Eta
*
      lEta=0
      If (PreScr) Then
         Do iEta = 1, mEta
            jEta = IndE(iEta)
            abcd= Data2(ip_ab(jEta,nEta)) * abMax
            If (Abs(abcd).ge.CutInt) Then
               lEta=lEta+1
               IndEta(lEta) = IndE(iEta)
               Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
               rKC(lEta)   = Data2(ip_Kappa(iEta,nEta))
               Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))
               Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))
               Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
               xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
               xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
               EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
               ip1 = ipPAO + mPAO*lZeta*(iEta-1)
               ip2 = ipPAO + mPAO*lZeta*(lEta-1)
               If (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,
     &                                         Scrtch(ip2),1)
            End If
         End Do
      Else
         Do iEta = 1, mEta
            lEta=lEta+1
            IndEta(lEta) = IndE(iEta)
            Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
            rKC(lEta)   = Data2(ip_Kappa(iEta,nEta))
            Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))
            Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))
            Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))
            xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
            xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
            EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
            ip1 = ipPAO + mPAO*lZeta*(iEta-1)
            ip2 = ipPAO + mPAO*lZeta*(lEta-1)
            if (ldot) call dcopy_(lZeta*mPAO,Scrtch(ip1),1,
     &                                      Scrtch(ip2),1)
         End Do
      End If
      If (lEta.eq.0) Go To 999
*
      If (iphX2.ne.1) Call DScal_(lEta,-One,Q(1,1),1)
      If (iphY2.ne.1) Call DScal_(lEta,-One,Q(1,2),1)
      If (iphZ2.ne.1) Call DScal_(lEta,-One,Q(1,3),1)
*
*-----Transpose mPAO,zeta,eta to zeta,eta,mPAO
*
      If (ldot) Call DGeTMO(Scrtch(ipPAO),mPAO,mPAO,lZeta*lEta,
     &                                          PAO,lZeta*lEta)
*
 999  Continue
*
      ij = 0
      Do iEta = 1,lEta
         Et   = Eta(iEta)
         rKCD = rkC(iEta)
         Do iZeta = 1,lZeta
            Zt   = Zeta(iZeta)
            rKAB = rkA(iZeta)
            ij = ij + 1
            xpre(ij) = rKAB*rKCD*Sqrt(1.0d0/(Zt+Et))
         End Do
      End Do
      If (ldot) Then
         jPAO = 0
         Do iPAO = 1, mPAO
            Do iZE = 0, lZeta*lEta-1
               jPAO=jPAO+1
               PAO(jPAO)=xpre(iZE+1)*PAO(jPAO)
            End Do
         End Do
      End If
#ifdef _DEBUG_
      If (iPrint.ge.39) Call RecPrt(' PAO',' ',
     &                              PAO,lZeta*lEta,mPAO)
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(ztmx)
         Call Unused_real(zexpmax)
         Call Unused_integer(nAlpha)
         Call Unused_integer(nBeta)
         Call Unused_real(etmx)
         Call Unused_real(eexpmax)
         Call Unused_integer(nGamma)
         Call Unused_integer(nDelta)
      End If
      End
#endif
