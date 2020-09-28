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
      SubRoutine Screen_g(PAO,Scrtch,mPAO,
     &                    nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                    Zeta,ZInv,P,xA,xB,Data1,nAlpha,nBeta,IndZ,
     &                    Eta, EInv,Q,xG,xD,Data2,nGamma,nDelta,IndE,
     &                    iphX1,iphY1,iphZ1,iphX2,iphY2,iphZ2,CutGrd,
     &                    l2DI,ab,abg,nab,cd,cdg,ncd,PreScr,nScrtch,
     &                    IsChi,ChiI2)
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
*             April '92 modified for gradient estimate                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "print.fh"
#include "real.fh"
      Real*8 PAO(mZeta*mEta*mPAO),Scrtch(nScrtch),
     &       Zeta(nZeta), ZInv(nZeta), P(nZeta,3),
     &        Eta(nEta),  EInv(nEta),  Q(nEta, 3),
     &       xA(nZeta), xB(nZeta), xG(nEta), xD(nEta),
     &       Data1(nZeta*(nDArray-1)), Data2(nEta*(nDArray-1)),
     &       ab(nZeta,nab), abg(nZeta,nab),
     &       cd(nEta ,ncd), cdg(nEta ,ncd),
     &       ChiI2
      Logical l2DI, PreScr, ZPreScr, EPreScr
      Integer IndZ(mZeta), IndE(mEta), IsChi
*
      iRout = 180
      iPrint = nPrint(iRout)
      LOut=6
      If (iPrint.ge.99) Then
         Call RecPrt(' Data1',' ',Data1,nZeta,nDArray-1)
         Call RecPrt(' Data2',' ',Data2,nEta ,nDArray-1)
         If (l2DI) Then
            Call RecPrt(' ab   ',' ',ab,nZeta,nab)
            Call RecPrt(' cd   ',' ',cd ,nEta ,ncd)
            Call RecPrt(' abg  ',' ',abg,nZeta,nab)
            Call RecPrt(' cdg  ',' ',cdg,nEta ,ncd)
         End If
         Call RecPrt('2nd order density matrix',' ',
     &               PAO,mZeta*mEta,mPAO)
      End If
      If (PreScr.and..Not.l2DI) Then
         Write (LOut,*) ' Screen: .Not.l2DI no activated',
     &      '  prescr=',prescr,'  l2di=',l2di
         Call Abend()
      End If
*
      Cut2=CutGrd
      lZeta=0
      lEta =0
      tOne = One
      ip=1
*
*-----Compute the prefactor. There are
*     four sections of code here which will take care of the
*     different cases of triangularization.
*
      ipFac = ip
      ip = ip + mZeta*mEta
      ij = ipFac-1
      Do iEta = 1, mEta
         Et   = Data2(ip_Z    (iEta,nEta))
         rKCD = Data2(ip_Kappa(iEta,nEta))
         Do iZeta = 1, mZeta
            Zt   = Data1(ip_Z    (iZeta,nZeta))
            rKAB = Data1(ip_Kappa(iZeta,nZeta))
            ij = ij + 1
            if (IsChi.eq.1) then
               Scrtch(ij)=rKAB*rKCD*
     &                    Sqrt(tOne/(Zt+Et+(Zt*Et*ChiI2)*Dble(IsChi)))
            else
               Scrtch(ij) = rKAB*rKCD*Sqrt(tOne/(Zt+Et))
            endif
         End Do
      End Do
      If (iPrint.ge.99) Call RecPrt(' Collected prefactors',' ',
     &                               Scrtch(ipFac),mZeta,mEta)
*
*-----Modify the 2nd order density matrix with the prefactor.
*
      ipPAO = ip
      ip = ip + mZeta*mEta*mPAO
      ipOAP = ip
      ip = ip + mZeta*mEta*mPAO
      jPAO = 1
      Do iPAO = 1, mPAO
         Do iZE = 0, mZeta*mEta-1
            PAO(jPAO+iZE)=Scrtch(ipFac+iZE)*PAO(jPAO+iZE)
         End Do
         jPAO = jPAO + mZeta*mEta
      End Do
      If (iPrint.ge.49)
     &   Call RecPrt(' Modified 2nd order density matrix',' ',
     &               PAO,mZeta*mEta,mPAO)
*
*-----Scan the modified 2nd order density matrix for the
*     largest absolute value.
*
      vMax = DNrm2_(mZeta*mEta*mPAO,PAO,1)
*-----Skip prescreening if too low.
      If (PreScr.and.Abs(vMax).lt.0.5D-4*CutGrd) Then
         lZeta = 0
         lEta  = 0
         Go To 999
      End If
*
*-----Assemble gradient estimate
*-----Multiply modified density matrix with the gradient
*     estimate.
*-----Compress all indices except zeta in ipZ
*-----Compress all indices except eta in ipE
*
      If (PreScr) Then
*
         ipZ = ip
         ip = ip + mZeta
         ipE = ip
         ip = ip + mEta
         Do i=1,mZeta+mEta
           Scrtch(ipZ+i-1)=Zero
         End Do
         Do icd = 1, ncd
            Do iEta = 1, mEta
               alpha = cd (iEta,icd)
               beta  = cdg(iEta,icd)
               iZE = (iEta-1)*mZeta
               Do iab = 1, nab
                  iabcd = (icd-1)*nab + iab
                  iOff = (iabcd-1)*mZeta*mEta+iZE
                  Do iZeta = 1, mZeta
                     temp=(alpha*abg(iZeta,iab)+beta*ab(iZeta,iab))*
     &                    PAO(iZeta+iOff)
                        Scrtch(ipZ+iZeta-1)=Scrtch(ipZ+iZeta-1)
     &                                     +Abs(temp)
                        Scrtch(ipE+iEta -1)=Scrtch(ipE+iEta -1)
     &                                     +Abs(temp)
                  End Do
               End Do
            End Do
         End Do
      Else
         ipZ = 1
         ipE = 1
      End If
*
      If (ip-1.gt.nScrtch) Then
         Write (LOut,*) ' Screen: ip-1.gt.nScrtch'
         Write (LOut,*) 'ip-1=',ip-1
         Write (LOut,*) 'nScrtch=',nScrtch
         Call Abend()
      End If
*
      rEta=DBLE(nEta*mPAO)
      qEta=DBLE(mEta*mPAO)
      rqEta=rEta/qEta
      ZPreScr=PreScr
      If (PreScr) Then
         iMin=iDMin(mZeta,Scrtch(ipZ),1)
         zMin=Scrtch(ipZ-1+iMin)
         If (zMin.ge.Cut2/rqEta) Then
            ZPreScr=.False.
         End If
         If (iPrint.ge.99) Call RecPrt(' Screening array(Eta)',' ',
     &                     Scrtch(ipZ),mZeta,1)
      End If
*
      rZeta=DBLE(nZeta*mPAO)
      qZeta=DBLE(mZeta*mPAO)
      rqZeta=rZeta/qZeta
      EPreScr=PreScr
      If (PreScr) Then
         iMin=iDMin(mEta,Scrtch(ipE),1)
         eMin=Abs(Scrtch(ipE-1+iMin))
         If (eMin.ge.Cut2/rqZeta) Then
            EPreScr=.False.
         End If
         If (iPrint.ge.99) Call RecPrt(' Screening array(Zeta)',' ',
     &                     Scrtch(ipE),mEta,1)
      End If
*
*-----Prescreen Zeta
*
      Px = One
      Py = One
      Pz = One
      If (iphX1.ne.1) Px = -Px
      If (iphY1.ne.1) Py = -Py
      If (iphZ1.ne.1) Pz = -Pz
      lZeta=0
      If (ZPreScr) Then
         Do iZeta = 1, mZeta
            If (Scrtch(ipZ+iZeta-1).ge.Cut2/rqEta) Then
               lZeta=lZeta+1
               Zeta(lZeta)  = Data1(ip_Z    (iZeta,nZeta))
               P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))*Px
               P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))*Py
               P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))*Pz
               xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
               xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
               ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
               ip2 = ipOAP + mEta*mPAO*(lZeta-1) - 1
               Do iEP = 1, mEta*mPAO
                  Scrtch(ip2+iEP)=PAO((iEP-1)*mZeta+iZeta)
               End Do
            End If
         End Do
      Else
         Do iZeta = 1, mZeta
            lZeta=lZeta+1
            Zeta(lZeta)  = Data1(ip_Z    (iZeta,nZeta))
            P(lZeta,1)   = Data1(ip_PCoor(iZeta        ,nZeta))*Px
            P(lZeta,2)   = Data1(ip_PCoor(iZeta+  nZeta,nZeta))*Py
            P(lZeta,3)   = Data1(ip_PCoor(iZeta+2*nZeta,nZeta))*Pz
            xA(lZeta)    = Data1(ip_Alpha(iZeta,nZeta,1))
            xB(lZeta)    = Data1(ip_Beta (iZeta,nZeta,2))
            ZInv(lZeta)  = Data1(ip_ZInv (iZeta,nZeta))
         End Do
         If (EPreScr) Call DGeTMO(PAO,mZeta,mZeta,mEta*mPAO,
     &                     Scrtch(ipOAP),mEta*mPAO)
      End If
      If (lZeta.eq.0) Go To 999
*
*
*-----Prescreen Eta
*
      Qx = One
      Qy = One
      Qz = One
      If (iphX2.ne.1) Qx = -Qx
      If (iphY2.ne.1) Qy = -Qy
      If (iphZ2.ne.1) Qz = -Qz
      lEta=0
      If (EPreScr) Then
         Do iEta = 1, mEta
            If (Scrtch(ipE+iEta-1).ge.Cut2/rqZeta) Then
               lEta=lEta+1
               Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
               Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))*Qx
               Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))*Qy
               Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))*Qz
               xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
               xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
               EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
               ip1 = ipOAP + iEta-1
               ip2 = ipPAO + (lEta-1)*mPAO*lZeta-1
               Do jPZ =1,mPAO*lZeta
                  Scrtch(ip2+jPZ)=Scrtch((jPZ-1)*mEta+ip1)
               End Do
            End If
         End Do
         ipP = ipPAO
         l1=mPAO
         l2=lZeta*lEta
      Else
         Do iEta = 1, mEta
            lEta=lEta+1
            Eta(lEta)   = Data2(ip_Z   (iEta,nEta))
            Q(lEta,1)   = Data2(ip_PCoor(iEta       ,nEta))*Qx
            Q(lEta,2)   = Data2(ip_PCoor(iEta+  nEta,nEta))*Qy
            Q(lEta,3)   = Data2(ip_PCoor(iEta+2*nEta,nEta))*Qz
            xG(lEta)    = Data2(ip_Alpha(iEta,nEta,1))
            xD(lEta)    = Data2(ip_Beta (iEta,nEta,2))
            EInv(lEta)  = Data2(ip_ZInv (iEta,nEta))
         End Do
         ipP = ipOAP
         l1=lEta*mPAO
         l2=lZeta
      End If
      If (lEta.eq.0) Go To 999
*
*---- Pick up the screened two-particle density.
*     .Not.PreScr            : density is already in PAO
*
*     Transpose mPAO,zeta,eta to zeta,eta,mPAO   or
*               eta,mPAO,zeta to zeta,eta,mPAO
*
      If (PreScr) Then
         If (ZPreScr.or.EPreScr) Then
            Call DGeTMO(Scrtch(ipP),l1,l1,l2,PAO,l2)
         End If
      End If
*
 999  Continue
      If (iPrint.ge.39) Call RecPrt(' PAO',' ',
     &                              PAO,lZeta*lEta,mPAO)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(IndZ)
         Call Unused_integer_array(IndE)
         Call Unused_integer(nAlpha)
         Call Unused_integer(nBeta)
         Call Unused_integer(nGamma)
         Call Unused_integer(nDelta)
      End If
      End
