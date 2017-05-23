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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine NAPrm(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,nComp,la,lb,A,RB,nRys,
     &                 Array,nArr,CCoor,nOrdOp)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              mHrr                                                    *
*              DCR                                                     *
*              Rys                                                     *
*              Hrr                                                     *
*              DaXpY   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, January 1991                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*     Used for normal nuclear attraction integrals
      External TNAI, Fake, XCff2D, XRys2D
*     Used for finite nuclei
      External TERI, ModU2, vCff2D, vRys2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,2),
     &       Array(nZeta*nArr)
*-----Local arrys
      Real*8 C(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga(4)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 151
      iPrint = nPrint(iRout)
C     Call qEnter('NAPrm')
*
      Call FZero(Final,nZeta*nElem(la)*nElem(lb)*nComp)
*
      lc=0
      ld=0
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = lc
      iAnga(4) = ld
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(2*3,Coora,1,Coori,1)
      mabMin = nabSz(Max(la,lb)-1)+1
      mabMax = nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
*
*     Compute FLOPs and size of work array which Hrr will use.
*
      Call mHrr(la,lb,nFLOP,nMem)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
       call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
       call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
*     Modify Zeta if the two-electron code will be used!
*
      If (Nuclear_Model.eq.Gaussian_Type .or.
     &    Nuclear_Model.eq.mGaussian_Type) Then
         Do iZeta = 1, nZeta
            rKappa(iZeta)=rKappa(iZeta)*(TwoP54/Zeta(iZeta))
         End Do
      End If
*
      iCnttp = INT(CCoor(1,2))
      Q_Nuc=Charge(iCnttp)

      If (Q_Nuc.eq.Zero) Go To 111
      call dcopy_(3,CCoor,1,C,1)
      If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
      Call DCopy_(3,C,1,CoorAC(1,2),1)
      Call DCopy_(3,C,1,Coori(1,3),1)
      Call DCopy_(3,C,1,Coori(1,4),1)
      Call DCopy_(3,C,1,Coora(1,3),1)
      Call DCopy_(3,C,1,Coora(1,4),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute integrals with the Rys quadrature.
*                                                                      *
************************************************************************
*                                                                      *
      nT = nZeta
*                                                                      *
************************************************************************
*                                                                      *
      If (Nuclear_Model.eq.Gaussian_Type) Then
*
*        Gaussian nuclear charge distribution
*
         NoSpecial=.False.
         Eta=ExpNuc(iCnttp)
         EInv=One/Eta
         rKappcd=TwoP54/Eta
*        Tag on the normalization
         rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
*        s-type function
         mcdMin=0
         mcdMax=0
         Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &            Eta,EInv,1,P,nZeta,
     &            C,1,rKappa,rKappcd,Coori,Coora,CoorAC,
     &            mabmin,mabmax,mcdMin,mcdMax,
     &            Array,nArr*nZeta,
     &            TERI,ModU2,vCff2D,vRys2D,NoSpecial)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Nuclear_Model.eq.mGaussian_Type) Then
*
*        Modified Gaussian nuclear charge distribution
*
         NoSpecial=.False.
         Eta=ExpNuc(iCnttp)
         EInv=One/Eta
         rKappcd=TwoP54/Eta
*        Tag on the normalization
         rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
     &          /(One+Three*w_mGauss(iCnttp)/(Two*Eta))
*        s type function
         mcdMin=0
         mcdMax=0
         Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &            Eta,EInv,1,P,nZeta,
     &            C,1,rKappa,rKappcd,Coori,Coora,CoorAC,
     &            mabmin,mabmax,mcdMin,mcdMax,
     &            Array,nArr*nZeta,
     &            TERI,ModU2,vCff2D,vRys2D,NoSpecial)
*
*        d type function w*(x**2+y**2+z**2)
         If (w_mGauss(iCnttp).gt.0.0D0) Then
            rKappcd = rKappcd*w_mGauss(iCnttp)
            iAnga(3)=2
            mcdMin=nabSz(2+ld-1)+1
            mcdMax = nabSz(2+ld)
*           tweak the pointers
            ipOff = 1 + nZeta * (la+1)*(la+2)/2
     &            * (lb+1)*(lb+2)/2
            mArr = nArr - (la+1)*(la+2)/2*(lb+1)*(lb+2)/2
            Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &               Eta,EInv,1,P,nZeta,
     &               C,1,rKappa,rKappcd,Coori,Coora,CoorAC,
     &               mabMin,mabMax,mcdMin,mcdMax,
     &               Array(ipOff),mArr*nZeta,
     &               TERI,ModU2,vCff2D,vRys2D,NoSpecial)
            iAnga(3)=0
*
*           Add the s and d contributions together!
*
            Call Assemble_mGauss(Array,Array(ipOff),
     &                           nZeta*(mabMax-mabMin+1))
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Nuclear_Model.eq.Point_Charge) Then
*
*        Point-like nuclear charge distribution
*
         NoSpecial=.True.
         Eta=One
         EInv=One
         rKappcd=One
         mcdMin=0
         mcdMax=0
         Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &            Eta,EInv,1,P,nZeta,
     &            C,1,rKappa,rKappcd,Coori,Coora,CoorAC,
     &            mabMin,mabMax,mcdMin,mcdMax,
     &            Array,nArr*nZeta,
     &            TNAI,Fake,XCff2D,XRys2D,NoSpecial)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Use the HRR to compute the required primitive integrals.
*
      Call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)
*
      Call DCopy_(nZeta*nElem(la)*nElem(lb)*nComp,Array(ipIn),1,Final,1)
      Call DScal_(nZeta*nElem(la)*nElem(lb)*nComp,-Q_Nuc,Final,1)
*
 111  Continue
*
      If (Nuclear_Model.eq.Gaussian_Type .or.
     &    Nuclear_Model.eq.mGaussian_Type) Then
         Do iZeta = 1, nZeta
            rKappa(iZeta)=rKappa(iZeta)/(TwoP54/Zeta(iZeta))
         End Do
      End If
*
C     Call qExit('NAPrm')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_integer(nOrdOp)
      End If
      End
