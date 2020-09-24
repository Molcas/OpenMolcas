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
* Copyright (C) 1991,1995,2002, Roland Lindh                           *
************************************************************************
      SubRoutine NAInt_GIAO(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                      Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                      Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                      iStabM,nStabM,
     &                      PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of electric field         *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy  (ESSL)                                           *
*              SOS                                                     *
*              DCR                                                     *
*              XRys                                                    *
*              Util1                                                   *
*              DaXpY  (ESSL)                                           *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, January '91                             *
*                                                                      *
* Modified for explicit code, R. Lindh, February '95.                  *
*                                                                      *
* Modified for GIAOs, R. Lindh, June 2002, Tokyo, Japan.               *
************************************************************************
      Use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake,  XCff2D, XRys2D
      External TERI, MODU2, vCff2D, vRys2D
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7),
     &        lOper(nComp), iChO(nComp)
*---- Local arrays
      Real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga_EF(4), iAnga_NA(4)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 200
      iPrint = nPrint(iRout)
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
      call dcopy_(3, A,1,Coori(1,1),1)
      call dcopy_(3,RB,1,Coori(1,2),1)
*
      iAnga_EF(1) = la
      iAnga_EF(2) = lb
      iAnga_NA(1) = la
      iAnga_NA(2) = lb
      mabMin=nabSz(Max(la,lb)-1)+1
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
      mabMax=nabSz(la+lb)
      lab=(mabMax-mabMin+1)
      kab=nElem(la)*nElem(lb)
*
      iAnga_EF(3) = nOrdOp
      iAnga_EF(4) = 0
      mcdMin_EF=nabSz(nOrdOp-1)+1
      mcdMax_EF=nabSz(nOrdop)
      lcd_EF=(mcdMax_EF-mcdMin_EF+1)
      labcd_EF=lab*lcd_EF
*
      iAnga_NA(3) = nOrdOp-1
      iAnga_NA(4) = 0
      mcdMin_NA=nabSz(nOrdOp-2)+1
      mcdMax_NA=nabSz(nOrdop-1)
      lcd_NA=(mcdMax_NA-mcdMin_NA+1)
      labcd_NA=lab*lcd_NA
*
*---- Compute Flop's and size of work array which HRR will Use.
*
      Call mHRR(la,lb,nFLOP,nMem)
      nHRR=Max(labcd_EF,labcd_NA,lcd_EF*nMem,lcd_NA*nMem)
*
*---- Distribute the work array
*
      mArr = nArr - labcd_EF - nHRR
      ipEFInt = 1
      ipRys   = ipEFInt + nZeta*labcd_EF
      ipHRR   = ipRys   + nZeta*mArr
*
*---- Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
         call dcopy_(3, A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
*
*     Modify Zeta if the two-electron code will be used!
*
      If (Nuclear_Model.eq.Gaussian_Type) Then
         Do iZeta = 1, nZeta
            rKappa(iZeta)=rKappa(iZeta)*(TwoP54/Zeta(iZeta))
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over nuclear centers
*
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
            If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
*-----------Find the DCR for M and S
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab, dc(kdc+kCnt)%nStab,
     &               iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do 102 lDCRT = 0, nDCRT-1
               Call OA(iDCRT(lDCRT),C,TC)
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1, Coori(1,3),1)
               call dcopy_(3,TC,1, Coori(1,4),1)
*                                                                      *
************************************************************************
*                                                                      *
*------------- Compute integrals with the Rys-Gauss quadrature.        *
*                                                                      *
************************************************************************
*1)                                                                    *
*              Do the EF integrals
*
               nT=nZeta
               If (Nuclear_Model.eq.Gaussian_Type) Then
                  NoSpecial=.False.
                  Eta=dbsc(kCnttp)%ExpNuc
                  EInv=One/Eta
                  rKappcd=TwoP54/Eta
*                 Tag on the normalization
                  rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
                  Call Rys(iAnga_EF,nT,Zeta,ZInv,nZeta,
     &                     [Eta],[EInv],1,P,nZeta,TC,1,
     &                     rKappa,[rKappcd],Coori,Coori,CoorAC,
     &                     mabMin,mabMax,mcdMin_EF,mcdMax_EF,
     &                     Array(ipRys),mArr*nZeta,
     &                     TERI,MODU2,vCff2D,vRys2D,NoSpecial)
               Else If (Nuclear_Model.eq.Point_Charge) Then
                  NoSpecial=.True.
                  Call Rys(iAnga_EF,nT,Zeta,ZInv,nZeta,
     &                     [One],[One],1,P,nZeta,TC,1,
     &                     rKappa,[One],Coori,Coori,CoorAC,
     &                     mabMin,mabMax,mcdMin_EF,mcdMax_EF,
     &                     Array(ipRys),mArr*nZeta,
     &                     TNAI,Fake,XCff2D,XRys2D,NoSpecial)
               Else
*...more to come...
               End If
*
*------------- The integrals are now ordered as ijkl,e,f
*
*              a) Change the order to f,ijkl,e
*              b) Unfold e to ab, f,ijkl,ab
*              c) Change the order back to ijkl,ab,f
*
*a)-----------
*
               Call DGetMO(Array(ipRys),nZeta*lab,nZeta*lab,lcd_EF,
     &                     Array(ipHRR),lcd_EF)
*
*b)----------- Use the HRR to unfold e to ab
*
               Call HRR(la,lb,A,RB,Array(ipHRR),lcd_EF*nZeta,nMem,ipIn)
               ip3=ipHRR-1+ipIn
*
*c)-----------
*
               Call DGetMO(Array(ip3),lcd_EF,lcd_EF,nZeta*kab,
     &                     Array(ipEFInt),nZeta*kab)
*
*              Stored as nZeta,iElem,jElem,iComp
*                                                                      *
************************************************************************
*2)                                                                    *
*              Do the NA integrals
*
               If (Nuclear_Model.eq.Gaussian_Type) Then
                  NoSpecial=.False.
                  Eta=dbsc(kCnttp)%ExpNuc
                  EInv=One/Eta
                  rKappcd=TwoP54/Eta
*                 Tag on the normalization
                  rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
                  Call Rys(iAnga_NA,nT,Zeta,ZInv,nZeta,
     &                     [Eta],[EInv],1,P,nZeta,TC,1,
     &                     rKappa,[rKappcd],Coori,Coori,CoorAC,
     &                     mabMin,mabMax,0,0,
     &                     Array(ipRys),mArr*nZeta,
     &                     TERI,MODU2,vCff2D,vRys2D,NoSpecial)
               Else If (Nuclear_Model.eq.Point_Charge) Then
                  NoSpecial=.True.
                  Call Rys(iAnga_NA,nT,Zeta,ZInv,nZeta,
     &                     [One],[One],1,P,nZeta,TC,1,
     &                     rKappa,[One],Coori,Coori,CoorAC,
     &                     mabMin,mabMax,0,0,
     &                     Array(ipRys),mArr*nZeta,
     &                     TNAI,Fake,XCff2D,XRys2D,NoSpecial)
               Else
*...more to come...
               End If
*
*--------------Use the HRR to compute the required primitive integrals
*
               Call HRR(la,lb,A,RB,Array(ipRys),nZeta,nMem,ipNAInt)
*                                                                      *
************************************************************************
*                                                                      *
*              Assemble dV/dB
*
               Call Assemble_dVdB(Array(ipNAInt),
     &                            Array(ipEFInt),
     &                            nZeta,la,lb,A,RB,TC)
*
*                                                                      *
************************************************************************
*                                                                      *
*------- Accumulate contributions
*
               nOp = NrOpr(iDCRT(lDCRT))
               Call SymAdO(Array(ipEFInt),nZeta,la,lb,nComp,Final,nIC,
     &                     nOp,lOper,iChO,-Fact*dbsc(kCnttp)%Charge)
*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*                                                                      *
************************************************************************
*                                                                      *
      If (Nuclear_Model.eq.Gaussian_Type) Then
         Do iZeta = 1, nZeta
            rKappa(iZeta)=rKappa(iZeta)*(TwoP54/Zeta(iZeta))
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Call GetMem(' Exit EFInt','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_real_array(Ccoor)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
