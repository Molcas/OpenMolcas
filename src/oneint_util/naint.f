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
      SubRoutine NAInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                 Array,nArr,CCoor,nOrdOp,lOper,iChO,
     &                 iStabM,nStabM,
     &                 PtChrg,nGrid,iAddPot)
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
      Use Basis_Info
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
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,nComp),
     &       Array(nZeta*nArr)
      Integer iStabM(0:nStabM-1), lOper(nComp)
*-----Local arrys
      Real*8 C(3), TC(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial, No3Cnt
      Integer iAnga(4), iDCRT(0:7), iChO(nComp)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 151
      iPrint = nPrint(iRout)
C     Call qEnter('NAInt')
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
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
      No3Cnt = .FALSE.
      If (EQ(A,RB)) Then
        mabMin=nabSz(la+lb-1)+1
      Else If (NDDO) Then
        No3Cnt = .TRUE.
      End If
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
      llOper = lOper(1)
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
*     Loop over nuclear centers.
*
      kdc = 0
      Do 100 kCnttp = 1, nCnttp
*
*        Change nuclear charge if this is a relativistic ECP-case. This
*        is used for the DKH transformation (see dkh_util/dkrelint.f)!
*
         If (DKroll.and.Primitive_Pass.and.lECP) Then
            Q_Nuc=DBLE(dbsc(kCnttp)%AtmNr)
         Else
            Q_Nuc=dbsc(kCnttp)%Charge
         End If

         If (Q_Nuc.eq.Zero) Go To 111
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
            If (iPrint.ge.99) Call RecPrt('C',' ',C,1,3)
*
*-----------Find the DCR for M and S
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt) ,nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            If (iPrint.ge.99) Then
               Write (6,*) ' m      =',nStabM
               Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &               ii = 0, nStabM-1)
               Write (6,*) ' s      =',nStab(kdc+kCnt)
               Write (6,'(9A)') '(S)=',(ChOper(jStab(ii,kdc+kCnt)),
     &               ii = 0, nStab(kdc+kCnt)-1)
               Write (6,*) ' LambdaT=',LmbdT
               Write (6,*) ' t      =',nDCRT
               Write (6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),
     &               ii = 0, nDCRT-1)
            End If

*
            Do 102 lDCRT = 0, nDCRT-1
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
*              switch (only two center NA matrix...)
               If (No3Cnt .AND. .NOT.(EQ(A,TC).OR.EQ(RB,TC))) Go To 102
*              switch
               call dcopy_(3,TC,1,CoorAC(1,2),1)
               call dcopy_(3,TC,1,Coori(1,3),1)
               call dcopy_(3,TC,1,Coori(1,4),1)
               call dcopy_(3,TC,1,Coora(1,3),1)
               call dcopy_(3,TC,1,Coora(1,4),1)
*                                                                      *
************************************************************************
*                                                                      *
*              Compute integrals with the Rys quadrature.
*                                                                      *
************************************************************************
*                                                                      *
               nT = nZeta
*                                                                      *
************************************************************************
*                                                                      *
               If (Nuclear_Model.eq.Gaussian_Type) Then
*
*                 Gaussian nuclear charge distribution
*
                  NoSpecial=.False.
                  Eta=dbsc(kCnttp)%ExpNuc
                  EInv=One/Eta
                  rKappcd=TwoP54/Eta
*                 Tag on the normalization
                  rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
*                 s-type function
                  mcdMin=0
                  mcdMax=0
                  Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                     [Eta],[EInv],1,P,nZeta,
     &                     TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,
     &                     mabmin,mabmax,mcdMin,mcdMax,
     &                     Array,nArr*nZeta,
     &                     TERI,ModU2,vCff2D,vRys2D,NoSpecial)
*                                                                      *
************************************************************************
*                                                                      *
               Else If (Nuclear_Model.eq.mGaussian_Type) Then
*
*                 Modified Gaussian nuclear charge distribution
*
                  NoSpecial=.False.
                  Eta=dbsc(kCnttp)%ExpNuc
                  EInv=One/Eta
                  rKappcd=TwoP54/Eta
*                 Tag on the normalization
                  rKappcd=rKappcd*(Eta/Pi)**(Three/Two)
     &                   /(One+Three*w_mGauss(kCnttp)/(Two*Eta))
*                 s type function
                  mcdMin=0
                  mcdMax=0
                  Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                     [Eta],[EInv],1,P,nZeta,
     &                     TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,
     &                     mabmin,mabmax,mcdMin,mcdMax,
     &                     Array,nArr*nZeta,
     &                     TERI,ModU2,vCff2D,vRys2D,NoSpecial)
*
*                 d type function w*(x**2+y**2+z**2)
                  If (w_mGauss(kCnttp).gt.0.0D0) Then
                     rKappcd = rKappcd*w_mGauss(kCnttp)
                     iAnga(3)=2
                     mcdMin=nabSz(2+ld-1)+1
                     mcdMax = nabSz(2+ld)
*                    tweak the pointers
                     ipOff = 1 + nZeta * (la+1)*(la+2)/2
     &                                 * (lb+1)*(lb+2)/2
                     mArr = nArr - (la+1)*(la+2)/2*(lb+1)*(lb+2)/2
                     Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                        [Eta],[EInv],1,P,nZeta,
     &                        TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,
     &                        mabMin,mabMax,mcdMin,mcdMax,
     &                        Array(ipOff),mArr*nZeta,
     &                        TERI,ModU2,vCff2D,vRys2D,NoSpecial)
                     iAnga(3)=0
*
*                    Add the s and d contributions together!
*
                     Call Assemble_mGauss(Array,Array(ipOff),
     &                                    nZeta*(mabMax-mabMin+1))
                  End If
*                                                                      *
************************************************************************
*                                                                      *
               Else If (Nuclear_Model.eq.Point_Charge) Then
*
*                 Point-like nuclear charge distribution
*
                  NoSpecial=.True.
                  Eta=One
                  EInv=One
                  rKappcd=One
                  mcdMin=0
                  mcdMax=0
                  Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                     [Eta],[EInv],1,P,nZeta,
     &                     TC,1,rKappa,[rKappcd],Coori,Coora,CoorAC,
     &                     mabMin,mabMax,mcdMin,mcdMax,
     &                     Array,nArr*nZeta,
     &                     TNAI,Fake,XCff2D,XRys2D,NoSpecial)
               End If
*                                                                      *
************************************************************************
*                                                                      *
*
*--------------Use the HRR to compute the required primitive integrals.
*
               Call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)
*
*--------------Accumulate contributions to the symmetry adapted operator
*
               nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
               Call SymAdO(Array(ipIn),nZeta,la,lb,nComp,Final,nIC,
     &                     nOp         ,lOper,iChO,-Fact*Q_Nuc)
               If (iPrint.ge.99) Then
                  Write (6,*) Fact*Q_Nuc
                  Call RecPrt('NaInt: Array(ipIn)',' ',Array(ipIn),
     &                    nZeta,nElem(la)*nElem(lb)*nComp)
                  Call RecPrt('NaInt: Final',' ',Final,
     &                    nZeta,nElem(la)*nElem(lb)*nIC)
               End If
*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*
      If (Nuclear_Model.eq.Gaussian_Type .or.
     &    Nuclear_Model.eq.mGaussian_Type) Then
         Do iZeta = 1, nZeta
            rKappa(iZeta)=rKappa(iZeta)/(TwoP54/Zeta(iZeta))
         End Do
      End If
*
*     Call GetMem(' Exit NAInt','LIST','REAL',iDum,iDum)
C     Call qExit('NAInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_real_array(CCoor)
         Call Unused_integer(nOrdOp)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
