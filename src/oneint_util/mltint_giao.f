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
      SubRoutine MltInt_GIAO(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,
     &                       P,Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                       Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                       iStabM,nStabM,
     &                       PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: to compute the multipole moments integrals with the          *
*         Gauss-Hermite quadrature.                                    *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              SOS                                                     *
*              DCR                                                     *
*              Assmbl                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnMP                                                  *
*              SymAdO                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), RAB(3),
     &       Array(nZeta*nArr), Ccoor(3), TC(3)
      Character*80 Label, ChOper(0:7)*3
      Integer lOper(nComp), iStabM(0:nStabM-1), iStabO(0:7),
     &          iDCRT(0:7), iChO(nComp)
      Logical ABeq(3), EQ
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 122
      iPrint = nPrint(iRout)
*     Call qEnter('MltInt_GIAO')
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
      If (EQ(A,RB)) Go To 999
*
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
      RAB(1)=A(1)-RB(1)
      RAB(2)=A(2)-RB(2)
      RAB(3)=A(3)-RB(3)
*     switch (only single center overlap matrix...)
      If (NDDO.AND.
     &    .NOT.(ABeq(1).AND.ABeq(2).AND.ABeq(3))) Then
        call dcopy_(nZeta*nIC*nElem(la)*nElem(lb),Zero,0,Final,1)
*       Call qExit('MltInt_GIAO')
        Return
      End If
*     switch
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+2)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+2)
      ipFnl = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'MltInt_GIAO: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in MltInt'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In MltInt: A',' ',A,1,3)
         Call RecPrt(' In MltInt: RB',' ',RB,1,3)
         Call RecPrt(' In MltInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In MltInt: Kappa',' ',rKappa,nAlpha,nBeta)
         Call RecPrt(' In MltInt: Zeta',' ',Zeta,nAlpha,nBeta)
         Call RecPrt(' In MltInt: P',' ',P,nZeta,3)
         Write (6,*) ' In MltInt: la,lb=',la,lb
      End If
*
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb,HerR(iHerR(nHer)),nHer,ABeq)
*
*
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
      If (iPrint.ge.99) Then
         Write (6,*) ' m      =',nStabM
         Write (6,'(9A)') '{M}=',(ChOper(iStabM(ii)),ii = 0, nStabM-1)
         Write (6,*) ' s      =',nStabO
         Write (6,'(9A)') '{S}=',(ChOper(iStabO(ii)),ii = 0, nStabO-1)
         Write (6,*) ' LambdaT=',LmbdT
         Write (6,*) ' t      =',nDCRT
         Write (6,'(9A)') '{T}=',(ChOper(iDCRT(ii)),ii = 0, nDCRT-1)
      End If
*
      Do 102 lDCRT = 0, nDCRT-1
         TC(1) = iPhase(1,iDCRT(lDCRT))*CCoor(1)
         TC(2) = iPhase(2,iDCRT(lDCRT))*CCoor(2)
         TC(3) = iPhase(3,iDCRT(lDCRT))*CCoor(3)
*
*        Compute the contribution from the multipole moment operator
*
         ABeq(1) = .False.
         ABeq(2) = .False.
         ABeq(3) = .False.
         Call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),
     &               nOrdOp+1,HerR(iHerR(nHer)),nHer,ABeq)
*
*        Compute the cartesian components for the multipole moment
*        integrals. The integrals are factorized into components.
*
          Call Assmbl(Array(ipQxyz),
     &                Array(ipAxyz),la,
     &                Array(ipRxyz),nOrdOp+1,
     &                Array(ipBxyz),lb,
     &                nZeta,HerW(iHerW(nHer)),nHer)
*
*        Combine the cartesian components to the full one electron
*        integral.
*
         nB=3
         Call CmbnMP_GIAO(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,
     &                    Array(ipFnl),nComp/nB,nB,RAB,TC)
*
*        Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 102  Continue
*
 999  Continue
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in MltInt'
         Do 100 ia = 1, (la+1)*(la+2)/2
            Do 200 ib = 1, (lb+1)*(lb+2)/2
               Do 300 iIC = 1, nIC
                  Write (Label,'(A,I2,A,I2,A,I2,A)')
     &               ' Final(a=',ia,',b=',ib,',iIC=',iIC,')'
                  Call RecPrt(Label,' ',Final(1,ia,ib,iIC),nAlpha,nBeta)
 300           Continue
 200        Continue
 100     Continue
      End If
*
*     Call qExit('MltInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
