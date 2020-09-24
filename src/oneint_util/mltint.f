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
      SubRoutine MltInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
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
*
#include "rmat_option.fh"
*
#include "WrkSpc.fh"
#include "oneswi.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variable
      Real*8 TC(3), Origin(3)
      Integer iStabO(0:7), iDCRT(0:7)
      Logical ABeq(3), EQ
      Character*80 Label, ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
      Data Origin/0.0D0,0.0D0,0.0D0/
*
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 122
      iPrint = nPrint(iRout)
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*     Call GetMem(' Enter MltInt','LIST','REAL',iDum,iDum)
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
*     switch (only single center overlap matrix...)
      If (NDDO.AND.
     &    .NOT.(ABeq(1).AND.ABeq(2).AND.ABeq(3))) Then
        call dcopy_(nZeta*nIC*nElem(la)*nElem(lb),[Zero],0,Final,1)
        Return
      End If
*     switch
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
      ipFnl = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
*                                                                      *
************************************************************************
*                                                                      *
      If (RMat_type_integrals) Then
         ipRnr  =  nip
         nip = nip + nZeta*(la+lb+nOrdOp+1)
      Else
         ipRnr=-1
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'MltInt: nip-1.gt.nArr*nZeta')
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
*                                                                      *
************************************************************************
*                                                                      *
      If (RMat_type_integrals) Then
*
         If (.Not.EQ(CCoor,Origin)) Then
            Call WarningMessage(2,'MltInt: R-matrix error')
            Write (6,*) 'MltInt: Wrong center of origin in case of',
     &                   ' R-matrix type of integrals!'
            Write (6,*) ' Origin should always be (0.0,0.0,0.0)!'
            Write (6,*) ' User the CENTER option to do this',
     &                  ' (see the SEWARD input sectio in the manual).'
            Write (6,'(A,I3)') 'nOrdOp=',nOrdOp
            Call Abend()
         End If
*
*        R-matrix calculations: continuum basis functions (A=B=P=0)
*        Compute the contributions of the basis functions and multipole
*        radial part
*
         lsum=la+lb+nOrdOp
         Call radlc(Zeta,nZeta,lsum,Array(ipRnr))
*
*       Combine the radial and angular component to the full one electron
*       integral.
*
         Call CmbnMPr(Array(ipRnr),nZeta,la,lb,nOrdOp,Zeta,
     &                Array(ipFnl),nComp)
*
         Call SOS(iStabO,nStabO,llOper)
         Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
         If (iPrint.ge.99) Then
            Write (6,*) ' m      =',nStabM
            Write (6,'(9A)') '{M}=',
     &                       (ChOper(iStabM(ii)),ii = 0, nStabM-1)
            Write (6,*) ' s      =',nStabO
            Write (6,'(9A)') '{S}=',
     &                       (ChOper(iStabO(ii)),ii = 0, nStabO-1)
            Write (6,*) ' LambdaT=',LmbdT
            Write (6,*) ' t      =',nDCRT
            Write (6,'(9A)') '{T}=',(ChOper(iDCRT(ii)),ii = 0, nDCRT-1)
         End If
*
         Do lDCRT = 0, nDCRT-1
*
*           Accumulate contributions
*
            nOp = NrOpr(iDCRT(lDCRT))
            Call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,Final,nIC,
     &                  nOp         ,lOper,iChO,One)
         End Do
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb,HerR(iHerR(nHer)),nHer,ABeq)
*
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
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
      Do lDCRT = 0, nDCRT-1
         Call OA(iDCRT(lDCRT),CCoor,TC)
*
*        Compute the contribution from the multipole moment operator
*
         ABeq(1) = .False.
         ABeq(2) = .False.
         ABeq(3) = .False.
         Call CrtCmp(Zeta,P,nZeta,TC,Array(ipRxyz),
     &               nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
*
*        Compute the cartesian components for the multipole moment
*        integrals. The integrals are factorized into components.
*
          Call Assmbl(Array(ipQxyz),
     &                Array(ipAxyz),la,
     &                Array(ipRxyz),nOrdOp,
     &                Array(ipBxyz),lb,
     &                nZeta,HerW(iHerW(nHer)),nHer)
*
*        Combine the cartesian components to the full one electron
*        integral.
*
         Call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,
     &               Array(ipFnl),nComp)
*
*        Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
      End If
*
*
      If (iPrint.ge.99) Then
         Write (6,*)
         Write (6,*) ' Result in MltInt'
         Write (6,*)
         Write (6,*)  'la,lb,nHer=',la,lb,nHer
         Write (6,*)  'nComp=',nComp
         Write (6,*)
         Do iIC = 1, nIC
            Write (Label,'(A,I2,A)')
     &         ' MltInt(iIC=',iIC,')'
            Call RecPrt(Label,'(10G15.8) ',Final(1,1,1,iIC),nZeta,
     &                  nElem(la)*nElem(lb))
         End Do
      End If
*
*     Call GetMem(' Exit MltInt','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
