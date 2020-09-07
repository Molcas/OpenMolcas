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
      SubRoutine KnEInt_GIAO(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,
     &                       P,Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                       Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                       iStabM,nStabM,
     &                       PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: to compute the kinetic energy integrals with the Gauss-      *
*         Hermite quadrature.                                          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              Kntc                                                    *
*              CmbnKE                                                  *
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
*
#include "rmat_option.fh"
#include "rmat.fh"
*
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Integer lOper(nComp), iStabM(0:nStabM-1), iStabO(0:7),
     &        iDCRT(0:7), iChO(nComp)
      Logical ABeq(3)
*
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
*     Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*     iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
*     Index(ixyz,ix,iz) = Ind(ixyz,ix,iz) + iOff(ixyz)
*
      iRout = 150
      iPrint = nPrint(iRout)
*     Call qEnter('KnEInt_GIAO')
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+2)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+2)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+2)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+2)*(lb+2)*(nOrdOp+2)
      ipTxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+2)
      ipWxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*2
      ipA = nip
      nip = nip + nZeta
      ipB = nip
      nip = nip + nZeta
      ipFnl = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
*                                                                      *
************************************************************************
*                                                                      *
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'KNEInt: nip-1.gt.nArr*nZeta')
         Write (6,*) 'nip=',nip
         Write (6,*) 'nArr,nZeta=',nArr,nZeta
         Call  Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In KnEInt: A',' ',A,1,3)
         Call RecPrt(' In KnEInt: RB',' ',RB,1,3)
         Call RecPrt(' In KnEInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In KnEInt: P',' ',P,nZeta,3)
         Write (6,*) ' In KnEInt: la,lb=',la,lb
      End If
*                                                                      *
************************************************************************
*                                                                      *
      llOper = lOper(1)
      Do iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la+1,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb+1,HerR(iHerR(nHer)),nHer,ABeq)
*
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
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
     &               nOrdOp+1,HerR(iHerR(nHer)),nHer,ABeq)
*
*        Compute the cartesian components for the multipole moment
*        integrals. The integrals are factorized into components.
*
         Call Assmbl(Array(ipQxyz),
     &               Array(ipAxyz),la+1,
     &               Array(ipRxyz),nOrdOp+1,
     &               Array(ipBxyz),lb+1,
     &               nZeta,HerW(iHerW(nHer)),nHer)
*
*        Compute the cartesian components for the kinetic energy
*        integrals. The kinetic energy components are linear
*        combinations of overlap components.
*
         ipAOff = ipA
         Do iBeta = 1, nBeta
            call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
            ipAOff = ipAOff + nAlpha
         End Do
*
         ipBOff = ipB
         Do iAlpha = 1, nAlpha
            call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
            ipBOff = ipBOff + 1
         End Do
*
         Call Kntc_GIAO(Array(ipTxyz),Array(ipQxyz),Array(ipWxyz),
     &                  la,lb,nOrdOp,
     &                  Array(ipA),Array(ipB),nZeta)
*
*        Combine the cartesian components to the full one electron
*        integral.
*
         nB=3
         Call CmbnKE_GIAO(Array(ipQxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,
     &                    Array(ipFnl),nComp/nB,nB,Array(ipTxyz),
     &                    Array(ipWxyz),A,RB,TC)
*
*        Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipFnl),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
*     Call qExit('KnEInt_GIAO')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
