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
      SubRoutine KnEPrm(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp)
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
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Logical ABeq(3)
*
*     Call qEnter('KnEInt')
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
      nip = nip + nZeta*3*nHer*(nOrdOp-1)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+2)*(lb+2)*(nOrdOp-1)
      ipTxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)
      ipA = nip
      nip = nip + nZeta
      ipB = nip
      nip = nip + nZeta
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
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt(' In KnEInt: A',' ',A,1,3)
      Call RecPrt(' In KnEInt: RB',' ',RB,1,3)
      Call RecPrt(' In KnEInt: Ccoor',' ',Ccoor,1,3)
      Call RecPrt(' In KnEInt: P',' ',P,nZeta,3)
      Write (6,*) ' In KnEInt: la,lb=',la,lb
#endif
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
*     Compute the contribution from the multipole moment operator
*
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),
     &            nOrdOp-2,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call Assmbl(Array(ipQxyz),
     &             Array(ipAxyz),la+1,
     &             Array(ipRxyz),nOrdOp-2,
     &             Array(ipBxyz),lb+1,
     &             nZeta,HerW(iHerW(nHer)),nHer)
*
*     Compute the cartesian components for the kinetic energy integrals.
*     The kinetic energy components are linear combinations of overlap
*     components.
*
      ipAOff = ipA
      Do 200 iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
         ipAOff = ipAOff + nAlpha
 200  Continue
*
      ipBOff = ipB
      Do 210 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
 210  Continue
*
      Call Kntc(Array(ipTxyz),Array(ipQxyz),la,lb,
     &          Array(ipA),Array(ipB),nZeta)
*
*     Combine the cartesian components to the full one electron
*     integral.
*
      Call CmbnKE(Array(ipQxyz),nZeta,la,lb,nOrdOp-2,Zeta,rKappa,Final,
     &            nComp,Array(ipTxyz))
*
*
*     Call qExit('KnEInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
      End If
      End
