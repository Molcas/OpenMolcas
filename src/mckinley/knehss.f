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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Knehss(
#define _CALLING_
#include "hss_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: to compute the gradients of the overlap matrix               *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnT2                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*     Author: Anders Bernhardsson, 1995                                *
************************************************************************
      use Her_RW
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
c#include "print.fh"

#include "hss_interface.fh"

*     Local variables
      Logical ABeq(3)
*
c     iRout = 122
c     iPrint = nPrint(iRout)
c     Call qEnter('KneHss')
*     Write (*,*) ' IfHss=',IfHss
*     Write (*,*) ' IndHss=',IndHss
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+3)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+3)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipTxyz = nip
*      nip = nip + nZeta*3*(la+3)*(lb+3)*(nOrdOp+1)
      ipRnxyz = nip
      nip = nip + nZeta*3*(la+3)*(lb+3)*(nOrdOp+1)
      ipAlph = nip
      nip = nip + nZeta
      ipBeta = nip
      nip = nip + nZeta
      If (nip-1.gt.nArr) Then
         Write (6,*) 'KneHss: nip-1.gt.nArr'
         Write (6,*) 'nip,nArr=',nip,nArr
         Call QTrace
         Call Abend()
      End If
*
c     If (iPrint.ge.49) Then
c        Call RecPrt(' In KneHss: A',' ',A,1,3)
c        Call RecPrt(' In KneHss: RB',' ',RB,1,3)
c        Call RecPrt(' In KneHss: Ccoor',' ',Ccoor,1,3)
c        Call RecPrt(' In KneHss: P',' ',P,nZeta,3)
c     End If
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la+2,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb+2,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the contribution from the multipole moment operator
*
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),
     &            nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call Assmbl(Array(ipRnxyz),
     &             Array(ipAxyz),la+2,
     &             Array(ipRxyz),nOrdOp,
     &             Array(ipBxyz),lb+2,
     &             nZeta,HerW(iHerW(nHer)),nHer)
*
*      Call Kntc(Array(ipTxyz),Array(ipRnxyz),la+2,lb,
*    &          Array(ipA),Array(ipB),nZeta)

*
*     Combine the cartesian components to the gradient of the one
*     electron integral and contract with the Fock matrix.
*
      ip = ipAlph
      Do 20 iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ip),1)
         ip = ip + nAlpha
 20   Continue
      ip = ipBeta
      Do 21 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ip),nAlpha)
         ip = ip + 1
 21   Continue
*
      Call CmbnT2(Array(ipRnxyz),nZeta,la,lb,Zeta,
     &            rKappa,Final,
     &            Array(ipAlph),Array(ipBeta),Hess,nHess,DAO,
     &            IfHss,IndHss,indgrd,dc(mdc)%nStab,dc(ndc)%nStab,nOp)

*
c     Call qExit('KneHss')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_logical_array(ifgrd)
         Call Unused_integer_array(iStabM)
         Call Unused_integer(nStabM)
      End If
      End
