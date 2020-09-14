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
      SubRoutine KnEGrd_mck(Alpha,nAlpha,Beta, nBeta,
     &                 Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,B,nHer,
     &                 Array,nArr,Ccoor,nOrdOp,
     &                 IfGrad,IndGrd,nOp,
     &                 lOper,iu,iv,nrOp,iDCar,iDCnt,iStabM,nStabM,trans)
************************************************************************
*                                                                      *
* Object: to compute the gradient of the kinetic energy integrals      *
*         with the Gauss-Hermite quadrature                            *
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
*              CmbnT1                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Anders Bernhardsson,1995                                 *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
c#include "print.fh"
      Integer IndGrd(0:nIrrep), nOp(2)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), B(3),
     &       Array(nArr), Ccoor(3)
      Logical ABeq(3), IfGrad(3,2),trans(2)
*
*     Statement function for Cartesian index
*
      nElem(li)=(li+1)*(li+2)/2
*
c     iRout = 150
c     iPrint = nPrint(iRout)
c     Call qEnter('KnEGrd')
      ABeq(1) = A(1).eq.B(1)
      ABeq(2) = A(2).eq.B(2)
      ABeq(3) = A(3).eq.B(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+3)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+3)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipRnxyz = nip
      nip = nip + nZeta*3*(la+3)*(lb+3)*(nOrdOp+1)
      ipTxyz = nip
      nip = nip + nZeta*3*(la+2)*(lb+2)
      ipA = nip
      nip = nip + nZeta
      ipB = nip
      nip = nip + nZeta
      ipSc=nip
      nip=nip+nElem(la)*nElem(lb)*nZeta
      If (nip-1.gt.nArr) Then
         Write (6,*) 'KneGrd_Mck: nip-1.gt.nArr'
         Write (6,*) 'nip,nArr=',nip,nArr
         Call QTrace
         Call Abend()
      End If
*
c     Call GetMem(' Beg KnEGrd','CHEC','REAL',iDum,iDum)
c     If (iPrint.ge.49) Then
c        Call RecPrt(' In KnEGrd: A',' ',A,1,3)
c        Call RecPrt(' In KnEGrd: B',' ',B,1,3)
c        Call RecPrt(' In KnEGrd: Ccoor',' ',Ccoor,1,3)
c        Call RecPrt(' In KnEGrd: P',' ',P,nZeta,3)
c        Write (*,*) ' In KnEGrd: la,lb=',la,lb
c     End If
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la+2,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,B,Array(ipBxyz),
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
      Call Kntc(Array(ipTxyz),Array(ipRnxyz),la+1,lb+1,
     &          Array(ipA),Array(ipB),nZeta)
*
*     Combine the cartesian components to the gradient of the kinetic
*     energy integral and trace with the variational density matrix.
*

      Call CmbnT1_mck(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,
     &            Array(ipSc),Array(ipTxyz),
     &            Array(ipA),Array(ipB),IfGrad)
*
*?
      call dcopy_(nElem(la)*nElem(lb)*nZeta*NrOp,[Zero],0,Final,1)
*
*
*     Symmetry adopt the gradient operator
*
      Call SymAdO_mck(Array(ipSc),nZeta*nElem(la)*nElem(lb),
     &                Final,nrOp,
     &                nop,loper,IndGrd,iu,iv,ifgrad,idCar,trans)

*
c     Call qExit('KnEGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer(iDCnt)
         Call Unused_integer(iStabM)
         Call Unused_integer(nStabM)
      End If
      End
