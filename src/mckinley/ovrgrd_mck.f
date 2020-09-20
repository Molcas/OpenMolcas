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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine OvrGrd_mck(
#define _CALLING_
#include "grd_mck_interface.fh"
     &                     )
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
*              CmbnS1_mck                                              *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
*     Author: Anders Bernhardsson                                      *
*             November '90                                             *
*                                                                      *
*             Modified to gradients of the overlap matrix. October     *
*             '91.                                                     *
*             Modified for respons calculation in May '95  By          *
*             Anders Bernhardsson                                      *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
c#include "print.fh"

#include "grd_mck_interface.fh"

*     Local variables
      Logical ABeq(3)
*
*     Statement function for Cartesian index
*
      nElem(la)=(la+2)*(la+1)/2
*
c     iRout = 122
      iprint=0
c     iPrint = nPrint(iRout)
*     Write (*,*) ' IfGrad=',IfGrad
*     Write (*,*) ' IndGrd=',IndGrd
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
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipRnxyz = nip
      nip = nip + nZeta*3*(la+2)*(lb+2)*(nOrdOp+1)
      ipAlph = nip
      nip = nip + nZeta
      ipBeta = nip
      nip = nip + nZeta
      ipScrt=nip
      nip=nip+nElem(la)*nElem(lb)*nZeta*2


      If (nip-1.gt.nArr) Then
         Write (6,*) 'OvrGrd_Mck: nip-1.gt.nArr'
         Write (6,*) 'nip,nArr=',nip,nArr
         Call QTrace
         Call Abend()
      End If
*
c     If (iPrint.ge.49) Then
c        Call RecPrt(' In OvrGrd: A',' ',A,1,3)
c        Call RecPrt(' In OvrGrd: RB',' ',RB,1,3)
c        Call RecPrt(' In OvrGrd: Ccoor',' ',Ccoor,1,3)
c        Call RecPrt(' In OvrGrd: P',' ',P,nZeta,3)
c        Write (*,*) ' In OvrGrd: la,lb=',la,lb
c     End If
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
     &            nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call Assmbl(Array(ipRnxyz),
     &             Array(ipAxyz),la+1,
     &             Array(ipRxyz),nOrdOp,
     &             Array(ipBxyz),lb+1,
     &             nZeta,HerW(iHerW(nHer)),nHer)
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
      Call CmbnS1_mck(Array(ipRnxyz),nZeta,la,lb,Zeta,
     &            rKappa,Array(ipScrt),
     &            Array(ipAlph),Array(ipBeta),IfGrad,nOp)
*
      If (iPrint.ge.49)
     &    Call RecPrt(' Primitive Integrals',' ',
     &                Array(ipScrt),nZeta,
     &                nElem(la)*nElem(lb))
*
*
*     Symmetry adopt the gradient operator
*

      Call SymAdO_mck(Array(ipScrt),nZeta*nElem(la)*nElem(lb),
     &            Final,nrOp,
     &            nop,loper,IndGrd,iu,iv,ifgrad,idcar,trans)
      If (iPrint.ge.49)
     &    Call RecPrt(' Primitive Integrals SO',' ',
     &                Final,nZeta,
     &                nElem(la)*nElem(lb)*nrOp)

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer(iDCnt)
         Call Unused_integer_array(iStabM)
         Call Unused_integer(nStabM)
      End If
      End
