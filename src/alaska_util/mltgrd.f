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
************************************************************************
      SubRoutine MltGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,B,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                  IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: to compute the gradients of the Multipole operator           *
*                                                                      *
* Called from: DrvH1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnMlt1                                                *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to gradients of the overlap matrix. October     *
*             '91.                                                     *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#ifdef _DEBUG_
#include "print.fh"
#endif
      Integer IndGrd(3,2), kOp(2), iStabM(0:nStabM-1), lOper(nComp)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), B(3),
     &       Array(nZeta*nArr), Ccoor(3), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Logical ABeq(3), IfGrad(3,2)
      parameter (lforce=20)
      real*8 Force(lforce)
      common /finfld/Force
*                                                                      *
************************************************************************
*
#ifdef _DEBUG_
      iRout = 122
      iPrint = nPrint(iRout)
#endif
      ABeq(1) = A(1).eq.B(1)
      ABeq(2) = A(2).eq.B(2)
      ABeq(3) = A(3).eq.B(3)
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
      If (nip-1.gt.nArr*nZeta) Then
         Write (6,*) ' nArr is Wrong! ', nip,' > ',nArr*nZeta
         Call ErrTra
         Write (6,*) ' Abend in MltGrd'
         Call Abend
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.49) Then
         Call RecPrt(' In MltGrd: RKappa',' ',rKappa,1,nZeta)
         Call RecPrt(' In MltGrd: A',' ',A,1,3)
         Call RecPrt(' In MltGrd: B',' ',B,1,3)
         Call RecPrt(' In MltGrd: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In MltGrd: P',' ',P,nZeta,3)
         Write (6,*) ' In MltGrd: la,lb=',la,lb
      End If
#endif
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la+1,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,B,Array(ipBxyz),
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
      Call CmbnMlt1(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,Final,
     &            Array(ipAlph),Array(ipBeta),Grad,nGrad,DAO,
     &            IfGrad,IndGrd,nStab(mdc),nStab(ndc),
     &            nIrrep,kOp,iChBas,MxFnc,nOrdOp,Force)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iStabM)
      End If
      End
