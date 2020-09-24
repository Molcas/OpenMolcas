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
      SubRoutine MltGrd(
#define _CALLING_
#include "grd_interface.fh"
     &                 )
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
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
<<<<<<< HEAD
#include "itmax.fh"
#include "info.fh"
#ifdef _DEBUGPRINT_
=======
#ifdef _DEBUG_
>>>>>>> ce6a32feb787e2ffaaec8e43f7b8110f0a0cbe02
#include "print.fh"
#endif

#include "grd_interface.fh"

*     Local variables
      Logical ABeq(3)
      parameter (lforce=20)
      real*8 Force(lforce)
      common /finfld/Force
*                                                                      *
************************************************************************
*
#ifdef _DEBUGPRINT_
      iRout = 122
      iPrint = nPrint(iRout)
#endif
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
      If (nip-1.gt.nArr*nZeta) Then
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Call ErrTra
         Write (6,*) ' Abend in MltGrd'
         Call Abend
      End If
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Call RecPrt(' In MltGrd: RKappa',' ',rKappa,1,nZeta)
         Call RecPrt(' In MltGrd: A',' ',A,1,3)
         Call RecPrt(' In MltGrd: RB',' ',RB,1,3)
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
      Call CmbnMlt1(Array(ipRnxyz),nZeta,la,lb,Zeta,rKappa,Final,
     &            Array(ipAlph),Array(ipBeta),Grad,nGrad,DAO,
     &            IfGrad,IndGrd,dc(mdc)%nStab,dc(ndc)%nStab,
     &            kOp,nOrdOp,Force)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iStabM)
      End If
      End
