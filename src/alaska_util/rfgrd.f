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
* Copyright (C) 1990,1992,1995, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RFGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,B,nHer,
     &                 Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                 IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                 iStabM,nStabM)
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
*              Assmbl                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnMP                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
*                                                                      *
*             Roland Lindh, Dept. of Theoratical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to reaction field calculations July '92         *
*             Modified to gradient calculations May '95                *
************************************************************************
      use Her_RW
      use PCM_arrays, only: MM
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "rctfld.fh"
      Integer IndGrd(3,2), kOp(2), iStabM(0:nStabM-1),lOper(nComp)
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), B(3),
     &       Array(nZeta*nArr), Ccoor(3), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+1)/2,(lb+1)*(lb+2)/2)
      Logical ABeq(3), IfGrad(3,2)
*
      iRout = 122
      iPrint = nPrint(iRout)
*     iPrint = 99
      Call qEnter('RFGrd')
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
      ipTemp1 = nip
      nip = nip + nZeta
      ipTemp2 = nip
      nip = nip + nZeta
      ipTemp3 = nip
      nip = nip + 3*nZeta*nHer
      ipAlph = nip
      nip = nip + nZeta
      ipBeta = nip
      nip = nip + nZeta
      If (nip-1.gt.nArr*nZeta) Then
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Call ErrTra
         Write (6,*) ' Abend in RFGrd'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In RFGrd: A',' ',A,1,3)
         Call RecPrt(' In RFGrd: B',' ',B,1,3)
         Call RecPrt(' In RFGrd: CCoor',' ',CCoor,1,3)
         Call RecPrt(' In RFGrd: P',' ',P,nZeta,3)
         Write (6,*) ' In RFGrd: la,lb=',la,lb
         Write (6,*) ' In RFGrd: nHer=',nHer
      End If
*
*     Compute the cartesian values of the basis functions angular part
*
      Do 10 iZeta = 1, nZeta
         Array(ipTemp1-1+iZeta) = Zeta(iZeta)**(-Half)
 10   Continue
*
      Call vCrtCmp(Array(ipTemp1),P,nZeta,A,Array(ipAxyz),
     &               la+1,HerR(iHerR(nHer)),nHer,ABeq)
      Call vCrtCmp(Array(ipTemp1),P,nZeta,B,Array(ipBxyz),
     &               lb+1,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the contribution from the multipole moment operator
*
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call vCrtCmp(Array(ipTemp1),P,nZeta,Ccoor,Array(ipRxyz),
     &            nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call vAssmbl(Array(ipRnxyz),
     &              Array(ipAxyz),la+1,
     &              Array(ipRxyz),nOrdOp,
     &              Array(ipBxyz),lb+1,
     &              nZeta,HerW(iHerW(nHer)),nHer,Array(ipTemp3))
*
*     Combine the cartesian components to the full one electron
*     integral.
*
      ip = ipAlph
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ip),1)
         ip = ip + nAlpha
      End Do
      ip = ipBeta
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ip),nAlpha)
         ip = ip + 1
      End Do
      nCav=(nOrdOp+1)*(nOrdOp+2)*(nOrdOp+3)/6
      Call CmbnRF1(Array(ipRnxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Final,
     &             nComp,Array(ipTemp1),Array(ipTemp2),
     &             Array(ipAlph),Array(ipBeta),Grad,nGrad,DAO,
     &             IfGrad,IndGrd,dc(mdc)%nStab,dc(ndc)%nStab,nIrrep,
     &             kOp,iChBas,MxFnc,MM(1,2))
*
      Call qExit('RFGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iStabM)
      End If
      End
