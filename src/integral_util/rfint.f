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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RFInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nComp,la,lb,A,B,nHer,
     &                  Array,nArr,Ccoor,nOrdOp)
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
*              GetMem                                                  *
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
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), B(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Logical ABeq(3)
*
      iRout = 122
      iPrint = nPrint(iRout)
*     iPrint = 99
      iQ = 0
      Call qEnter('RFInt')
      ABeq(1) = A(1).eq.B(1)
      ABeq(2) = A(2).eq.B(2)
      ABeq(3) = A(3).eq.B(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipRnxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
      ipTemp1 = nip
      nip = nip + nZeta
      ipTemp2 = nip
      nip = nip + nZeta
      ipTemp3 = nip
      nip = nip + 3*nZeta*nHer
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'RFInt: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in RFInt'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In RFInt: A',' ',A,1,3)
         Call RecPrt(' In RFInt: B',' ',B,1,3)
         Call RecPrt(' In RFInt: CCoor',' ',CCoor,1,3)
         Call RecPrt(' In RFInt: P',' ',P,nZeta,3)
         Write (6,*) ' In RFInt: la,lb=',la,lb
         Write (6,*) ' In RFInt: nHer=',nHer
      End If
*
*     Compute the cartesian values of the basis functions angular part
*
      Do 10 iZeta = 1, nZeta
         Array(ipTemp1-1+iZeta) = 1/Sqrt(Zeta(iZeta))
*        Array(ipTemp1-1+iZeta) = Zeta(iZeta)**(-Half)
 10   Continue
      Call vCrtCmp(Array(ipTemp1),P,nZeta,A,Array(ipAxyz),
     &               la,HerR(iHerR(nHer)),nHer,ABeq)
      Call vCrtCmp(Array(ipTemp1),P,nZeta,B,Array(ipBxyz),
     &               lb,HerR(iHerR(nHer)),nHer,ABeq)
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
     &              Array(ipAxyz),la,
     &              Array(ipRxyz),nOrdOp,
     &              Array(ipBxyz),lb,
     &              nZeta,HerW(iHerW(nHer)),nHer,Array(ipTemp3))
*
*     Combine the cartesian components to the full one electron
*     integral.
*
      Call CmbnRF(Array(ipRnxyz),nZeta,la,lb,nOrdOp,Zeta,rKappa,Final,
     &          nComp,Array(ipTemp1),Array(ipTemp2))
*
*     Call GetMem(' Exit RFInt','LIST','REAL',iDum,iDum)
      Call qExit('RFInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
      End If
      End
