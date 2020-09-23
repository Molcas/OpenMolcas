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
      SubRoutine ElGrddot
     &           (Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 nZeta,la,lb,A,B,nHer,
     &                 Array,nArr,Ccoor,nOrdOp,rout,
     &                 IndGrd,DAO,mdc,ndc,nOp,
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
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnMP                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
*                                                                      *
*             Modified to reaction field calculations July '92         *
*             Modified to gradient calculations May '95                *
************************************************************************
      use Her_RW
      use Center_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Integer IndGrd(2,3,3,0:nirrep-1), nOp(2), iStabM(0:nStabM-1)
      Real*8
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), B(3),
     &       Array(nZeta*nArr), Ccoor(3), rout(*),
     &       DAO(nZeta,(la+1)*(la+1)/2,(lb+1)*(lb+2)/2)
      Logical ABeq(3)
*
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2

*
      iprint=0
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
      ipFinal=nip
      nip=nip+nzeta*nElem(la)*nElem(lb)*4*6
      If (nip-1.gt.nArr*nZeta) Then
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Call ErrTra
         Write (6,*) ' Abend in RFGrd'
         Call Abend
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
      ncomp=4
      Call Cmbneldot(Array(ipRnxyz),nZeta,la,lb,nOrdOp,Zeta,
     &            rKappa,Array(ipFinal),
     &             ncomp,Array(ipTemp1),Array(ipTemp2),
     &             Array(ipAlph),Array(ipBeta),DAO,
     &             dc(mdc)%nStab,dc(ndc)%nStab,nOp,rout,indgrd)
*o
*     Call GetMem(' Exit RFGrd','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(iStabM)
      End If
      End
