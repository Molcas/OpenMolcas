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
      SubRoutine KnEInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
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
*
#include "rmat_option.fh"
#include "rmat.fh"
*
#include "print.fh"

#include "int_interface.fh"
*     Local variables
      Logical ABeq(3)
*
      iRout = 150
      iPrint = nPrint(iRout)
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
      If (RMat_type_integrals) Then
         ipRnr=nip
         nip  = nip + nZeta*(la+lb+3)
         ipqC = nip
         nip  = nip + nZeta*(la+lb+1)
         ipDi =nip
         nip  = nip + nZeta*(la+lb+1)
      Else
         ipRnr=-1
         ipqC =-1
         ipDi =-1
      End If
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
*
      If (RMat_type_integrals) Then
*                                                                      *
************************************************************************
*                                                                      *
*     R-matrix calculations: continuum basis functions (A=B=P=0)
*      Compute the contributions of the basis functions and multipole
*      radial part
*
       lsum=la+lb+2
       Call radlc(Zeta,nZeta,lsum,Array(ipRnr))
*
*     Optional for photoionization:
*     R-matrix calculations: continuum basis functions (A=B=P=0)
*      Compute the contributions of the Coulomb operator times qCoul
*      outside the sphere Omega
*
       if(abs(qCoul).gt.Epsq) then
        lsum=la+lb
        icop=1
        Call radlq(Zeta,nzeta,lsum,Array(ipqC),icop)
       endif
*
       if(abs(dipol1).gt.Epsq) then
        lsum=la+lb
        icop=2
        Call radlq(Zeta,nzeta,lsum,Array(ipDi),icop)
       endif
*
*     Combine the radial and angular component to the full one electron
*     integral.
*
       Call CmbnKEr(Array(ipRnr),Array(ipqC),Array(ipDi),
     &              nZeta,la,lb,Zeta,Final,
     &              nComp,Alpha,nAlpha,Beta,nBeta)
*
      Else
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
     &          nComp,Array(ipTxyz))
*
      End If
*
*     Call GetMem(' Exit KnEInt','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_integer_array(iStabM)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
