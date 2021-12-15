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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine VeInt(
#define _CALLING_
#include "int_interface.fh"
     &                )
************************************************************************
*                                                                      *
* Object: to compute the velocity integrals with the Gauss-Hermite     *
*         quadrature.                                                  *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden, January '91                  *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables.
      Logical ABeq(3)
      Integer iStabO(0:7), iDCRT(0:7)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 195
      iPrint = nPrint(iRout)
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+2)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+2)
      ipVxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)
      ipB = nip
      nip = nip + nZeta
      ipRes = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'VeInt: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in VeInt'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In VeInt: A',' ',A,1,3)
         Call RecPrt(' In VeInt: RB',' ',RB,1,3)
         Call RecPrt(' In VeInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In VeInt: P',' ',P,nZeta,3)
         Write (6,*) ' In VeInt: la,lb=',la,lb
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*     Compute the cartesian values of the basis functions angular part
*
      Call CrtCmp(Zeta,P,nZeta,A,Array(ipAxyz),
     &               la,HerR(iHerR(nHer)),nHer,ABeq)
      Call CrtCmp(Zeta,P,nZeta,RB,Array(ipBxyz),
     &               lb+1,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the contribution from the multipole moment operator
*
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call CrtCmp(Zeta,P,nZeta,Ccoor,Array(ipRxyz),
     &            0,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call Assmbl(Array(ipQxyz),
     &             Array(ipAxyz),la,
     &             Array(ipRxyz),0,
     &             Array(ipBxyz),lb+1,
     &             nZeta,HerW(iHerW(nHer)),nHer)
*
*     Compute the cartesian components for the velocity integrals.
*     The velocity components are linear combinations of overlap
*     components.
*
      ipBOff = ipB
      Do 210 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
         ipBOff = ipBOff + 1
 210  Continue
*
      Call VelInt(Array(ipVxyz),Array(ipQxyz),la,lb,
     &          Array(ipB),nZeta)
*
*     Combine the cartesian components to the full one electron
*     integral.
*
      Call CmbnVe(Array(ipQxyz),nZeta,la,lb,0,Zeta,rKappa,Array(ipRes),
     &          nComp,Array(ipVxyz))
*
      llOper=lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
      Do 102 lDCRT = 0, nDCRT-1
*
*--------Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 102  Continue
*
*     Call GetMem(' Exit VeInt','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(ZInv)
         Call Unused_integer(nOrdOp)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
