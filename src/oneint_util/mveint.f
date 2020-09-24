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
      SubRoutine MVeInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: to compute the mass-velocity integrals with the Gauss-       *
*         Hermite quadrature.                                          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              GetMem                                                  *
*              MVe                                                     *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden. February '91                 *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Logical ABeq(3)
      Character*80 Label
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 190
      iPrint = nPrint(iRout)
      Call qEnter('MVeInt')
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
      nip = nip + nZeta*3*nHer*(nOrdOp-3)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+3)*(lb+3)*(nOrdOp-3)
      iprV2  = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)* 2
      iprV4 = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)
      ipA = nip
      nip = nip + nZeta
      ipB = nip
      nip = nip + nZeta
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'MVeInt: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in MVeInt'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In MVeInt: A',' ',A,1,3)
         Call RecPrt(' In MVeInt: RB',' ',RB,1,3)
         Call RecPrt(' In MVeInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In MVeInt: P',' ',P,nZeta,3)
         Call RecPrt(' In MVeInt: Zeta',' ',Zeta,nZeta,1)
         Call RecPrt(' In MVeInt: Roots',' ',
     &               HerR(iHerR(nHer)),nHer,1)
         Call GetMem(' In MVeInt','LIST','REAL',iDum,iDum)
         Write (6,*) ' In MVeInt: la,lb=',la,lb
      End If
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
     &            nOrdOp-4,HerR(iHerR(nHer)),nHer,ABeq)
*
*     Compute the cartesian components for the multipole moment
*     integrals. The integrals are factorized into components.
*
       Call Assmbl(Array(ipQxyz),
     &             Array(ipAxyz),la+2,
     &             Array(ipRxyz),nOrdOp-4,
     &             Array(ipBxyz),lb+2,
     &             nZeta,HerW(iHerW(nHer)),nHer)
*
*     Compute the cartesian components for the mass-velocity integrals.
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
      Call MVe(Array(iprV2),Array(iprV4), Array(ipQxyz),la,lb,
     &          Array(ipA),Array(ipB),nZeta)
*
*     Combine the cartesian components to the full one electron
*     integral.
*
      Call CmbnMV(Array(ipQxyz),nZeta,la,lb,nOrdOp-4,Zeta,rKappa,Final,
     &          nComp,Array(iprV2),Array(iprV4))
*
      If (iPrint.ge.99) Then
         Do 300 ia = 1, nElem(la)
            Do 310 ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &               'Mass-Velocity(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,1,ia,ib),nZeta,nComp)
 310        Continue
 300     Continue
      End If
*
*     Call GetMem(' Exit MVeInt','LIST','REAL',iDum,iDum)
      Call qExit('MVeInt')
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
