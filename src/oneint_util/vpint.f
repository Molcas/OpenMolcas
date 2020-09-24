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
* Copyright (C) 1993, Bernd Artur Hess                                 *
************************************************************************
      SubRoutine VPInt(
#define _CALLING_
#include "int_interface.fh"
     &                )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of  pV integrals          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : qEnter                                                  *
*              RecPrt                                                  *
*              Util1                                                   *
*              DCopy  (ESSL)                                           *
*              NSOInt                                                  *
*              GetMem                                                  *
*              qExit                                                   *
*                                                                      *
*     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
*             Chemie, University of Bonn, Germany, April 1993          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake, XCff2D, XRys2D
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2
*
      iRout = 221
      iPrint = nPrint(iRout)
*
c     Call qEnter('vpint')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In vpint: Alpha','(5D20.13)',Alpha,nAlpha,1)
         Call RecPrt(' In vpint: Beta','(5D20.13)',Beta,nBeta,1)
      End If
*
      nRys=nHer
*
      nip = 1
      ipB = nip
      nip = nip + nZeta
      ipS1 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)
      If (lb.gt.0) Then
         ipS2 = nip
         nip = nip + nZeta*nElem(la)*nElem(lb-1)
      Else
         ipS2=ipS1
      End If
      ipArr = nip
      mArr = nArr - (nip-1)/nZeta
      If (mArr.lt.0) Then
         Call WarningMessage(2,'VpInt: mArr<0!')
         Call Abend()
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
      call dcopy_(nZeta*nArr,[Zero],0,Array,1)
*     Compute contribution from a,b+1
*
      kRys = ((la+1)+lb+2)/2
*
      kIC=1
      kComp=1
      Call NAInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &           Array(ipS1),nZeta,nIC,nComp,la,lb+1,A,RB,kRys,
     &           Array(ipArr),mArr,CCoor,nOrdOp,lOper,iChO,
     &           iStabM,nStabM,
     &           PtChrg,nGrid,iAddPot)

      ipOff = ipB
      Do 100 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
100   Continue
*
*     Compute contribution from a,b-1
*
      If (lb.gt.0) Then
         kRys = ((la-1)+lb+2)/2
*
         Call NAInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,nRys,
     &              Array(ipArr),mArr,CCoor,nOrdOp,lOper,iChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot)
      End If
*
*     Assemble final integral from the derivative integrals
*
      If (iPrint.ge.99) Call RecPrt(' In vpint: Beta (expanded)',
     &                  '(5D20.13)',Array(ipB),nZeta,1)
*
      Call Util8(Array(ipB),nZeta,Final,la,lb,Array(ipS1),Array(ipS2))
*
      If (iPrint.ge.49) Then
         Do i=1,3
            Call RecPrt('VpInt: Final',' ',Final(1,1,1,i),nZeta,
     &                 nElem(la)*nElem(lb))
         End Do
      End If
c     Call qExit('vpint')
      Return
      End
