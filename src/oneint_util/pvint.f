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
      SubRoutine PVInt(
#define _CALLING_
#include "int_interface.fh"
     &                , Kernel)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of  pX integrals          *
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
      External Kernel
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"

#include "int_interface.fh"
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2
*                                                                      *
************************************************************************
*                                                                      *
*      Interface
*      Subroutine Kernel(
*#define _CALLING_
*#include "int_interface.fh"
*     &                 )
*#include "int_interface.fh"
*      End Subroutine Kernel
*      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 221
      iPrint = nPrint(iRout)
      Call qEnter('pvint')
*
      If (iPrint.ge.99) Then
         Write (6,*) 'PVInt: nIC,nComp=',nIC,nComp
         Call RecPrt(' In pvint: Alpha','(5D20.13)',Alpha,nAlpha,1)
         Call RecPrt(' In pvint: Beta','(5D20.13)',Beta,nBeta,1)
      End If
*
      nip = 1
      ipA = nip
      nip = nip + nZeta
      ipS1 = nip
      nip = nip + nZeta*nElem(la+1)*nElem(lb)*nIC
      ipS2 =  1
      If (la.gt.0) Then
         ipS2 = nip
         nip = nip + nZeta*nElem(la-1)*nElem(lb)*nIC
      Else
         ipS2=ipS1
      End If
      ipArr = nip
      mArr = nArr-(nip-1)/nZeta
      If (mArr.lt.0) Then
         Call WarningMessage(2,'pVInt: mArr<0!')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Compute contribution from a+1,b
*
      kRys = ((la+1)+lb+2)/2
      Call Kernel(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &           Array(ipS1),nZeta,nIC,nComp,la+1,lb,A,RB,kRys,
     &           Array(ipArr),mArr,CCoor,nOrdOp,lOper,iChO,
     &           iStabM,nStabM,
     &           PtChrg,nGrid,iAddPot)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute contribution from a-1,b
*
      If (la.gt.0) Then
         kRys = ((la-1)+lb+2)/2
         Call Kernel(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Array(ipS2),nZeta,nIC,nComp,la-1,lb,A,RB,kRys,
     &              Array(ipArr),mArr,CCoor,nOrdOp,lOper,iChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ipOff = ipA
      Do iBeta = 1, nBeta
         call dcopy_(nAlpha,Alpha,1,Array(ipOff),1)
         ipOff = ipOff + nAlpha
      End Do
      If (iPrint.ge.99) Then
         Call RecPrt(' In pvint: Alpha (expanded)','(5D20.13)',
     &         Array(ipA),nZeta,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Assemble final integral from the derivative integrals
*
      Call Ass_pX(Array(ipA),nZeta,Final,la,lb,Array(ipS1),Array(ipS2),
     &           nIC)
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.49) Then
         Do i=1,3
           Call RecPrt('pVInt: Final',' ',Final(1,1,1,i),
     &                  nZeta,nElem(la)*nElem(lb))
         End Do
      End If
      Call qExit('pvint')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nHer)
      End
