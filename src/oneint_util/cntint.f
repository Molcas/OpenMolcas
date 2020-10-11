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
* Copyright (C) 1991,2008, Roland Lindh                                *
************************************************************************
      SubRoutine CntInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: to compute contact integrals.                                *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden, February '91                 *
*             Modified from D1Int January 2008.                        *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables

      Character*80 Label
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 150
      iPrint = nPrint(iRout)
*
      Call FZero(Final,nZeta*((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*nIC)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipArr=nip
      na=(la+1)*(la+2)/2
      nb=(lb+1)*(lb+2)/2
      mArr=na*nb
      nip = nip + nZeta*mArr
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'CntInt: nip-1.gt.nArr*nZeta')
         Write (6,*) 'nip=',nip
         Write (6,*) 'nArr,nZeta=',nArr,nZeta
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In CntInt: A',' ',A,1,3)
         Call RecPrt(' In CntInt: RB',' ',RB,1,3)
         Call RecPrt(' In CntInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In CntInt: P',' ',P,nZeta,3)
         Write (6,*) ' In CntInt: la,lb=',la,lb
      End If
*
*     Compute the contact terms.
*
      Call Contact(Zeta,P,nZeta,
     &             A,Array(ipAxyz),la,
     &             RB,Array(ipBxyz),lb,
     &             Ccoor,lOper,iCho,nIC,
     &             Array(ipArr),mArr,
     &             Final,iStabM,nStabM,nComp,rKappa)
*
      If (iPrint.ge.99) Then
         Do iIC = 1, nIC
         Do ia = 1, nElem(la)
            Do ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &               'Contact term(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,iIC),1,nZeta)
            End Do
         End Do
         End Do
      End If
*
*     Call GetMem(' Exit CntInt','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_integer(nOrdOp)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
