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
      SubRoutine D1Int(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: to compute the 1-electron Darwin contact term.               *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              Darwin                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden, February '91                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Character*80 Label
      Integer lOper(nComp), iStabM(0:nStabM-1), iChO(nComp)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 150
      iPrint = nPrint(iRout)
      Call qEnter('D1Int')
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'D1Int: nip-1.gt.nArr*nZeta')
         Write (6,*) 'nip=',nip
         Write (6,*) 'nArr,nZeta=',nArr,nZeta
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In D1Int: A',' ',A,1,3)
         Call RecPrt(' In D1Int: RB',' ',RB,1,3)
         Call RecPrt(' In D1Int: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In D1Int: P',' ',P,nZeta,3)
         Write (6,*) ' In D1Int: la,lb=',la,lb
      End If
*
*     Compute the contact terms.
*
      Call Darwin(Zeta,P,nZeta,A,Array(ipAxyz),la,
     &                       RB,Array(ipBxyz),lb,
     &                       Final,iStabM,nStabM,nComp,rKappa)
*
      If (iPrint.ge.99) Then
         Do 300 ia = 1, nElem(la)
            Do 310 ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &               'Darwin contact(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,1,ia,ib),nZeta,nComp)
 310        Continue
 300     Continue
      End If
*
*     Call GetMem(' Exit D1Int','LIST','REAL',iDum,iDum)
      Call qExit('D1Int')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
