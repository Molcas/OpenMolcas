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
*               2016, Lasse Kragh Soerensen                            *
************************************************************************
      SubRoutine OAMInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of orbital angular        *
*         momentum integrals.                                          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              MltPrm                                                  *
*              Util2                                                   *
*              DCopy   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, February '91                            *
*             Placed restrictions on the differentiation. Lasse '16    *
************************************************************************
      use Logical_Info, only: lUPONLY, lDOWNONLY
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Real*8 TC(3)
      Integer  iStabO(0:7), iDCRT(0:7)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 210
      iPrint = nPrint(iRout)
*     Call qEnter('OAMInt')
*
      nip = 1
      ipB = nip
      nip = nip + nZeta
      ipS1 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*3
      ipS2 =  1
      If (lb.gt.0) Then
         ipS2 = nip
         nip = nip + nZeta*nElem(la)*nElem(lb-1)*3
      End If
      ipRes = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,' OAMInt: nip-1.gt.nZeta*nArr')
         Call Abend()
      End If
      ipArr = nip
      mArr = (nArr*nZeta - (nip-1))/nZeta
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
      iComp=3
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
      ipOff = ipB
      Do 100 iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
 100  Continue
*
      Do 102 lDCRT = 0, nDCRT-1
         Call OA(iDCRT(lDCRT),CCoor,TC)
*
         IF (.NOT.lDOWNONLY) THEN
            nHer = (la + (lb+1) + (nOrdOp-1) + 2) / 2
            Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &                Array(ipS1),nZeta,iComp,la,lb+1,A,RB,nHer,
     &                Array(ipArr),mArr,TC,nOrdOp-1)
         END IF
*
         If (lb.gt.0) Then
            IF (.NOT.lUPONLY) THEN
               nHer = (la + (lb-1) + (nOrdOp-1) + 2) / 2
               Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &                   Array(ipS2),nZeta,iComp,la,lb-1,A,RB,nHer,
     &                   Array(ipArr),mArr,TC,nOrdOp-1)
            END IF
         End If
*
*        Combine derivatives of dipole integrals to generate the
*        orbital angular momentum integrals.
*
         Call Util2(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipS1),
     &                                           Array(ipS2))
*
*--------Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 102  Continue
*
*     Call qExit('OAMInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
