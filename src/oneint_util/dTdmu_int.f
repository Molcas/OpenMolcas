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
* Copyright (C) 2002, Roland Lindh                                     *
************************************************************************
      SubRoutine dTdmu_int(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of diamagnetic shielding  *
*         integrals.                                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              EFPrm                                                   *
*              Util4                                                   *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics, University      *
*             of Lund, Sweden, September 2002.                         *
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
     &       Array(nZeta*nArr), Ccoor(3,2), TC(3,2)
      Integer lOper(nComp), iStabM(0:nStabM-1), iDCRT(0:7),
     &          iStabO(0:7), iChO(nComp)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 230
      iPrint = nPrint(iRout)
      Call qEnter('dTdmu_int')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In dTdmu_int: Alpha',' ',Alpha,nAlpha,1)
         Call RecPrt(' In dTdmu_int: Beta',' ',Beta,nBeta,1)
      End If
*
      nip = 1
      ipS1 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*3
      ipS2 = nip
      If (lb.ge.1) nip = nip + nZeta*nElem(la)*nElem(lb-1)*3
      ipRes = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      ipB = nip
      nip = nip + nZeta
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,'dTdmu_int: nip-1.gt.nZeta*nArr')
         Write (6,*) 'nip=',nip
         Write (6,*) 'nZeta,nArr=',nZeta,nArr
         Call Abend()
      End If
      ipArr = nip
      mArr = nZeta*nArr - nip + 1
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
*
      ipOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
      End Do
*
      iComp = 1
      llOper = lOper(1)
      Do iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
      End Do
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
*
      Do lDCRT = 0, nDCRT-1
         TC(1,1) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1,1)
         TC(2,1) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2,1)
         TC(3,1) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3,1)
         TC(1,2) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1,2)
         TC(2,2) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2,2)
         TC(3,2) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3,2)
*
*-------Compute contribution from a,b+1
*
         Call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &               Array(ipS1),nZeta,nComp,la,lb+1,A,RB,nRys,
     &               Array(ipArr),mArr,TC,nOrdOp)
*
*--------Compute contribution from a,b-1
*
         If (lb.ge.1)
     &      Call EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &                 Array(ipS2),nZeta,nComp,la,lb-1,A,RB,nRys,
     &                 Array(ipArr),mArr,TC,nOrdOp)
*
*--------Assemble final integral from the derivative integrals
*
         Call Assemble_dTdmu(nZeta,Array(ipRes),la,lb,
     &                       Array(ipS1),Array(ipS2),Array(ipB))
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
*     Call GetMem(' Exit dTdmu_int','LIST','REAL',iDum,iDum)
      Call qExit('dTdmu_int')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
