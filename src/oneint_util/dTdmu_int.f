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
      SubRoutine dTdmu_int(
#define _CALLING_
#include "int_interface.fh"
     &                    )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of diamagnetic shielding  *
*         integrals.                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics, University      *
*             of Lund, Sweden, September 2002.                         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Real*8 TC(3,2)
      Integer iDCRT(0:7), iStabO(0:7)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 230
      iPrint = nPrint(iRout)
*
      nRys=nHer
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
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
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
      Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
      Do lDCRT = 0, nDCRT-1
         Call OA(iDCRT(lDCRT),Ccoor(1:3,1),TC(1:3,1))
         Call OA(iDCRT(lDCRT),Ccoor(1:3,2),TC(1:3,2))
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
         nOp = NrOpr(iDCRT(lDCRT))
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
*     Call GetMem(' Exit dTdmu_int','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
