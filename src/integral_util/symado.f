************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine SymAdO(ArrIn,nZeta,la,lb,nComp,ArrOut,nIC,iDCRT,
     &                  lOper,iChO,Factor)
      use Symmetry_Info, only: iChTbl, iOper
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Real*8 ArrIn (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       ArrOut(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Prmt(0:7)
      Integer iDCRT ,iTwoj(0:7), lOper(nComp), iChO(nComp)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 200
      iPrint = nPrint(iRout)
C     Call qEnter('SymAdO')
C     nA = (la+1)*(la+2)/2
C     nB = (lb+1)*(lb+2)/2
C     Call RecPrt('SymAdO: ArrIn',' ',ArrIn,nZeta*nA*nB, nComp)
*
*--------Accumulate contributions
*
      iIC = 0
      Do 103 iComp = 1, nComp
         pO = xPrmt(iOper(iDCRT),iChO(iComp))
         Do 104 iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).eq.0) Go To 104
            iIC = iIC + 1
            Xg = DBLE(iChTbl(iIrrep,iDCRT))
            Call DaXpY_(nZeta*nElem(la)*nElem(lb),Xg*pO*Factor,
     &                 ArrIn(1,1,1,iComp),1,ArrOut(1,1,1,iIC),1)
 104     Continue
 103  Continue
      If (iIC.ne.nIC) Then
         Call WarningMessage(2,' Abend in SymAdO: iIC.ne.nIC')
         Write (6,*) 'iIC,nIC=',iIC,nIC
         Call Abend()
      End If
C     Call RecPrt('SymAdO: ArrOut',' ',ArrOut,nZeta*nA*nB, nIC)
*
*     Call GetMem(' Exit SymAdO','LIST','REAL',iDum,iDum)
C     Call qExit('SymAdO')
      Return
      End
