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
      Subroutine SOAdpt_NQ(AOValue,mAO,nCoor,mBas,mBas_Eff,
     &                     nCmp,nOp,SOValue,nDeg,IndShl)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
      Real*8 AOValue(mAO,nCoor,mBas_Eff,nCmp),
     &       SOValue(mAO,nCoor,mBas,nCmp*nDeg), Aux(8)
      Integer   iTwoj(0:7)
      Character*80 Label
      Data iTwoj/1,2,4,8,16,32,64,128/
*
      Call QEnter('SOAdpt')
      iRout=133
      iPrint=nPrint(iRout)
*     Call GetMem('SOAdpt_E','CHEC','REAL',iDum,iDum)
*
      If (MolWgh.eq.0) Then
         Fact=One/DBLE(nDeg)
      Else If (MolWgh.eq.1) Then
         Fact=One
      Else
         Fact=One/Sqrt(DBLE(nDeg))
      End If
      iSO=1
      iAdd=mBas-mBas_Eff
      Do i1 = 1, nCmp
         iaux=0
         Do j1 = 0, nIrrep-1
            If (iAnd(IrrCmp(IndShl+i1),iTwoj(j1)).ne.0) Then
               iaux=iaux+1
               xa= DBLE(iChTbl(j1,nOp))
               Aux(iAux)=Fact*xa
            End If
         End Do
         If (iPrint.ge.49) Call RecPrt('Aux',' ',Aux,1,iAux)
         Call DnaXpY(iAux,mAO*nCoor*mBas_Eff,Aux,1,
     &               AOValue(1,1,1,i1),1,0,
     &               SOValue(1,1,1+iAdd,iSO),1,mAO*nCoor*mBas)
         iSO=iSO+iAux
      End Do
*
      If (iPrint.ge.49) Then
         Do iCmp = 1, nCmp*nDeg
            Write (Label,'(A,I2,A)') 'SOValue(mAO,nCoor,mBas,',iCmp,')'
            Call RecPrt(Label,' ',SOValue(1,1,1,iCmp),mAO*nCoor,mBas)
         End Do
      End If
*
*     Call GetMem('SOAdpt_X','CHEC','REAL',iDum,iDum)
      Call QExit ('SOAdpt')
      Return
      End
