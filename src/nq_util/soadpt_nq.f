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
     &                     nCmp,nOp,SOValue,nMOs,nDeg,iAO)
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: MolWgh, nBas
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 AOValue(mAO*nCoor,mBas_Eff,nCmp),
     &       SOValue(mAO*nCoor,nMOs)
      Integer   iOff_MO(0:7)
#ifdef _DEBUGPRINT_
      Character*80 Label
#endif
*
      If (MolWgh.eq.0) Then
         Fact=One/DBLE(nDeg)
      Else If (MolWgh.eq.1) Then
         Fact=One
      Else
         Fact=One/Sqrt(DBLE(nDeg))
      End If
*
*---- Compute some offsets
*
      itmp1=1
      Do iIrrep = 0, nIrrep-1
         iOff_MO(iIrrep)=itmp1
         itmp1=itmp1+nBas(iIrrep)
      End Do

      iAdd=mBas-mBas_Eff
      Do i1 = 1, nCmp
         Do i2 = 1, mBas_Eff
            Do iIrrep = 0, nIrrep-1
               iSO0=iAOtSO(iAO+i1,iIrrep)
*              Write (6,*) 'iAO,i1,i2,iIrrep:',iAO,i1,i2,iIrrep
*              Write (6,*) 'iSO0:',iSO0
               If (iSO0<0) Cycle
               iMO =iOff_MO(iIrrep)
               xa= DBLE(iChTbl(iIrrep,nOp))
               iSO = iSO0 + i2 - 1
               iSO1=iMO+iSO-1+iAdd
*              Write (6,*) 'Fact*xa:',Fact*xa
*              Write (6,*) 'iMO,iSO,iAdd=',iMO,iSO,iAdd
*              Write (6,*) 'iSO1=',iSO1
*              Call RecPrt('TabAO',' ',AOValue(:,i2,i1),mAO,nCoor)
*              Call RecPrt('TabSO',' ',SOValue(:,iSO1),mAO,nCoor)

               Call DaXpY_(mAO*nCoor,Fact*xa,
     &                     AOValue(:,i2,i1),1,
     &                     SOValue(:,iSO1),1)
*              Call RecPrt('TabSO',' ',SOValue(:,iSO1),mAO,nCoor)
            End Do
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Do iCmp = 1, nCmp*nDeg
         Write (Label,'(A,I2,A)') 'SO_tmp(mAO,nCoor,mBas,',iCmp,')'
         Call RecPrt(Label,' ',SO_tmp(1,1,1,iCmp),mAO*nCoor,mBas)
      End Do
#endif
*
      Return
      End
