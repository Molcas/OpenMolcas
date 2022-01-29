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
     &                     nCmp,nOp,SOValue,nDeg,iAO)
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: MolWgh
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 AOValue(mAO,nCoor,mBas_Eff,nCmp),
     &       SOValue(mAO,nCoor,mBas,nCmp*nDeg)
      Character*80 Label
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
         Do j1 = 0, nIrrep-1
            If (iAOtSO(iAO+i1,j1)>0) Then
               xa= DBLE(iChTbl(j1,nOp))
               Call DaXpY_(mAO*nCoor*mBas_Eff,Fact*xa,
     &                   AOValue(:,:,:,i1),1,
     &                   SOValue(1,1,1+iAdd,iSO),1)
               iSO = iSO + 1
            End If
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Do iCmp = 1, nCmp*nDeg
         Write (Label,'(A,I2,A)') 'SOValue(mAO,nCoor,mBas,',iCmp,')'
         Call RecPrt(Label,' ',SOValue(1,1,1,iCmp),mAO*nCoor,mBas)
      End Do
#endif
*
      Return
      End
