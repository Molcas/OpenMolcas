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
     &                     nCmp,nOp,SO_tmp,nDeg,iAO)
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: MolWgh
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 AOValue(mAO*nCoor,mBas_Eff,nCmp),
     &       SO_tmp(mAO*nCoor,mBas,nCmp*nDeg)
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

      iOff=1
      iAdd=mBas-mBas_Eff
      Do i1 = 1, nCmp
         Do j1 = 0, nIrrep-1
            iSO=iAOtSO(iAO+i1,j1)

            If (iSO<0) Cycle
            xa= DBLE(iChTbl(j1,nOp))
            Do i2 = 1, mBas_Eff
               Call DaXpY_(mAO*nCoor,Fact*xa,
     &                     AOValue(:,i2,i1),1,
     &                     SO_tmp(:,i2+iAdd,iOff),1)
            End Do
            iOff = iOff + 1
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
