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
* Copyright (C) 2022, Roland Lindh                                     *
************************************************************************
      Subroutine mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s)
      use iSD_data
      use Center_Info
      use Symmetry_Info, only: nIrrep, iChTbl
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: MolWgh, nBas
      use nq_Grid, only: iBfn_Index, TabAO
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 TabSO(mAO*mGrid,nMOs)
      Integer :: list_s(2,nList_s), list_bas(2,nlist_s)
      Integer   iOff_MO(0:7)
*
*---- Compute some offsets
*
      itmp1=1
      Do iIrrep = 0, nIrrep-1
         iOff_MO(iIrrep)=itmp1
         itmp1=itmp1+nBas(iIrrep)
      End Do

      nBfn=Size(iBfn_Index,2)
      Do iBfn = 1, nBfn
         ilist_s=iBfn_Index(2,iBfn)
         i1     =iBfn_Index(3,iBfn)
         i2     =iBfn_Index(4,iBfn)
         iSh    =list_s(1,ilist_s)
         kDCRE = list_s(2,ilist_s)
         mBas_Eff=List_Bas(1,ilist_s)
         mBas   =iSD( 3,iSh)
         iAO    =iSD( 7,iSh)
         mdci   =iSD(10,iSh)
         nDeg   =nIrrep/dc(mdci)%nStab
         nOp = NrOpr(kDCRE)

         If (MolWgh.eq.0) Then
            Fact=One/DBLE(nDeg)
         Else If (MolWgh.eq.1) Then
            Fact=One
         Else
            Fact=One/Sqrt(DBLE(nDeg))
         End If

         iAdd=mBas-mBas_Eff
         Do iIrrep = 0, nIrrep-1
            iSO0=iAOtSO(iAO+i1,iIrrep)
            If (iSO0<0) Cycle

            iMO=iOff_MO(iIrrep)

            xa= DBLE(iChTbl(iIrrep,nOp))
            iSO = iSO0 + i2 - 1
            iSO1=iMO+iSO-1+iAdd
            Call DaXpY_(mAO*mGrid,Fact*xa,
     &                  TabAO(:,:,iBfn),1,
     &                  TabSO(:,iSO1),1)
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('mk_SOs: TabSO',' ',TabSO,mAO*mGrid,nMOs)
#endif
*
      Return
      End Subroutine mk_SOs
