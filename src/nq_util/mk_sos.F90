!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

subroutine mk_SOs(TabSO,mAO,mGrid,nMOs,List_s,List_Bas,nList_s,jList_s)

use iSD_data
use Center_Info
use Symmetry_Info, only: nIrrep, iChTbl
use SOAO_Info, only: iAOtSO
use Basis_Info, only: MolWgh, nBas
use nq_Grid, only: iBfn_Index, TabAO

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 TabSO(mAO*mGrid,nMOs)
integer :: list_s(2,nList_s), list_bas(2,nlist_s)
integer iOff_MO(0:7)
integer :: jList_s

! Compute some offsets

itmp1 = 1
do iIrrep=0,nIrrep-1
  iOff_MO(iIrrep) = itmp1
  itmp1 = itmp1+nBas(iIrrep)
end do

nBfn = size(iBfn_Index,2)
do iBfn=1,nBfn
  ilist_s = iBfn_Index(2,iBfn)
  if ((jlist_s /= 0) .and. (ilist_s /= jlist_s)) cycle
  i1 = iBfn_Index(3,iBfn)
  i2 = iBfn_Index(4,iBfn)
  iSh = list_s(1,ilist_s)
  kDCRE = list_s(2,ilist_s)
  mBas_Eff = List_Bas(1,ilist_s)
  mBas = iSD(3,iSh)
  iAO = iSD(7,iSh)
  mdci = iSD(10,iSh)
  nDeg = nIrrep/dc(mdci)%nStab
  nOp = NrOpr(kDCRE)

  if (MolWgh == 0) then
    Fact = One/dble(nDeg)
  else if (MolWgh == 1) then
    Fact = One
  else
    Fact = One/sqrt(dble(nDeg))
  end if

  iAdd = mBas-mBas_Eff
  do iIrrep=0,nIrrep-1
    iSO0 = iAOtSO(iAO+i1,iIrrep)
    if (iSO0 < 0) cycle

    iMO = iOff_MO(iIrrep)

    xa = dble(iChTbl(iIrrep,nOp))
    iSO = iSO0+i2-1
    iSO1 = iMO+iSO-1+iAdd
    call DaXpY_(mAO*mGrid,Fact*xa,TabAO(:,:,iBfn),1,TabSO(:,iSO1),1)
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('mk_SOs: TabSO',' ',TabSO,mAO*mGrid,nMOs)
#endif

return

end subroutine mk_SOs
