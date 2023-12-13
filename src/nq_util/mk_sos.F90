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

use iSD_data, only: iSD
use Center_Info, only: dc
use Symmetry_Info, only: iChTbl, nIrrep
use SOAO_Info, only: iAOtSO
use Basis_Info, only: MolWgh, nBas
use nq_Grid, only: iBfn_Index, TabAO
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, mGrid, nMOs, nlist_s, list_s(2,nlist_s), list_bas(2,nlist_s), jList_s
real(kind=wp), intent(inout) :: TabSO(mAO,mGrid,nMOs)
integer(kind=iwp) :: i1, i2, iAdd, iAO, iBfn, iIrrep, ilist_s, iMO, iOff_MO(0:7), iSh, iSO, iSO0, iSO1, itmp1, kDCRE, mBas, &
                     mBas_Eff, mdci, nBfn, nDeg, nOp
real(kind=wp) :: Fact, xa
integer(kind=iwp), external :: NrOpr

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
    Fact = One/real(nDeg,kind=wp)
  else if (MolWgh == 1) then
    Fact = One
  else
    Fact = One/sqrt(real(nDeg,kind=wp))
  end if

  iAdd = mBas-mBas_Eff
  do iIrrep=0,nIrrep-1
    iSO0 = iAOtSO(iAO+i1,iIrrep)
    if (iSO0 < 0) cycle

    iMO = iOff_MO(iIrrep)

    xa = real(iChTbl(iIrrep,nOp),kind=wp)
    iSO = iSO0+i2-1
    iSO1 = iMO+iSO-1+iAdd
    TabSO(:,:,iSO1) = TabSO(:,:,iSO1)+Fact*xa*TabAO(:,:,iBfn)
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('mk_SOs: TabSO',' ',TabSO,mAO*mGrid,nMOs)
#endif

return

end subroutine mk_SOs
