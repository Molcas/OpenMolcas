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
! Copyright (C) 1991,2021, Roland Lindh                                *
!***********************************************************************

subroutine SymAdp_Full(SOIntegrals,nSOInt,list_s,nlist_s,Fact,ndc,nD)
!***********************************************************************
!                                                                      *
! Object: to transform the one-electron matrix elements from AO basis  *
!         to SO basis.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1991                                             *
!***********************************************************************

use iSD_data, only: iSD
use Basis_Info, only: nBas
use Symmetry_Info, only: iChTbl, nIrrep
use SOAO_Info, only: iAOtSO
use nq_Grid, only: Dens_AO, iBfn_Index
use Index_Functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSOInt, nlist_s, list_s(2,nlist_s), ndc, nD
real(kind=wp), intent(inout) :: SOIntegrals(nSOInt,nD)
real(kind=wp), intent(in) :: Fact(ndc,ndc)
integer(kind=iwp) :: iAO, iBfn, iBfn_, iCmp, ilist_s, indAO1, Indij, iPnt, iShell, iSkal, iSO, iSO1, j1, jBfn, jBfn_, jlist_s, &
                     jShell, jSkal, jSO, kDCRE, kDCRR, loper, mBfn, mdci, mdcj, nBfn, nOp(2)
real(kind=wp) :: xa, xaxb, xb
integer(kind=iwp), allocatable :: BasList(:,:)
integer(kind=iwp), external :: iPntSO, NrOpr

!                                                                      *
!***********************************************************************
!                                                                      *
nBfn = size(iBfn_Index,2)
call mma_Allocate(BasList,2,nBfn,Label='BasList')
loper = 1
do j1=0,nIrrep-1
  iPnt = iPntSO(j1,j1,lOper,nbas)

  ! Pick up only basis functions which contribute to (j1,j1)
  mBfn = 0
  do iBfn=1,nBfn
    ilist_s = iBfn_Index(2,iBfn)
    iCmp = iBfn_Index(3,iBfn)
    indAO1 = iBfn_Index(6,iBfn)
    iSkal = list_s(1,ilist_s)
    iAO = iSD(7,iSkal)
    iSO1 = iAOtSO(iAO+iCmp,j1)
    if (iSO1 < 0) cycle
    mBfn = mBfn+1
    BasList(1,mBfn) = iBfn
    BasList(2,mBfn) = iSO1+IndAO1-1
  end do

  do iBfn_=1,mBfn
    iBfn = BasList(1,iBfn_)
    iSO = BasList(2,iBfn_)

    ilist_s = iBfn_Index(2,iBfn)
    iSkal = list_s(1,ilist_s)
    kDCRE = list_s(2,ilist_s)
    mdci = iSD(10,iSkal)
    iShell = iSD(11,iSkal)
    nOp(1) = NrOpr(kDCRE)
    xa = real(iChTbl(j1,nOp(1)),kind=wp)

    do jBfn_=1,iBfn_
      jBfn = BasList(1,jBfn_)
      jSO = BasList(2,jBfn_)

      jlist_s = iBfn_Index(2,jBfn)
      jSkal = list_s(1,jlist_s)
      kDCRR = list_s(2,jlist_s)
      mdcj = iSD(10,jSkal)
      jShell = iSD(11,jSkal)
      nOp(2) = NrOpr(kDCRR)
      xb = real(iChTbl(j1,nOp(2)),kind=wp)

      xaxb = xa*xb
      if ((iShell == jShell) .and. (nOp(1) /= nOp(2)) .and. (iSO == jSO)) xaxb = xaxb*Two

      Indij = iPnt+iTri(iSO,jSO)

      SOIntegrals(Indij,:) = SOIntegrals(Indij,:)+Fact(mdci,mdcj)*xaxb*Dens_AO(iBfn,jBfn,:)

    end do ! jBfn
  end do   ! iBfn
end do     ! j1
call mma_deAllocate(BasList)

return

end subroutine SymAdp_Full
