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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine ChoMP2_Setup_Index(iFirst,iFirstS,NumOcc,LnOcc,NumBatOrb,LnBatOrb,LnT1am,LiT1am,LnPQprod,LiPQprod,LnMatij,LiMatij,mSym, &
                              mBatch)
!
! Thomas Bondo Pedersen, Nov. 2004 / Feb. 2005.
!
! Purpose: set local index arrays and counters.

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Cholesky, only: nSym
use ChoMP2, only: ChoAlg, iBatOrb, nBatch, nBatOrbT, nDel, nFro, nOcc, nVir
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mBatch, mSym
integer(kind=iwp), intent(out) :: iFirst(mBatch), iFirstS(mSym,mBatch), NumOcc(mBatch), LnOcc(mSym,mBatch), NumBatOrb(mBatch), &
                                  LnBatOrb(mSym,mBatch), LnT1am(mSym,mBatch), LiT1am(mSym,mSym,mBatch), LnPQprod(mSym,mBatch), &
                                  LiPQprod(mSym,mSym,mBatch), LnMatij(mSym,mBatch), LiMatij(mSym,mSym,mBatch)
integer(kind=iwp) :: i, iBatch, iSym, iSyma, iSymi, iSymj, Left, Num
integer(kind=iwp), external :: Cho_iRange

if (mBatch /= nBatch) call SysAbendMsg('ChoMP2_Setup_Index','mBatch !=  nBatch','Error')
if (mSym /= nSym) call SysAbendMsg('ChoMP2_Setup_Index','mSym !=  nSym','Error')

iFirst(:) = 0
iFirstS(:,:) = 0
NumOcc(:) = 0
NumBatOrb(:) = 0
LnOcc(:,:) = 0
LnBatOrb(:,:) = 0
LnT1am(:,:) = 0
LiT1am(:,:,:) = 0
if (.false.) then
  LnPQprod(:,:) = 0
  LiPQprod(:,:,:) = 0
end if
if (ChoAlg == 2) then
  LnMatij(:,:) = 0
  LiMatij(:,:,:) = 0
end if

Num = nBatOrbT/nBatch

! I am not sure if NumOcc is used somewhere else so I will
! define it as before even if Im using NumInBat for setting up
! indices in this routine. //Jonas
do iBatch=1,nBatch
  if (.false.) then
    NumBatOrb(iBatch) = Num
  else
    NumOcc(iBatch) = Num
    NumBatOrb(iBatch) = Num
  end if
end do

Left = nBatOrbT-nBatch*Num
do iBatch=nBatch,nBatch-Left+1,-1
  if (.false.) then
    NumBatOrb(iBatch) = NumBatOrb(iBatch)+1
  else
    NumOcc(iBatch) = NumOcc(iBatch)+1
    NumBatOrb(iBatch) = NumBatOrb(iBatch)+1
  end if
end do

iFirst(1) = 1
do i=1,NumBatOrb(1)
  iSym = Cho_iRange(i,iBatOrb,nSym,.false.)
  if (.false.) then
    LnBatOrb(iSym,1) = LnBatOrb(iSym,1)+1
  else
    LnBatOrb(iSym,1) = LnBatOrb(iSym,1)+1
    LnOcc(iSym,1) = LnOcc(iSym,1)+1
  end if
  if (iFirstS(iSym,1) < 1) iFirstS(iSym,1) = i-iBatOrb(iSym)
end do

do iBatch=2,nBatch
  iFirst(iBatch) = iFirst(iBatch-1)+NumBatOrb(iBatch-1)
  do i=iFirst(iBatch),iFirst(iBatch)+NumBatOrb(iBatch)-1
    iSym = Cho_iRange(i,iBatOrb,nSym,.false.)
    if (.false.) then
      LnBatOrb(iSym,iBatch) = LnBatOrb(iSym,iBatch)+1
    else
      LnBatOrb(iSym,iBatch) = LnBatOrb(iSym,iBatch)+1
      LnOcc(iSym,iBatch) = LnOcc(iSym,iBatch)+1
    end if
    if (iFirstS(iSym,iBatch) < 1) iFirstS(iSym,iBatch) = i-iBatOrb(iSym)
  end do
end do

do iBatch=1,nBatch
  do iSym=1,nSym
    do iSymi=1,nSym
      iSyma = Mul(iSymi,iSym)
      LiT1am(iSyma,iSymi,iBatch) = LnT1am(iSym,iBatch)
      LnT1am(iSym,iBatch) = LnT1am(iSym,iBatch)+nVir(iSyma)*LnOcc(iSymi,iBatch)
      if (.false.) then
        LiPQprod(iSyma,iSymi,iBatch) = LnPQprod(iSym,iBatch)
        LnPQprod(iSym,iBatch) = LnPQprod(iSym,iBatch)+(nOcc(iSymA)+nVir(iSymA)+nFro(iSymA)+nDel(iSymA))*(LnBatOrb(iSymI,iBatch))
      end if
    end do
  end do
end do

if (ChoAlg == 2) then
  do iBatch=1,nBatch
    do iSym=1,nSym
      do iSymj=1,nSym
        iSymi = Mul(iSymj,iSym)
        if (iSymi == iSymj) then
          LiMatij(iSymi,iSymi,iBatch) = LnMatij(iSym,iBatch)
          LnMatij(iSym,iBatch) = LnMatij(iSym,iBatch)+nTri_Elem(LnOcc(iSymi,iBatch))
        else if (iSymi < iSymj) then
          LiMatij(iSymi,iSymj,iBatch) = LnMatij(iSym,iBatch)
          LiMatij(iSymj,iSymi,iBatch) = LnMatij(iSym,iBatch)
          LnMatij(iSym,iBatch) = LnMatij(iSym,iBatch)+LnOcc(iSymi,iBatch)*LnOcc(iSymj,iBatch)
        end if
      end do
    end do
  end do
end if

end subroutine ChoMP2_Setup_Index
