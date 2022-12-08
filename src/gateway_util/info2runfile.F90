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
! Copyright (C) 2006, Roland Lindh                                     *
!               2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Info2Runfile()
!***********************************************************************
!                                                                      *
!     Object: dump misc. information to the runfile.                   *
!                                                                      *
!     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
!             October 2006                                             *
!***********************************************************************

use Period, only: AdCell, Cell_l, ispread, lthCell, VCell
use Basis_Info, only: dbsc, nBas, nCnttp
use Center_Info, only: dc
use External_Centers, only: iXPolType, XF, nXF
use Gateway_global, only: Expert, DirInt
use Sizes_of_Seward, only: S
use RICD_Info, only: Do_RI, Cholesky, Cho_OneCenter, LocalDF
use Gateway_Info, only: CoC, CoM, DoFMM
use Symmetry_Info, only: nIrrep, VarR, VarT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
#include "cholesky.fh"
#include "rctfld.fh"
#include "embpcharg.fh"
#include "localdf.fh"
integer(kind=iwp) :: i, iCnt, iCnttp, iFMM, iGO, iLocalDF, iNTC, iNuc, iOption, iter_S, mdc, nData, nNuc, nDel(8)
logical(kind=iwp) :: Found, Pseudo
integer(kind=iwp), allocatable :: ICh(:), IsMM(:), nStab(:), NTC(:)
real(kind=wp), allocatable :: DCh(:), DCh_Eff(:), DCo(:,:)
character(len=LenIn), allocatable :: xLblCnt(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nDel(:) = 0
call Put_iArray('nFro',nDel,nIrrep) ! put to 0
call qpg_iArray('nDel',Found,nData)
if (.not. Found) then
  call Put_iArray('nDel',nDel,nIrrep)
end if
do i=1,nIrrep ! note that in Basis_Info is nBas(0:7)
  nDel(i) = nBas(i-1)-nDel(i)
end do
call Put_iArray('nOrb',nDel,nIrrep) ! nDel is corrupted here!

! VarR and VarT

if (lRF .and. (.not. PCM)) VarT = .true.
Pseudo = .false.
do iCnttp=1,nCnttp
  Pseudo = Pseudo .or. (dbsc(iCnttp)%pChrg .and. dbsc(iCnttp)%Fixed)
end do
if (.not. DoEMPC) then
  if (allocated(XF) .or. Pseudo) then
    VarR = .true.
    VarT = .true.
  end if
end if

! Manipulate the option flag

iOption = 0
if (DirInt) iOption = ibset(iOption,0)
if (Expert) iOption = ibset(iOption,1)
if (lRF) iOption = ibset(iOption,2)
if (lLangevin .or. (iXPolType > 0)) iOption = ibset(iOption,3)
if (PCM) then
  iOption = ibset(iOption,4)
  nPCM_Info = 0
  call Put_iScalar('PCM info length',nPCM_Info)
end if
iOption = ibset(iOption,5)
! 2el-integrals from the Cholesky vectors
if (Cholesky .or. Do_RI) iOption = ibset(iOption,9)
! RI-Option
if (Do_RI) then
  iOption = ibset(iOption,10)
  ! Local or non-local
  if (LocalDF) then
    call Put_dScalar('LDF Accuracy',Thr_Accuracy)
    call Put_iScalar('LDF Constraint',LDF_Constraint)
    iLocalDF = 1
  else
    iLocalDF = 0
  end if
  call Put_iScalar('DF Mode',iLocalDF)
end if
! 1C-CD
if (Cholesky .and. Cho_1Center) iOption = ibset(iOption,12)
Cho_OneCenter = Cho_1Center
call Put_iScalar('System BitSwitch',iOption)

call Put_iScalar('Highest Mltpl',S%nMltpl)
iGO = 0
call Put_iScalar('Grad ready',iGO)
iFMM = 0
if (DoFMM) iFMM = 1
call Put_iScalar('FMM',iFMM)

call qpg_iScalar('Saddle Iter',Found)
if (.not. Found) then
  iter_S = 0
  call Put_iScalar('Saddle Iter',iter_S)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of unique centers (atoms+pseudo)

nNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc,label='DCo')
call mma_allocate(ICh,nNuc,label='ICh')
call mma_allocate(DCh_Eff,nNuc,label='DCh_Eff')
call mma_allocate(xLblCnt,MxAtom,label='xLblCnt')
mdc = 0
iNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
      DCh_Eff(iNuc) = dbsc(iCnttp)%Charge
      ICh(iNuc) = dbsc(iCnttp)%AtmNr
      xLblCnt(iNuc) = dc(mdc)%LblCnt(1:LenIn)
    end do
  else
    mdc = mdc+dbsc(iCnttp)%nCntr
  end if
end do

call Put_iScalar('Unique centers',nNuc)
call Put_dArray('Un_cen Coordinates',DCo,3*nNuc)
call Put_iArray('Un_cen charge',ICh,nNuc)
call Put_dArray('Un_cen effective charge',DCh_Eff,nNuc)
call Put_cArray('Un_cen Names',xLblCnt(1),LenIn*nNuc)

call mma_deallocate(DCh_Eff)
call mma_deallocate(ICh)
call mma_deallocate(DCo)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of unique atoms

nNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc,label='DCo')
call mma_allocate(DCh,nNuc,label='DCh')
call mma_allocate(DCh_Eff,nNuc,label='DCh_Eff')
call mma_allocate(nStab,nNuc,label='nStab')
mdc = 0
iNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
      DCh_Eff(iNuc) = dbsc(iCnttp)%Charge
      DCh(iNuc) = real(dbsc(iCnttp)%AtmNr,kind=wp)
      xLblCnt(iNuc) = dc(mdc)%LblCnt(1:LenIn)
      nStab(iNuc) = dc(mdc)%nStab
    end do
  else
    mdc = mdc+dbsc(iCnttp)%nCntr
  end if
end do

call Put_iScalar('Unique atoms',nNuc)

call mma_allocate(IsMM,size(dbsc),Label='IsMM')
do i=1,size(dbsc)
  IsMM(i) = dbsc(i)%IsMM
end do
call Put_iArray('IsMM',IsMM,size(dbsc))
call mma_deallocate(IsMM)

call Put_dArray('Unique Coordinates',DCo,3*nNuc)
call Put_dArray('Center of Mass',CoM,3)
call Put_dArray('Center of Charge',CoC,3)
call Put_dArray('Nuclear charge',DCh,nNuc)
call Put_dArray('Effective nuclear Charge',DCh_Eff,nNuc)
call Put_cArray('Unique Atom Names',xLblCnt(1),LenIn*nNuc)
call Put_iArray('nStab',nStab,nNuc)
if (allocated(XF)) call Put_iScalar('nXF',nXF)
if (Cell_l) then
  call Put_dArray('Unit Cell Vector',VCell,9)
  call Put_iArray('Spread of Coord.',ispread,3)
  call Put_iScalar('Unit Cell NAtoms',lthCell)
  call Put_iArray('Unit Cell Atoms',AdCell,lthCell)
end if

! Initiate entry to zero.
DCh_Eff(:) = Zero
call Put_dArray('Mulliken Charge',DCh_Eff,nNuc)

call mma_deallocate(xLblCnt)
call mma_deallocate(nStab)
call mma_deallocate(DCh_Eff)
call mma_deallocate(DCh)
call mma_deallocate(DCo)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate a translation array iNuc -> iCnttp (for espf)

call mma_allocate(NTC,nNuc,label='NTC')

iNTC = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iNuc=1,dbsc(iCnttp)%nCntr
      NTC(iNTC+iNuc) = iCnttp
    end do
    iNTC = iNTC+dbsc(iCnttp)%nCntr
  end if
end do

call Put_iArray('Atom -> Basis',NTC,nNuc)
call mma_deallocate(NTC)
!                                                                      *
!***********************************************************************
!                                                                      *
! Coordinate list as above but for all centers with proper basis
! functions.

nNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc,label='DCo')
iNuc = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
    end do
  end if
end do

call Put_iScalar('Bfn atoms',nNuc)
call Put_dArray('Bfn Coordinates',DCo,3*nNuc)

call mma_deallocate(DCo)
!                                                                      *
!***********************************************************************
!                                                                      *
! Coordinate list as above but for only pseudo centers

nNuc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%pChrg) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc,label='DCo')
call mma_allocate(DCh,nNuc,label='DCh')
iNuc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%pChrg) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
      DCh(iNuc) = real(dbsc(iCnttp)%AtmNr,kind=wp)
    end do
  else
  end if
end do

call Put_iScalar('Pseudo atoms',nNuc)
call Put_dArray('Pseudo Coordinates',DCo,3*nNuc)
call Put_dArray('Pseudo charge',DCh,nNuc)

call mma_deallocate(DCh)
call mma_deallocate(DCo)
!                                                                      *
!***********************************************************************
!                                                                      *
call Mk_ChDisp()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Info2Runfile
