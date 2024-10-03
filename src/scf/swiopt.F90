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
! Copyright (C) 1996, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine SwiOpt(AllCnt,OneHam,Ovrlp,mBT,CMO,mBB,nD)
!***********************************************************************
!                                                                      *
!     purpose: Switch from HF AO to MO optimization                    *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use Gateway_Info, only: ThrInt
use InfSCF, only: DltNth, DThr, EThr, FThr, nBO, nBT, nIterP, PotNuc
use SCFFiles, only: FnDel, FnDGd, FnDSt, FnGrd, FnOSt, FnTSt, Fnx, Fny, LuDel, LuDGd, LuDSt, LuGrd, LuOSt, LuTSt, Lux, Luy
use NDDO, only: twoel_NDDO
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: AllCnt
integer(kind=iwp), intent(in) :: mBT, mBB, nD
real(kind=wp), intent(out) :: OneHam(mBT), Ovrlp(mBT)
real(kind=wp), intent(inout) :: CMO(mBB,nD)
integer(kind=iwp) :: iComp, iD, iOpt, iRC, lOper
real(kind=wp) :: DNTh_o = Zero, DThr_o = Zero, EThr_o = Zero, FThr_o = Zero, ThrInt_o = Zero
character(len=8) :: Label

if (AllCnt .and. twoel_NDDO) then
  nIterP = 1
  ! read full overlap matrix from ONEINT file
  Label = 'Mltpl  0'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  iComp = 1
  call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
  call Error_check()
  ! read full one-electron Hamiltonian from ONEINT file
  Label = 'OneHam  '
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
  call Error_check()
  !call Get_PotNuc(PotNuc)
  !call Get_dScalar('PotNuc',PotNuc)
  call Peek_dScalar('PotNuc',PotNuc)
  ! orthonormalize CMO
  do iD=1,nD
    call Ortho(CMO(:,iD),nBO,Ovrlp,nBT)
  end do
  ! restore threshold values in WfCtl
  EThr = EThr_o
  call Put_dScalar('EThr',EThr)
  FThr = FThr_o
  DThr = DThr_o
  DltNTh = DNTh_o
  ThrInt = ThrInt_o
  ! set twoel to AllCnt...
  twoel_NDDO = .false.
  ! close and reopen some DA files...
  call DaClos(LuDSt)
  call DaClos(LuOSt)
  call DaClos(LuTSt)
  call DaClos(LuGrd)
  call DaClos(LuDGd)
  call DaClos(Lux)
  call DaClos(LuDel)
  call DaClos(Luy)
  call DAName(LuDSt,FnDSt)
  call DAName(LuOSt,FnOSt)
  call DAName(LuTSt,FnTSt)
  call DAName(LuGrd,FnGrd)
  call DAName(LuDGd,FnDGd)
  call DAName(Lux,Fnx)
  call DAName(LuDel,FnDel)
  call DAName(Luy,Fny)
else
  nIterP = 0
  ! read kinetic energy matrix, use space for overlap matrix,
  ! since that one is reread afterwards anyway...
  Label = 'Kinetic '
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  iComp = 1
  call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
  call Error_check()
  ! read NDDO NA matrix from ONEINT file...
  Label = 'AttractS'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
  call Error_check()
  ! and form NDDO one-electron Hamiltonian...
  OneHam(:) = OneHam(:)+Ovrlp(:)
  ! read NDDO overlap matrix from ONEINT file...
  Label = 'MltplS 0'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
  call Error_check()
  ! save threshold values in WfCtl
  EThr_o = EThr
  FThr_o = FThr
  DThr_o = DThr
  DNTh_o = DltNTh
  ThrInt_o = ThrInt
  ! and set new thresholds
  !EThr = EThr*1.0e4_wp
  !call Put_dScalar('EThr',EThr)
  !FThr = FThr*1.0e4_wp
  !DThr = DThr*1.0e4_wp
  !DltNTh = DltNTh*1.0e4_wp
  !ThrInt = ThrInt_o*1.0e4_wp
  ! set twoel to OneCnt...
  twoel_NDDO = .true.
end if

return

contains

subroutine Error_check()

  use Definitions, only: u6

  if (iRC /= 0) then
    write(u6,*) 'SwiOpt: Error reading ONEINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if

end subroutine Error_check

end subroutine SwiOpt
