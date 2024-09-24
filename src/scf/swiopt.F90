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
use InfSO, only: DltNth
use InfSCF, only: DThr, EThr, FThr, nBO, nBT, nIterP, PotNuc
use Files, only: FnDel, FnDGd, FnDSt, FnGrd, FnOSt, FnTSt, Fnx, Fny, LuDel, LuDGd, LuDSt, LuGrd, LuOSt, LuTSt, Lux, Luy
use NDDO, only: twoel_NDDO
use Constants, only: One

implicit none
! declaration of subroutine parameter...
logical AllCnt
integer mBT, mBB, nD
real*8 OneHam(mBT), Ovrlp(mBT), CMO(mBB,nD)
! declaration of some local variables...
real*8 EThr_o, FThr_o, DThr_o, DNTh_o, ThrInt_o
save EThr_o, FThr_o, DThr_o, DNTh_o, ThrInt_o
character(len=8) Label
integer iComp, iD, iOpt, iRC, lOper
real*8, external :: Get_ThrInt

if (AllCnt .and. twoel_NDDO) then
  nIterP = 1
  ! read full overlap matrix from ONEINT file
  Label = 'Mltpl  0'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  iComp = 1
  call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
  if (iRC /= 0) goto 9999
  ! read full one-electron Hamiltonian from ONEINT file
  Label = 'OneHam  '
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
  if (iRc /= 0) goto 9999
  !call Get_PotNuc(PotNuc)
  !call Get_dScalar('PotNuc',PotNuc)
  call Peek_dScalar('PotNuc',PotNuc)
  ! orthonormalize CMO
  do iD=1,nD
    call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
  end do
  ! restore threshold values in WfCtl
  EThr = EThr_o
  call Put_dScalar('EThr',EThr)
  FThr = FThr_o
  DThr = DThr_o
  DltNTh = DNTh_o
  call xSet_ThrInt(ThrInt_o)
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
  if (iRC /= 0) goto 9999
  ! read NDDO NA matrix from ONEINT file...
  Label = 'AttractS'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
  if (iRC /= 0) goto 9999
  ! and form NDDO one-electron Hamiltonian...
  call DaXpY_(nBT,One,Ovrlp,1,OneHam,1)
  ! read NDDO overlap matrix from ONEINT file...
  Label = 'MltplS 0'
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iRC = -1
  call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
  if (iRC /= 0) goto 9999
  ! save threshold values in WfCtl
  EThr_o = EThr
  FThr_o = FThr
  DThr_o = DThr
  DNTh_o = DltNTh
  ThrInt_o = Get_ThrInt()
  ! and set new thresholds
  !EThr = EThr*1.0D+04
  !call Put_dScalar('EThr',EThr)
  !FThr = FThr*1.0D+04
  !DThr = DThr*1.0D+04
  !DltNTh = DltNTh*1.0D+04
  !call xSet_ThrInt(ThrInt_o*1.0D+04)
  ! set twoel to OneCnt...
  twoel_NDDO = .true.
end if

return

! Error exit
9999 continue
write(6,*) 'SwiOpt: Error reading ONEINT'
write(6,'(A,A)') 'Label=',Label
call Abend()

end subroutine SwiOpt
