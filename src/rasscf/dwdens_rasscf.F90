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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************

! Construct a dynamically-weighted density which is used for reaction field

subroutine DWDens_RASSCF(CMO,D1A,RCT_FS,IFINAL)

use Constants, only: Zero, One
use rasscf_global, only: DoDMRG, ITER, NAC, NACPAR, NACPR2, nRoots, IADR15, Ener
use DWSol, only: DWSol_wgt, W_SOLV
use gas_data, only: iDoGAS
use general_data, only: JOBIPH, NACTEL, NCONF
use gugx, only: SGS
use lucia_data, only: PAtmp, Pscr, PTmp, DStmp, Dtmp
use Lucia_Interface, only: Lucia_Util
use sxci, only: IDXSX
#ifdef _DMRG_
use lucia_data, only: RF1, RF2
use rasscf_global, only: TwoRDM_qcm
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: D1A(*), RCT_FS(*)
integer(kind=iwp), intent(in) :: IFINAL
integer(kind=iwp) :: i, iDisk, iOpt, ITERcurr, jDisk
real(kind=wp) :: rdum(1), wgt
real(kind=wp), allocatable :: CIVEC(:), DA_ave(:), DS_ave(:), DX(:)

call mma_allocate(DA_ave,NAC**2,Label='DA_ave')
call mma_allocate(DS_ave,NAC**2,Label='DS_ave')
call mma_allocate(DX,NACPAR,Label='DX')
DA_ave(1:NACPAR) = Zero
DS_ave(1:NACPAR) = Zero

ITERcurr = 1
if (ITER /= 1) ITERcurr = ITER-1
call DWSol_wgt(2,ENER(:,ITERcurr))

if ((iFinal == 0) .or. (iFinal == 1)) then
  jDisk = IADR15(3)
  do i=1,nRoots
    wgt = W_SOLV(i)
    call DDaFile(JOBIPH,2,DX,NACPAR,jDisk)
    call DDaFile(JOBIPH,2,RCT_FS,NACPAR,jDisk)
    call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
    call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
    if (wgt < 1.0e-10_wp) cycle
    DA_ave(1:NACPAR) = DA_ave(1:NACPAR)+wgt*DX(1:NACPAR)
    DS_ave(1:NACPAR) = DS_ave(1:NACPAR)+wgt*RCT_FS(1:NACPAR)
  end do
else if (iFinal == 2) then
  call mma_allocate(CIVEC,NCONF,Label='CIVEC')
  call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
  call mma_allocate(DStmp,NAC**2,Label='DStmp')
  call mma_allocate(Ptmp,NACPR2,Label='Ptmp')

  iDisk = IADR15(4)
  do i=1,nRoots
    wgt = W_SOLV(i)
    if (NACTEL == 0) then
      CIVEC(1) = One
    else
      if (.not. doDMRG) then
        iOpt = 2
        ! load back one CI vector at the time
        call DDafile(JOBIPH,iOpt,CIVEC,nConf,iDisk)
      end if
    end if

    ! compute density matrices

    if (NAC >= 1) then
      if (NACTEL == 0) then
        Dtmp(:) = Zero
        DStmp(:) = Zero
        Ptmp(:) = Zero
      else
        if (doDMRG) then
#         ifdef _DMRG_
          ! copy the DMs from d1rf/d2rf for ipcmroot
          Dtmp(1:NACPAR) = rf1(1:NACPAR)
          if (twordm_qcm) Ptmp(1:NACPR2) = rf2(1:NACPR2)
          DStmp(:) = Zero
#         endif
        else
          call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
          call mma_allocate(Pscr,NACPR2,Label='Pscr')
          call Lucia_Util('Densi',CI_Vector=CIVEC(:))
          if ((SGS%IFRAS > 2) .or. (iDoGAS)) call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
          call mma_deallocate(Pscr)
          call mma_deallocate(PAtmp)
        end if ! doDMRG/doBLOK or CI
      end if
    else
      Dtmp(:) = Zero
      DStmp(:) = Zero
      Ptmp(:) = Zero
    end if
    DA_ave(1:NACPAR) = DA_ave(1:NACPAR)+wgt*Dtmp(1:NACPAR)
    DS_ave(1:NACPAR) = DS_ave(1:NACPAR)+wgt*DStmp(1:NACPAR)
  end do
  call mma_deallocate(DStmp)
  call mma_deallocate(Dtmp)
  call mma_deallocate(Ptmp)
  call mma_deallocate(CIVEC)
end if

! Construct D-ACTIVE AND D-INACTIVE IN AO BASIS

DX(1:NACPAR) = DS_ave(1:NACPAR)
call DBLOCK(DX)
call Get_D1A_RASSCF(CMO,DX,RCT_FS)

DX(1:NACPAR) = DA_ave(1:NACPAR)
call DBLOCK(DX)
call Get_D1A_RASSCF(CMO,DX,D1A)

call mma_deallocate(DA_ave)
call mma_deallocate(DS_ave)
call mma_deallocate(DX)

return

end subroutine DWDens_RASSCF
