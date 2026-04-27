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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************

subroutine dens2file(array1,array2,array3,adim,lu,adr,iEmpty,iOpt,iGo,iState,jState)
!***********************************************************************
!   This is the generalized disk interface for TDMs.
!   Some of the parameters and variables are explained here.
!
!   AO_Mode: true if TDMs are in the AO basis, otherwise the TDMs are
!            stored in the basis of the active orbitals only (no sym).
!   iEmpty: the three lowest bits are set if the TDMAB, TSDMAB, and
!           WDMAB, are stored on disk, respectively. That is, for
!           example, if iEmpty=5 the code will write only the first
!           and the last matrix. On read of the second matrix the
!           routine will generate a zero matrix.
!   iOpt: values are 1, or 2, for write and read, respectively.
!   iGo: the three lowest bits are set to tell which matrices are
!        requested. For example, iGo=2, means that only spin-densities
!        are requested.
!
!***********************************************************************

use rassi_aux, only: AO_Mode, CMO1, CMO2, DMAB, Job_Index, mTRA, nasht_save
use stdalloc, only: mma_Allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: adim, lu, iEmpty, iOpt, iGo, iState, jState
real(kind=wp), intent(inout) :: array1(adim), array2(adim), array3(adim)
integer(kind=iwp), intent(inout) :: adr
integer(kind=iwp) :: bdim, IRC, J1, J2, JOB1, JOB1_Old = -1, JOB2, JOB2_Old = -1
real(kind=wp), allocatable :: TRA1(:), TRA2(:)

bdim = adim
if ((.not. AO_Mode) .and. (iOpt == 2)) bdim = nasht_save**2+1
if (bdim > adim) then
  write(u6,*) 'Dens2file: bdim > adim'
  call Abend()
end if

!write(u6,*) 'iState,jState=',iState,jState
if (btest(iGo,0)) then
  if (btest(iEmpty,0)) then
    !if (iOpt == 1) write(u6,*) 'array1=',DDot_(bdim-1,Array1,1,Array1,1),Array1(bdim)
    call ddafile(lu,iOpt,array1,bdim,adr)
  else if (iOpt == 2) then
    array1(:) = Zero
  end if
else if (btest(iEmpty,0)) then
  call ddafile(lu,0,array1,bdim,adr)
end if
if (btest(iGo,1)) then
  if (btest(iEmpty,1)) then
    !if (iOpt == 1) write(u6,*) 'array2=',DDot_(bdim-1,Array2,1,Array2,1),Array2(bdim)
    call ddafile(lu,iOpt,array2,bdim,adr)
  else if (iOpt == 2) then
    array2(:) = Zero
  end if
else if (btest(iEmpty,1)) then
  call ddafile(lu,0,array2,bdim,adr)
end if
if (btest(iGo,2)) then
  if (btest(iEmpty,2)) then
    !if (iOpt == 1) write(u6,*) 'array3=',DDot_(bdim-1,Array3,1,Array3,1),Array3(bdim)
    call ddafile(lu,iOpt,array3,bdim,adr)
  else if (iOpt == 2) then
    array3(:) = Zero
  end if
end if

if ((.not. AO_Mode) .and. (iOpt == 2)) then

  ! Expand the TDMs to AO basis.

  JOB1 = JOB_Index(iState)
  JOB2 = JOB_Index(jState)
  J1 = max(JOB1,JOB2)
  J2 = min(JOB1,JOB2)
  if ((J1 /= JOB1_Old) .or. (J2 /= JOB2_Old)) then
    call RDCMO_RASSI(J1,CMO1)
    call RDCMO_RASSI(J2,CMO2)
    JOB1_OLD = J1
    JOB2_OLD = J2
    call mma_Allocate(TRA1,mTra,Label='TRA1')
    call mma_Allocate(TRA2,mTra,Label='TRA2')
    call FINDT(CMO1,CMO2,TRA1,TRA2)
    call mma_deallocate(TRA2)
    call mma_deallocate(TRA1)
  end if

  if (btest(iGo,0) .and. btest(iEmpty,0)) then
    !write(u6,*) 'array1=',DDot_(nasht_save**2,Array1,1,Array1,1),Array1(bdim)
    call MKTDAB(array1(bdim),array1,DMAB,iRC)
    call MKTDZZ(CMO1,CMO2,DMAB,array1,iRC)
  end if

  if (btest(iGo,1) .and. btest(iEmpty,1)) then
    !write(u6,*) 'array2=',DDot_(nasht_save**2,Array2,1,Array2,1),Array2(bdim)
    call MKTDAB(array2(bdim),array2,DMAB,iRC)
    call MKTDZZ(CMO1,CMO2,DMAB,array2,iRC)
  end if

  if (btest(iGo,2) .and. btest(iEmpty,2)) then
    !write(u6,*) 'array3=',DDot_(nasht_save**2,Array3,1,Array3,1),Array3(bdim)
    call MKTDAB(array3(bdim),array3,DMAB,iRC)
    call MKTDZZ(CMO1,CMO2,DMAB,array3,iRC)
  end if

end if

end subroutine dens2file
