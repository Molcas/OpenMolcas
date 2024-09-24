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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine Start2(FName,LuOrb,CMO,mBB,nD,Ovrlp,mBT,EOrb,OccNo,mmB)
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals INPORB                            *
!                                                                      *
!     called from: SOrb                                                *
!                                                                      *
!     calls to: Ortho                                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

#ifdef _MSYM_
use, intrinsic :: iso_c_binding, only: c_ptr
#endif
#ifdef _HDF5_
use mh5, only: mh5_exists_dset
#endif
use InfSCF, only: Aufb, FileOrb_id, isHDF5, nBO, nBT, nSym, OnlyProp, VTitle, nOcc, nOrb, nBas, nnB, nDel
use InfSCF, only: mSymON
use Files, only: LuOut
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Two

implicit none
character(len=*) FName
integer LuOrb, mBB, nD, mBT, mmB
real*8 CMO(mBB,nD), Ovrlp(mBT), EOrb(mmB,nD), OccNo(mmB,nD)
integer iBas, iD, iErr, indx, iOff, iOrb, isUHF, iSym, Lu_
integer nTmp(8), iWFtype
character(len=6) OrbName
! Pam 2012 Changed VECSORT arg list, need dummy array:
integer iDummy(1)
integer, dimension(:,:), allocatable :: IndT
#ifdef _MSYM_
type(c_ptr) msym_ctx
#endif
integer iDum(7,8)
real*8 Dummy(1)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Read vectors (vectors MUST !!!!! be written in symmetry blocks
! of dimension nBas(i)*nBas(i))

call mma_allocate(IndT,nnB,nD,Label='IndT')

Lu_ = LuOrb
if (nD == 1) then
  if (isHDF5) then
    call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
  else
    call RdVec_(FName,Lu_,'COEI',nD-1,nSym,nBas,nOrb,CMO,Dummy,OccNo,Dummy,EOrb(1,1),Dummy,IndT(1,1),VTitle,1,iErr,iWFtype)
  end if
  call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT(1,1),0,iDummy,iErr)
  indx = 1
  do iSym=1,nSym
    nTmp(iSym) = 0
    do iBas=1,nBas(iSym)
      if (IndT(indx,1) == 7) nTmp(iSym) = nTmp(iSym)+1
      indx = indx+1
    end do
    if (nOrb(iSym) > nBas(iSym)-nTmp(iSym)) then
      nOrb(iSym) = nBas(iSym)-nTmp(iSym)
      nDel(iSym) = nTmp(iSym)
    end if
  end do
  call TrimCMO(CMO,CMO,nSym,nBas,nOrb)
  call TrimEor(EOrb(1,1),EOrb(1,1),nSym,nBas,nOrb)

  call Setup_SCF()
  if ((.not. Aufb) .and. (.not. OnlyProp)) then
    iOff = 0
    do iSym=1,nSym
      do iOrb=1,nOcc(iSym,1)
        OccNo(iOrb+iOff,1) = Two
      end do
      do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
        OccNo(iOrb+iOff,1) = Zero
      end do
      iOff = iOff+nOrb(iSym)
    end do
  end if
else
  if (isHDF5) then
    isUHF = 0
#   ifdef _HDF5_
    if (mh5_exists_dset(fileorb_id,'MO_ALPHA_VECTORS')) isUHF = 1
#   endif
  else
    call Chk_Vec_UHF(FNAME,Lu_,isUHF)
  end if
  if (isUHF == 1) then
    if (isHDF5) then
      call RdVec_HDF5(fileorb_id,'COEIA',nSym,nBas,CMO(1,1),OccNo(1,1),EOrb(1,1),IndT(1,1))
      call RdVec_HDF5(fileorb_id,'COEIB',nSym,nBas,CMO(1,2),OccNo(1,2),EOrb(1,2),IndT(1,2))
    else
      call RdVec_(FName,Lu_,'COEI',nD-1,nSym,nBas,nOrb,CMO(:,1),CMO(:,2),OccNo(:,1),OccNo(:,2),EOrb(:,1),EOrb(:,2),IndT(:,1), &
                  VTitle,1,iErr,iWFtype)
      call iCopy(nnB,IndT(1,1),1,IndT(1,2),1)
    end if
    call VecSort(nSym,nBas,nBas,CMO(1,1),OccNo(1,1),IndT(1,1),0,iDummy,iErr)
    call VecSort(nSym,nBas,nBas,CMO(1,2),OccNo(1,2),IndT(1,2),0,iDummy,iErr)
    indx = 1
    do iSym=1,nSym
      nTmp(iSym) = 0
      do iBas=1,nBas(iSym)
        if (IndT(indx,1) == 7) nTmp(iSym) = nTmp(iSym)+1
        indx = indx+1
      end do
      if (nOrb(iSym) > nBas(iSym)-nTmp(iSym)) then
        nOrb(iSym) = nBas(iSym)-nTmp(iSym)
        nDel(iSym) = nTmp(iSym)
      end if
    end do
    call TrimCMO(CMO(1,1),CMO(1,1),nSym,nBas,nOrb)
    call TrimEor(EOrb(1,1),EOrb(1,2),nSym,nBas,nOrb)
    call TrimCMO(CMO(1,2),CMO(1,2),nSym,nBas,nOrb)
    call TrimEor(EOrb(1,2),EOrb(1,2),nSym,nBas,nOrb)
    call Setup_SCF()
  else
    if (isHDF5) then
      call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
    else
      call RdVec_(FName,Lu_,'COEI',0,nSym,nBas,nOrb,CMO,Dummy,OccNo,Dummy,EOrb(1,1),Dummy,IndT(1,1),VTitle,1,iErr,iWFtype)
    end if
    call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT(1,1),0,iDummy,iErr)
    indx = 1
    do iSym=1,nSym
      nTmp(iSym) = 0
      do iBas=1,nBas(iSym)
        if (IndT(indx,1) == 7) nTmp(iSym) = nTmp(iSym)+1
        indx = indx+1
      end do
      if (nOrb(iSym) > nBas(iSym)-nTmp(iSym)) then
        nOrb(iSym) = nBas(iSym)-nTmp(iSym)
        nDel(iSym) = nTmp(iSym)
      end if
    end do
    call TrimCMO(CMO,CMO,nSym,nBas,nOrb)
    call TrimEor(EOrb(1,1),EOrb(1,1),nSym,nBas,nOrb)
    call Setup_SCF()
    call dCopy_(nBO,CMO(1,1),1,CMO(1,2),1)
    call dCopy_(nnB,OccNo(1,1),1,OccNo(1,2),1)
    call dCopy_(nnB,EOrb(1,1),1,EOrb(1,2),1)
    call dScal_(nnB,Half,OccNo(1,1),1)
    call dScal_(nnB,Half,OccNo(1,2),1)
  end if
  if (.not. Aufb) then
    iOff = 0
    do iSym=1,nSym
      do iOrb=1,nOcc(iSym,1)
        OccNo(iOrb+iOff,1) = One
      end do
      do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
        OccNo(iOrb+iOff,1) = Zero
      end do
      iOff = iOff+nOrb(iSym)
    end do
    iOff = 0
    do iSym=1,nSym
      do iOrb=1,nOcc(iSym,2)
        OccNo(iOrb+iOff,2) = One
      end do
      do iOrb=nOcc(iSym,2)+1,nOrb(iSym)
        OccNo(iOrb+iOff,2) = Zero
      end do
      iOff = iOff+nOrb(iSym)
    end do
  end if
end if
if (allocated(IndT)) call mma_deallocate(IndT)

if (MSYMON) then
# ifdef _MSYM_
  write(6,*) 'Symmetrizing start orbitals'
  call fmsym_create_context(msym_ctx)
  call fmsym_set_elements(msym_ctx)
  call fmsym_find_symmetry(msym_ctx)
  do iD=1,nD
    call fmsym_symmetrize_orbitals(msym_ctx,CMO(1,iD))
  end do
# else
  write(6,*) 'No msym support, skipping symmetrization of start orbitals...'
# endif
end if

do iD=1,nD
  call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
end do

#ifdef _MSYM_
if (MSYMON) call fmsym_release_context(msym_ctx)
#endif

! Dump orbitals

if (nD == 1) then
  OrbName = 'SCFORB'
  call WrVec_(OrbName,LuOut,'COE',nD-1,nSym,nBas,nBas,CMO,Dummy,OccNo,Dummy,EOrb(1,1),Dummy,iDum,VTitle,iWFtype)
else
  OrbName = 'UHFORB'
  call WrVec_(OrbName,LuOut,'COE',nD-1,nSym,nBas,nBas,CMO(:,1),CMO(:,2),OccNo(:,1),OccNo(:,2),EOrb(:,1),EOrb(:,2),iDum,VTitle, &
              iWFtype)
end if

return

end subroutine Start2
