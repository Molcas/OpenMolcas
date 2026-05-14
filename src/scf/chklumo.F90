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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ChkLumo(OccSet,FermSet,SpinSet)
!***********************************************************************
!                                                                      *
! This routine figure out the status of the lumo file, i.e. should it  *
! trigger keywords OCCUpied or FERMi?                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

use InfSCF, only: FileOrb_ID, iAU_AB, isHDF5, nBas, nD, nOcc, NoFerm, nOrb, nSym, nSym, SCF_FileOrb, Tot_EL_Charge, vTitle
#ifdef _HDF5_
use mh5, only: mh5_exists_dset
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use InfSCF, only: Tot_Charge, Tot_Nuc_Charge
use Definitions, only: u6
#endif

implicit none
logical(kind=iwp), intent(inout) :: OccSet, FermSet
logical(kind=iwp), intent(in) :: SpinSet
integer(kind=iwp) :: I, iD, iDiff, iDummy(1), iErr, iOff, isUHF, iSym, iWFType, J, LU, LU_, N, ntmp, nVec
real(kind=wp) :: Dummy(1), qA, qB, Tmp, Tmp1
logical(kind=iwp) :: Idem, Skip
character(len=512) :: FNAME
integer(kind=iwp), allocatable :: Ind(:,:)
real(kind=wp), allocatable :: EpsVec(:,:), OccVec(:,:)

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
call Peek_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
nVec = sum(nBas(1:nSym))
!----------------------------------------------------------------------*
! Allocate fields                                                      *
!----------------------------------------------------------------------*
call mma_Allocate(OccVec,nVec,nD,Label='OccVec')
call mma_Allocate(EpsVec,nVec,nD,Label='EpsVec')
!----------------------------------------------------------------------*
! Read occupation numbers and orbital energies                         *
!----------------------------------------------------------------------*
Lu = 17
FName = SCF_FileOrb
if (nD == 1) then
  if (isHDF5) then
    call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(:,1),EpsVec(:,1),iDummy)
  else
    call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(:,1),Dummy,EpsVec(:,1),Dummy,iDummy,VTitle,1,iErr,iWFtype)
  end if
else
  isUHF = 0
  if (isHDF5) then
#   ifdef _HDF5_
    if (mh5_exists_dset(fileorb_id,'MO_ALPHA_VECTORS')) isUHF = 1
#   endif
  else
    Lu_ = 18
    isUHF = -1
    call Chk_Vec_UHF(FNAME,Lu_,isUHF)
  end if
  if (isUHF == 1) then
    if (isHDF5) then
      call RdVec_HDF5(fileorb_id,'OEA',nSym,nBas,Dummy,OccVec(:,1),EpsVec(:,1),iDummy)
      call RdVec_HDF5(fileorb_id,'OEB',nSym,nBas,Dummy,OccVec(:,2),EpsVec(:,2),iDummy)
    else
      call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(:,1),OccVec(:,2),EpsVec(:,1),EpsVec(:,2),iDummy,VTitle,1, &
                  iErr,iWFtype)
    end if
  else
    if (isHDF5) then
      call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(:,1),EpsVec(:,1),iDummy)
    else
      call RdVec_(FNAME,Lu,'OE',0,nSym,nBas,nOrb,Dummy,Dummy,OccVec(:,1),Dummy,EpsVec(:,1),Dummy,iDummy,VTitle,1,iErr,iWFtype)
    end if
    OccVec(:,2) = OccVec(:,1)
    EpsVec(:,2) = EpsVec(:,1)
    OccVec(:,:) = Half*OccVec(:,:)
  end if
end if
#ifdef _DEBUGPRINT_
if (nD == 1) then
  write(u6,*) 'Orbital energies'
  write(u6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
  write(u6,*) 'Occupation numbers'
  write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
else
  write(u6,*) 'Alpha orbital energies'
  write(u6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
  write(u6,*) 'Alpha occupation numbers'
  write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(u6,*) 'Beta orbital energies'
  write(u6,'(10f12.6)') (EpsVec(i,2),i=1,nVec)
  write(u6,*) 'Beta occupation numbers'
  write(u6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
end if
#endif
!----------------------------------------------------------------------*
! What are the charges                                                 *
!----------------------------------------------------------------------*
qa = Zero
qb = Zero
if (nD == 1) then
  if (nSym /= 1) then
    tmp1 = Zero
    do i=1,nVec
      tmp1 = tmp1+OccVec(i,1)
      if (tmp1 >= Two) then
        qa = qa+Two
        tmp1 = tmp1-Two
        OccVec(i,1) = Two
      else if (tmp1 >= One) then
        qa = qa+One
        tmp1 = tmp1-One
        OccVec(i,1) = One
      else
        OccVec(i,1) = Zero
      end if
    end do
  else
    ntmp = nint(sum(OccVec(:,1)))
    OccVec(:,:) = Zero
    tmp1 = real(ntmp,kind=wp)
    do i=1,(ntmp+1)/2
      if (tmp1 >= Two) then
        qa = qa+Two
        OccVec(i,1) = Two
        tmp1 = tmp1-Two
      else if (tmp1 >= One) then
        qa = qa+One
        OccVec(i,1) = One
        tmp1 = tmp1-One
      end if
    end do
  end if
  qa = Half*qa
  qb = qa
else
  if (nSym /= 1) then
    tmp1 = Zero
    do i=1,nVec
      tmp1 = tmp1+OccVec(i,1)+OccVec(i,2)
      if (tmp1 >= Two) then
        qa = qa+One
        qb = qb+One
        tmp1 = tmp1-Two
        OccVec(i,1) = One
        OccVec(i,2) = One
      else if (tmp1 >= One) then
        qa = qa+One
        qb = qb+Zero
        OccVec(i,1) = One
        OccVec(i,2) = Zero
        tmp1 = tmp1-One
      else
        OccVec(i,1) = Zero
        OccVec(i,2) = Zero
      end if
    end do
  else
    ntmp = nint(sum(OccVec(:,1:2)))
    OccVec(:,:) = Zero
    tmp1 = real(ntmp,kind=wp)
    do i=1,(ntmp+1)/2
      if (tmp1 >= Two) then
        qa = qa+One
        qb = qb+One
        OccVec(i,1) = One
        OccVec(i,2) = One
        tmp1 = tmp1-Two
      else if (tmp1 >= One) then
        qa = qa+One
        qb = qb+Zero
        OccVec(i,1) = One
        OccVec(i,2) = Zero
        tmp1 = tmp1-One
      else
        OccVec(i,1) = Zero
        OccVec(i,2) = Zero
      end if
    end do
  end if
end if
#ifdef _DEBUGPRINT_
if (nD == 1) then
  write(u6,*) 'chklumo: Idempotency'
  write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(u6,'(a,f12.6)') 'Tot charge         ',Tot_charge
  write(u6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
  write(u6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
  write(u6,'(a,f12.6)') 'Electron count     ',Two*qa
else
  write(u6,*) 'chklumo: Alpha idempotency'
  write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(u6,*) 'chklumo: Beta idempotency'
  write(u6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
  write(u6,'(a,f12.6)') 'Tot charge         ',Tot_charge
  write(u6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
  write(u6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
  write(u6,'(a,f12.6)') 'Alpha count        ',qa
  write(u6,'(a,f12.6)') 'Beta count         ',qb
end if
#endif
!----------------------------------------------------------------------*
! Is it the same charge.                                               *
!----------------------------------------------------------------------*
Skip = .false.
if (abs(qa+qb+Tot_el_charge) > half) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'chklumo: System has changed charge!'
# endif
  Occset = .false.
  FermSet = .true.
  Skip = .true.
end if
!----------------------------------------------------------------------*
! Is it the same spin.                                                 *
!----------------------------------------------------------------------*
if (SpinSet .and. (.not. Skip)) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'chklumo: System might have changed spin!'
  write(u6,*) '   iAu_ab=',iAu_ab
  write(u6,*) '   qa-qb=',qa-qb
# endif
  idiff = iAu_ab-int(qa-qb)
  if (idiff /= 0) then
#   ifdef _DEBUGPRINT_
    write(u6,*) '   yes indeed, spin has changed!'
#   endif
    Occset = .false.
    FermSet = .true.
    Skip = .true.
  end if
end if
if (.not. Skip) then
  !----------------------------------------------------------------------*
  ! Is it idempotent     D^2 = 2 D                                       *
  !----------------------------------------------------------------------*
  if (nD == 1) then
    Idem = .true.
    do i=1,nVec
      tmp = half*OccVec(i,1)*(One-half*OccVec(i,1))
      if (abs(tmp) > 0.1_wp) Idem = .false.
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'chklumo: Idempotency'
    write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
#   endif
  else
    Idem = .true.
    do i=1,nVec
      tmp = OccVec(i,1)*(One-OccVec(i,1))
      if (abs(tmp) > 0.1_wp) Idem = .false.
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'chklumo: Alpha idempotency'
    write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
#   endif
    do i=1,nVec
      tmp = OccVec(i,2)*(One-OccVec(i,2))
      if (abs(tmp) > 0.1_wp) Idem = .false.
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'chklumo: Beta idempotency'
    write(u6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
#   endif
  end if
  !----------------------------------------------------------------------*
  ! Was it idempotent?                                                   *
  !----------------------------------------------------------------------*
  if (Idem) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'chklumo: Idempotent'
#   endif
    if (nD == 1) then
      iOff = 0
      do iSym=1,nSym
        n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),1)))
        nOcc(iSym,1) = n/2
        iOff = iOff+nBas(iSym)
      end do
#     ifdef _DEBUGPRINT_
      write(u6,'(a,8i5)') 'Occupation       ',(nOcc(i,1),i=1,nSym)
#     endif
    else
      iOff = 0
      do iSym=1,nSym
        n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),1)))
        nOcc(iSym,1) = n
        iOff = iOff+nBas(iSym)
      end do
      iOff = 0
      do iSym=1,nSym
        n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),2)))
        nOcc(iSym,2) = n
        iOff = iOff+nBas(iSym)
      end do
#     ifdef _DEBUGPRINT_
      write(u6,'(a,8i5)') 'Alpha occupation ',(nOcc(i,1),i=1,nSym)
      write(u6,'(a,8i5)') 'Beta occupation  ',(nOcc(i,2),i=1,nSym)
#     endif
    end if
    Occset = .true.
    FermSet = .false.
  else
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Not idempotent'
#   endif
    Occset = .false.
    FermSet = .true.
  end if
end if
! If Fermi aufbau is explicitly disabled, force the plain occupation
if (FermSet .and. NoFerm) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'Fermi aufbau disabled by the user'
# endif
  call mma_allocate(Ind,nVec,nD,label='Ind')
  Ind(:,:) = -1
  ! Sort orbital energies
  do iD=1,nD
    do i=1,nVec
      Tmp = huge(Tmp)
      ntmp = -1
      do j=1,nVec
        if (EpsVec(j,iD) < Tmp) then
          ntmp = j
          Tmp = EpsVec(j,iD)
        end if
      end do
      if (ntmp < 0) exit
      Ind(i,iD) = ntmp
      EpsVec(ntmp,iD) = huge(Tmp)
    end do
  end do
  OccVec(:,:) = Zero
  if (nD == 1) then
    ntmp = nint(-Tot_el_charge)/2
    do i=1,ntmp
      if (Ind(i,1) < 0) then
        call WarningMessage(2,'chklumo: Cannot find meaningful occupations')
        call Abend()
      end if
      OccVec(Ind(i,1),1) = Two
    end do
    iOff = 0
    do iSym=1,nSym
      n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),1)))
      nOcc(iSym,1) = n/2
      iOff = iOff+nBas(iSym)
    end do
#   ifdef _DEBUGPRINT_
    write(u6,'(a,8i5)') 'Forced occupation       ',(nOcc(i,1),i=1,nSym)
    write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
#   endif
  else
    ntmp = (nint(-Tot_el_charge)-iAU_AB)/2
    do i=1,ntmp+iAU_AB
      if (Ind(i,1) < 0) then
        call WarningMessage(2,'chklumo: Cannot find meaningful occupations')
        call Abend()
      end if
      OccVec(Ind(i,1),1) = One
      if (i > ntmp) cycle
      if (Ind(i,2) < 0) then
        call WarningMessage(2,'chklumo: Cannot find meaningful occupations')
        call Abend()
      end if
      OccVec(Ind(i,2),2) = One
    end do
    iOff = 0
    do iSym=1,nSym
      n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),1)))
      nOcc(iSym,1) = n
      iOff = iOff+nBas(iSym)
    end do
    iOff = 0
    do iSym=1,nSym
      n = sum(int(OccVec(iOff+1:iOff+nBas(iSym),2)))
      nOcc(iSym,2) = n
      iOff = iOff+nBas(iSym)
    end do
#   ifdef _DEBUGPRINT_
    write(u6,'(a,8i5)') 'Forced alpha occupation ',(nOcc(i,1),i=1,nSym)
    write(u6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
    write(u6,'(a,8i5)') 'Forced beta occupation  ',(nOcc(i,2),i=1,nSym)
    write(u6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
#   endif
  end if
  call mma_deallocate(Ind)
  OccSet = .true.
  FermSet = .false.
end if
!----------------------------------------------------------------------*
! Deallocate fields                                                    *
!----------------------------------------------------------------------*
call mma_deallocate(EpsVec)
call mma_deallocate(OccVec)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine ChkLumo
