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

#ifdef _HDF5_
use mh5, only: mh5_exists_dset
#endif
use InfSCF, only: nSym, nD, SCF_FileOrb, isHDF5, Tot_EL_Charge, iAU_AB, nOcc, nBas, nOrb, vTitle, FileOrb_ID, nSym
#ifdef _DEBUGPRINT_
use InfSCF, only: Tot_Charge, Tot_Nuc_Charge
#endif
use Constants, only: Zero, Half, One, Two
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
logical OccSet
logical FermSet
logical SpinSet
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
character(len=512) FNAME
logical Idem
real*8, dimension(:,:), allocatable :: OccVec, EpsVec
real*8 Dummy(1), qA, qB, Tmp, Tmp1
integer nVec, iSym, LU, isUHF, LU_, I, iDiff, iOff, N, iBas, iDummy(1), iErr, iWFType

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
call Peek_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
nVec = 0
do iSym=1,nSym
  nVec = nVec+nBas(iSym)
end do
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
    call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
  else
    call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),Dummy,EpsVec(1,1),Dummy,iDummy,VTitle,1,iErr,iWFtype)
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
      call RdVec_HDF5(fileorb_id,'OEA',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
      call RdVec_HDF5(fileorb_id,'OEB',nSym,nBas,Dummy,OccVec(1,2),EpsVec(1,2),iDummy)
    else
      call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),OccVec(1,2),EpsVec(1,1),EpsVec(1,2),iDummy,VTitle,1, &
                  iErr,iWFtype)
    end if
  else
    if (isHDF5) then
      call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
    else
      call RdVec_(FNAME,Lu,'OE',0,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),Dummy,EpsVec(1,1),Dummy,iDummy,VTitle,1,iErr,iWFtype)
    end if
    call dCopy_(nVec,OccVec(1,1),1,OccVec(1,2),1)
    call dCopy_(nVec,EpsVec(1,1),1,EpsVec(1,2),1)
    call dScal_(nVec*nD,half,OccVec,1)
  end if
end if
#ifdef _DEBUGPRINT_
if (nD == 1) then
  write(6,*) 'Orbital energies'
  write(6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
  write(6,*) 'Occupation numbers'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
else
  write(6,*) 'Alpha orbital energies'
  write(6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
  write(6,*) 'Alpha occupation numbers'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(6,*) 'Beta orbital energies'
  write(6,'(10f12.6)') (EpsVec(i,2),i=1,nVec)
  write(6,*) 'Beta occupation numbers'
  write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
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
      else
        OccVec(i,1) = Zero
      end if
    end do
  else
    tmp1 = Zero
    do i=1,nVec
      tmp1 = tmp1+OccVec(i,1)
    end do
    OccVec(:,:) = Zero
    tmp1 = dble(nint(tmp1))
    do i=1,(nint(tmp1)+1)/2
      if (tmp1 >= Two) then
        qa = qa+Two
        OccVec(i,1) = Two
        tmp1 = tmp1-Two
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
    tmp1 = Zero
    do i=1,nVec
      tmp1 = tmp1+OccVec(i,1)+OccVec(i,2)
    end do
    OccVec(:,:) = Zero
    tmp1 = dble(nint(tmp1))
    do i=1,(nint(tmp1)+1)/2
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
  write(6,*) 'chklumo: Idempotency'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(6,'(a,f12.6)') 'Tot charge         ',Tot_charge
  write(6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
  write(6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
  write(6,'(a,f12.6)') 'Electron count     ',Two*qa
else
  write(6,*) 'chklumo: Alpha idempotency'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
  write(6,*) 'chklumo: Beta idempotency'
  write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
  write(6,'(a,f12.6)') 'Tot charge         ',Tot_charge
  write(6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
  write(6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
  write(6,'(a,f12.6)') 'Alpha count        ',qa
  write(6,'(a,f12.6)') 'Beta count         ',qb
end if
#endif
!----------------------------------------------------------------------*
! Is it the same charge.                                               *
!----------------------------------------------------------------------*
if (abs(qa+qb+Tot_el_charge) > half) then
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: System has changed charge!'
# endif
  Occset = .false.
  FermSet = .true.
  goto 999
end if
!----------------------------------------------------------------------*
! Is it the same spin.                                                 *
!----------------------------------------------------------------------*
if (SpinSet) then
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: System might have changed spin!'
  write(6,*) '   iAu_ab=',iAu_ab
  write(6,*) '   qa-qb=',qa-qb
# endif
  idiff = iAu_ab-int(qa-qb)
  if (idiff /= 0) then
#   ifdef _DEBUGPRINT_
    write(6,*) '   yes indeed, spin has changed!'
#   endif
    Occset = .false.
    FermSet = .true.
    goto 999
  end if
end if
!----------------------------------------------------------------------*
! Is it idempotent     D^2 = 2 D                                       *
!----------------------------------------------------------------------*
if (nD == 1) then
  Idem = .true.
  do i=1,nVec
    tmp = half*OccVec(i,1)*(One-half*OccVec(i,1))
    if (abs(tmp) > 0.25d0) Idem = .false.
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: Idempotency'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
# endif
else
  Idem = .true.
  do i=1,nVec
    tmp = OccVec(i,1)*(One-OccVec(i,1))
    if (abs(tmp) > 0.25d0) Idem = .false.
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: Alpha idempotency'
  write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
# endif
  do i=1,nVec
    tmp = OccVec(i,2)*(One-OccVec(i,2))
    if (abs(tmp) > 0.25d0) Idem = .false.
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: Beta idempotency'
  write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
# endif
end if
!----------------------------------------------------------------------*
! Was it idempotent?                                                   *
!----------------------------------------------------------------------*
if (Idem) then
# ifdef _DEBUGPRINT_
  write(6,*) 'chklumo: Idempotent'
# endif
  if (nD == 1) then
    iOff = 0
    do iSym=1,nSym
      n = 0
      do iBas=1,nBas(iSym)
        n = n+int(OccVec(iOff+iBas,1))
      end do
      nOcc(iSym,1) = n/2
      iOff = iOff+nBas(iSym)
    end do
#   ifdef _DEBUGPRINT_
    write(6,'(a,8i5)') 'Occupation       ',(nOcc(i,1),i=1,nSym)
#   endif
  else
    iOff = 0
    do iSym=1,nSym
      n = 0
      do iBas=1,nBas(iSym)
        n = n+int(OccVec(iOff+iBas,1))
      end do
      nOcc(iSym,1) = n
      iOff = iOff+nBas(iSym)
    end do
    iOff = 0
    do iSym=1,nSym
      n = 0
      do iBas=1,nBas(iSym)
        n = n+int(OccVec(iOff+iBas,2))
      end do
      nOcc(iSym,2) = n
      iOff = iOff+nBas(iSym)
    end do
#   ifdef _DEBUGPRINT_
    write(6,'(a,8i5)') 'Alpha occupation ',(nOcc(i,1),i=1,nSym)
    write(6,'(a,8i5)') 'Beta occupation  ',(nOcc(i,2),i=1,nSym)
#   endif
  end if
  Occset = .true.
  FermSet = .false.
else
# ifdef _DEBUGPRINT_
  write(6,*) 'Not idempotent'
# endif
  Occset = .false.
  FermSet = .true.
end if
!----------------------------------------------------------------------*
! Deallocate fields                                                    *
!----------------------------------------------------------------------*
999 continue
call mma_deallocate(EpsVec)
call mma_deallocate(OccVec)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine ChkLumo
