!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine TR1CTL(Ovlp,HOne,Kine,CMO)
! Objective: Control section for transformation of one-electron
!            integrals (effective one electron Hamiltonian and
!            kinetic energy )

#include "intent.fh"

#ifdef _HDF5_QCM_
use hdf5_utils, only: datadim, datadim_bound, file_id, hdf5_put_data
use motra_global, only: ihdf5
#endif
use motra_global, only: BsLbl, Debug, FnOneMO, iPrint, LuOneMO, n2max, nBas, nDel, nFro, nOrb, nOrbtt, nSym, nTot1, nTot2, PotNuc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Ovlp(*), HOne(*), Kine(*)
real(kind=wp), intent(_IN_) :: CMO(*)
integer(kind=iwp) :: IDISK, ISTLT, ISYM, TCONEMO(64)
real(kind=wp) :: ECOR
real(kind=wp), allocatable :: DLT(:), DSQ(:), FLT(:), FMO(:), FSQ(:), KAO(:), KMO(:), OVP(:), TMP(:)
#ifdef _HDF5_QCM_
integer(kind=iwp) :: msym
real(kind=wp), allocatable :: writebuf(:,:)
#endif

! Initialize LUONEMO

call DANAME(LUONEMO,FNONEMO)
IDISK = 0
! Provisionally initialize ECOR to prevent false alarms from
! automatic detection of uninitialized variables.
ECOR = Zero
TCONEMO(:) = 0
call WR_MOTRA_Info(LUONEMO,1,iDisk,TCONEMO,64,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,8,BSLBL,size(BSLBL)*len(BSLBL))

! Write Mo coefficients to disc

TCONEMO(1) = IDISK
call dDAFILE(LUONEMO,1,CMO,NTOT2,IDISK)

! Generate Fock-matrix for inactive orbitals
! and compute the total core energy

call mma_allocate(FLT,NTOT1,label='FLT')
call mma_allocate(DLT,NTOT1,label='DLT')
call mma_allocate(FSQ,NTOT2,label='FSQ')
call mma_allocate(DSQ,NTOT2,label='DSQ')
FLT(:) = HONE(1:NTOT1)
FSQ(:) = Zero
DLT(:) = Zero
DSQ(:) = Zero
ECOR = Zero
call FCIN(FLT,NTOT1,DLT,FSQ,DSQ,ECOR,CMO)
call mma_deallocate(DLT)
call mma_deallocate(DSQ)
call mma_deallocate(FSQ)

ECOR = POTNUC+ECOR
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A,ES20.10)') 'TOTAL CORE ENERGY:',ECOR
end if

! Transform one-electron Fock matrix

call mma_allocate(FMO,NORBTT,label='FMO')
call mma_allocate(TMP,2*N2MAX,label='TMP')
FMO(:) = Zero
TMP(:) = Zero
call TRAONE_MOTRA(FLT,FMO,TMP,CMO)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Fock matrix in MO basis'
  ISTLT = 1
  do ISYM=1,NSYM
    if (NORB(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' symmetry species:',ISYM
      call TRIPRT(' ',' ',FMO(ISTLT),NORB(ISYM))
      ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
    end if
  end do
end if

#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  msym = nsym
  !> put data to file
  datadim(1) = 1
  datadim_bound = 1
  call hdf5_put_data(file_id(1),'ecore ',datadim,ecor)
  call hdf5_put_data(file_id(1),'norbtt',datadim,norbtt)
  call hdf5_put_data(file_id(1),'nsym  ',datadim,msym)
  datadim(1) = nsym
  call mma_allocate(writebuf,nsym,3,label='writebuf')
  writebuf(:,1) = norb(1:nsym)
  writebuf(:,2) = nfro(1:nsym)
  writebuf(:,3) = ndel(1:nsym)
  call hdf5_put_data(file_id(1),'norb  ',datadim,writebuf(1,1))
  call hdf5_put_data(file_id(1),'nfro  ',datadim,writebuf(1,2))
  call hdf5_put_data(file_id(1),'ndel  ',datadim,writebuf(1,3))
  call mma_deallocate(writebuf)
  datadim(1) = norbtt
  call hdf5_put_data(file_id(1),'FockMO',datadim,fmo)
end if
#endif

TCONEMO(2) = IDISK
call dDAFILE(LUONEMO,1,FMO,NORBTT,IDISK)
call mma_deallocate(FMO)
call mma_deallocate(FLT)

! Transform kinetic energy matrix

call mma_allocate(KAO,NTOT1,label='KAO')
call mma_allocate(KMO,NORBTT,label='KMO')
KAO(:) = KINE(1:NTOT1)
KMO(:) = Zero
TMP(:) = Zero
call TRAONE_MOTRA(KAO,KMO,TMP,CMO)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Kinetic integrals in MO basis'
  ISTLT = 1
  do ISYM=1,NSYM
    if (NORB(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' symmetry species:',ISYM
      call TRIPRT(' ',' ',KMO(ISTLT),NORB(ISYM))
      ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
    end if
  end do
end if
TCONEMO(3) = IDISK
call dDAFILE(LUONEMO,1,KMO,NORBTT,IDISK)
call mma_deallocate(KAO)
call mma_deallocate(KMO)
call mma_deallocate(TMP)

! Copy overlap matrix to luonem

call mma_allocate(OVP,NTOT1,label='OVP')
OVP(:) = OVLP(1:NTOT1)
TCONEMO(4) = IDISK
call dDAFILE(LUONEMO,1,OVP,NORBTT,IDISK)
call mma_deallocate(OVP)

! Create CIDATA and molecular orbital sections on LUONEM

TCONEMO(5) = IDISK
IDISK = 0
call WR_MOTRA_Info(LUONEMO,1,iDisk,TCONEMO,64,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,8,BSLBL,size(BSLBL)*len(BSLBL))
call DACLOS(LUONEMO)

call mma_deallocate(BSLBL)

return

end subroutine TR1CTL
