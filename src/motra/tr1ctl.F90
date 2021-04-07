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

#ifdef _HDF5_QCM_
use hdf5_utils
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Ovlp(*), HOne(*), Kine(*), CMO(*)
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
integer(kind=iwp) :: IDISK, ISTLT, ISYM, LWDLT, LWDSQ, LWFLT, LWFMO, LWFSQ, LWKAO, LWKMO, LWOVP, LWTMP
real(kind=wp) :: ECOR
#ifdef _HDF5_QCM_
real(kind=wp), allocatable :: writebuf(:,:)
#endif

! Initialize LUONEMO

call DANAME(LUONEMO,FNONEMO)
IDISK = 0
! Provisionally initialize ECOR to prevent false alarms from
! automatic detection of uninitialized variables.
ECOR = Zero
call WR_MOTRA_Info(LUONEMO,1,iDisk,TCONEMO,64,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,8,BSLBL,LENIN8*mxOrb)

! Write Mo coefficients to disc

TCONEMO(1) = IDISK
call dDAFILE(LUONEMO,1,CMO,NTOT2,IDISK)

! Generate Fock-matrix for inactive orbitals
! and compute the total core energy

call GETMEM('FLT','ALLO','REAL',LWFLT,NTOT1)
call GETMEM('DLT','ALLO','REAL',LWDLT,NTOT1)
call GETMEM('FSQ','ALLO','REAL',LWFSQ,NTOT2)
call GETMEM('DSQ','ALLO','REAL',LWDSQ,NTOT2)
call DCOPY_(NTOT1,HONE,1,WORK(LWFLT),1)
call DCOPY_(NTOT2,[Zero],0,WORK(LWFSQ),1)
call DCOPY_(NTOT1,[Zero],0,WORK(LWDLT),1)
call DCOPY_(NTOT2,[Zero],0,WORK(LWDSQ),1)
ECOR = Zero
call FCIN(WORK(LWFLT),NTOT1,WORK(LWDLT),WORK(LWFSQ),WORK(LWDSQ),ECOR,CMO)
call GETMEM('DSQ','FREE','REAL',LWDSQ,NTOT2)
call GETMEM('FSQ','FREE','REAL',LWFSQ,NTOT2)
call GETMEM('DLT','FREE','REAL',LWDLT,NTOT1)

ECOR = POTNUC+ECOR
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A,E20.10)') 'TOTAL CORE ENERGY:',ECOR
end if

! Transform one-electron Fock matrix

call GETMEM('FMO','ALLO','REAL',LWFMO,NORBTT)
call GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
call DCOPY_(NORBTT,[Zero],0,WORK(LWFMO),1)
call DCOPY_(2*N2MAX,[Zero],0,WORK(LWTMP),1)
call TRAONE_MOTRA(WORK(LWFLT),WORK(LWFMO),WORK(LWTMP),CMO)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Fock matrix in MO basis'
  ISTLT = 0
  do ISYM=1,NSYM
    if (NORB(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' symmetry species:',ISYM
      call TRIPRT(' ',' ',WORK(LWFMO+ISTLT),NORB(ISYM))
      ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
    end if
  end do
end if

#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  msym = nsym
  !> put data to file
  datadim(1) = 1; datadim_bound = 1
  call hdf5_put_data(file_id(1),"ecore ",datadim,ecor)
  call hdf5_put_data(file_id(1),"norbtt",datadim,norbtt)
  call hdf5_put_data(file_id(1),"nsym  ",datadim,msym)
  datadim(1) = nsym
  allocate(writebuf(nsym,3)); writebuf = -1
  do i=1,nsym
    writebuf(i,1) = norb(i)
    writebuf(i,2) = nfro(i)
    writebuf(i,3) = ndel(i)
  end do
  call hdf5_put_data(file_id(1),"norb  ",datadim,writebuf(1,1))
  call hdf5_put_data(file_id(1),"nfro  ",datadim,writebuf(1,2))
  call hdf5_put_data(file_id(1),"ndel  ",datadim,writebuf(1,3))
  deallocate(writebuf)
  datadim(1) = norbtt
  call hdf5_put_data(file_id(1),"FockMO",datadim,work(lwfmo:lwfmo+norbtt-1))
end if
#endif

TCONEMO(2) = IDISK
call dDAFILE(LUONEMO,1,WORK(LWFMO),NORBTT,IDISK)
call GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
call GETMEM('FMO','FREE','REAL',LWFMO,NORBTT)
call GETMEM('FLT','FREE','REAL',LWFLT,NTOT1)

! Transform kinetic energy matrix

call GETMEM('KAO','ALLO','REAL',LWKAO,NTOT1)
call GETMEM('KMO','ALLO','REAL',LWKMO,NORBTT)
call GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
call DCOPY_(NORBTT,[Zero],0,WORK(LWKMO),1)
call DCOPY_(2*N2MAX,[Zero],0,WORK(LWTMP),1)
call DCOPY_(NTOT1,KINE,1,WORK(LWKAO),1)
call TRAONE_MOTRA(WORK(LWKAO),WORK(LWKMO),WORK(LWTMP),CMO)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Kinetic integrals in MO basis'
  ISTLT = 0
  do ISYM=1,NSYM
    if (NORB(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' symmetry species:',ISYM
      call TRIPRT(' ',' ',WORK(LWKMO+ISTLT),NORB(ISYM))
      ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
    end if
  end do
end if
TCONEMO(3) = IDISK
call dDAFILE(LUONEMO,1,WORK(LWKMO),NORBTT,IDISK)
call GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
call GETMEM('KMO','FREE','REAL',LWKMO,NORBTT)
call GETMEM('KAO','FREE','REAL',LWKAO,NTOT1)

! Copy overlap matrix to luonem

call GETMEM('OVP','ALLO','REAL',LWOVP,NTOT1)
call DCOPY_(NTOT1,OVLP,1,WORK(LWOVP),1)
TCONEMO(4) = IDISK
call dDAFILE(LUONEMO,1,WORK(LWOVP),NORBTT,IDISK)
call GETMEM('OVP','FREE','REAL',LWOVP,NTOT1)

! Create CIDATA and molecular orbital sections on LUONEM

TCONEMO(5) = IDISK
IDISK = 0
call WR_MOTRA_Info(LUONEMO,1,iDisk,TCONEMO,64,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,8,BSLBL,LENIN8*mxOrb)
call DACLOS(LUONEMO)

return

end subroutine TR1CTL
