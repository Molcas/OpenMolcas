************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE TR1CTL(Ovlp,HOne,Kine,CMO)
*
*     Objective: Control section for transformation of one-electron
*                integrals (effective one electron Hamiltonian and
*                kinetic energy )
*
      !> module dependencies
#ifdef _HDF5_QCM_
      use hdf5_utils
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#ifdef _HDF5_QCM_
      real*8,  allocatable :: writebuf(:,:)
#endif
*
      Real*8 Ovlp(*), HOne(*), Kine(*), CMO(*)
*
*
*     Initialize LUONEMO
*
      CALL DANAME(LUONEMO,FNONEMO)
      IDISK=0
* Provisionally initialize ECOR to prevent false alarms from
* automatic detection of uninitialized variables.
      ECOR=0.0D0
      CALL WR_MOTRA_Info(LUONEMO,1,iDisk,
     &                   TCONEMO,64,ECOR,NSYM,
     &                   NBAS,NORB,NFRO,NDEL,8,BSLBL,LENIN8*mxOrb)
*
*     Write Mo coefficients to disc
*
      TCONEMO(1)=IDISK
      CALL dDAFILE(LUONEMO,1,CMO,NTOT2,IDISK)
*
*     Generate Fock-matrix for inactive orbitals
*     and compute the total core energy
*
      CALL GETMEM('FLT','ALLO','REAL',LWFLT,NTOT1)
      CALL GETMEM('DLT','ALLO','REAL',LWDLT,NTOT1)
      CALL GETMEM('FSQ','ALLO','REAL',LWFSQ,NTOT2)
      CALL GETMEM('DSQ','ALLO','REAL',LWDSQ,NTOT2)
      CALL DCOPY_(NTOT1,HONE,1,WORK(LWFLT),1)
      CALL DCOPY_(NTOT2,[0.0D0],0,WORK(LWFSQ),1)
      CALL DCOPY_(NTOT1,[0.0D0],0,WORK(LWDLT),1)
      CALL DCOPY_(NTOT2,[0.0D0],0,WORK(LWDSQ),1)
      ECOR=0.0D0
      CALL FCIN(WORK(LWFLT),NTOT1,WORK(LWDLT),
     &          WORK(LWFSQ),WORK(LWDSQ),ECOR,CMO)
      CALL GETMEM('DSQ','FREE','REAL',LWDSQ,NTOT2)
      CALL GETMEM('FSQ','FREE','REAL',LWFSQ,NTOT2)
      CALL GETMEM('DLT','FREE','REAL',LWDLT,NTOT1)
*
      ECOR=POTNUC+ECOR
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
         WRITE(6,'(6X,A,E20.10)') 'TOTAL CORE ENERGY:',ECOR
      ENDIF
*
*     Transform one-electron Fock matrix
*
      CALL GETMEM('FMO','ALLO','REAL',LWFMO,NORBTT)
      CALL GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
      CALL DCOPY_(NORBTT,[0.0D0],0,WORK(LWFMO),1)
      CALL DCOPY_(2*N2MAX,[0.0D0],0,WORK(LWTMP),1)
      CALL TRAONE_MOTRA(WORK(LWFLT),WORK(LWFMO),WORK(LWTMP),CMO)
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,A)') 'Fock matrix in MO basis'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' symmetry species:',ISYM
            CALL TRIPRT(' ',' ',WORK(LWFMO+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
      END IF

#ifdef _HDF5_QCM_
      if(ihdf5 == 1)then
        msym = nsym
        !> put data to file
        datadim(1)    = 1; datadim_bound = 1
        call hdf5_put_data(file_id(1),"ecore ",datadim,ecor)
        call hdf5_put_data(file_id(1),"norbtt",datadim,norbtt)
        call hdf5_put_data(file_id(1),"nsym  ",datadim,msym)
        datadim(1)   = nsym
        allocate(writebuf(nsym,3)); writebuf = -1
        do i = 1, nsym
          writebuf(i,1) = norb(i)
          writebuf(i,2) = nfro(i)
          writebuf(i,3) = ndel(i)
        end do
        call hdf5_put_data(file_id(1),"norb  ",datadim,writebuf(1,1))
        call hdf5_put_data(file_id(1),"nfro  ",datadim,writebuf(1,2))
        call hdf5_put_data(file_id(1),"ndel  ",datadim,writebuf(1,3))
        deallocate(writebuf)
        datadim(1)   = norbtt
        call hdf5_put_data(file_id(1),"FockMO",datadim,
     &                     work(lwfmo:lwfmo+norbtt-1))
      end if
#endif

      TCONEMO(2)=IDISK
      CALL dDAFILE(LUONEMO,1,WORK(LWFMO),NORBTT,IDISK)
      CALL GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
      CALL GETMEM('FMO','FREE','REAL',LWFMO,NORBTT)
      CALL GETMEM('FLT','FREE','REAL',LWFLT,NTOT1)
*
*     Transform kinetic energy matrix
*
      CALL GETMEM('KAO','ALLO','REAL',LWKAO,NTOT1)
      CALL GETMEM('KMO','ALLO','REAL',LWKMO,NORBTT)
      CALL GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
      CALL DCOPY_(NORBTT,[0.0D0],0,WORK(LWKMO),1)
      CALL DCOPY_(2*N2MAX,[0.0D0],0,WORK(LWTMP),1)
      CALL DCOPY_(NTOT1,KINE,1,WORK(LWKAO),1)
      CALL TRAONE_MOTRA(WORK(LWKAO),WORK(LWKMO),WORK(LWTMP),CMO)
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,A)') 'Kinetic integrals in MO basis'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' symmetry species:',ISYM
            CALL TRIPRT(' ',' ',WORK(LWKMO+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
      END IF
      TCONEMO(3)=IDISK
      CALL dDAFILE(LUONEMO,1,WORK(LWKMO),NORBTT,IDISK)
      CALL GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
      CALL GETMEM('KMO','FREE','REAL',LWKMO,NORBTT)
      CALL GETMEM('KAO','FREE','REAL',LWKAO,NTOT1)
*
*     Copy overlap matrix to luonem
*
      CALL GETMEM('OVP','ALLO','REAL',LWOVP,NTOT1)
      CALL DCOPY_(NTOT1,OVLP,1,WORK(LWOVP),1)
      TCONEMO(4)=IDISK
      CALL dDAFILE(LUONEMO,1,WORK(LWOVP),NORBTT,IDISK)
      CALL GETMEM('OVP','FREE','REAL',LWOVP,NTOT1)
*
*     Create CIDATA and molecular orbital sections on LUONEM
*
      TCONEMO(5)=IDISK
      IDISK=0
      CALL WR_MOTRA_Info(LUONEMO,1,iDisk,
     &                   TCONEMO,64,ECOR,NSYM,
     &                   NBAS,NORB,NFRO,NDEL,8,BSLBL,LENIN8*mxOrb)
      CALL DACLOS(LUONEMO)
*
      RETURN
      END
