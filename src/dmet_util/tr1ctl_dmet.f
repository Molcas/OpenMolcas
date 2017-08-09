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
      SUBROUTINE TR1CTL_DMET(Ovlp,HOne,Kine,CMO)
*
*     Objective: Control section for transformation of one-electron
*                integrals (effective one electron Hamiltonian and
*                kinetic energy )
*
      !> module dependencies

      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "fciqmc.fh"
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
*
      Real*8 Ovlp(*), HOne(*), Kine(*), CMO(*)
      logical okay
*
      Call qEnter('Tr1Ctl')
*
*     Initialize LUONEMO
*
      FnInpOrb='INPORB'
      FnJobIph='JOBIPH'
      FnOneAO='ONEINT'
      FnTwoAO='ORDINT'
      FnOneMO='TRAONE'
      FnTwoMO='TRAINT'
      FnHalf='TEMP1'
      FnExt='EXTRACT'
      FnCom='COMFILE'
      Debug=0
      iPrint=0
      iOneOnly=0
      iVecTyp=2
      iAutoCut=0
      iRFpert=0
      LuInpOrb=10
      LuJobIph=15
      LuOneAO=20
      LuTwoAO=40
      LuOneMO=30
      LuTwoMO=50
      LuHalf=60
      LuExt=18
      LuCom=22

      write(6,*) "TR1CTL_DMET 0"
      CALL DANAME(LUONEMO,FNONEMO)
      write(6,*) "TR1CTL_DMET 1"
      IDISK=0
* Provisionally initialize ECOR to prevent false alarms from
* automatic detection of uninitialized variables.
      ECOR=0.0D0
      write(6,*) "TR1CTL_DMET 2"
      CALL WR_MOTRA_Info(LUONEMO,1,iDisk,
     &                   TCONEMO,64,ECOR,NSYM,
     &                   NBAS,NORB,NFRO,NDEL,8,BSLBL,2*4*mxOrb)
      write(6,*) "TR1CTL_DMET 3"
*
*     Write Mo coefficients to disc
*
      write(6,*) "TR1CTL_DMET4"
      TCONEMO(1)=IDISK
      CALL dDAFILE(LUONEMO,1,CMO,NTOT2,IDISK)
      write(6,*) "TR1CTL_DMET5"
*
*     Generate Fock-matrix for inactive orbitals
*     and compute the total core energy
*
      CALL GETMEM('FLT','ALLO','REAL',LWFLT,NTOT1)
      write(6,*) "TR1CTL_DMET6"
      CALL GETMEM('DLT','ALLO','REAL',LWDLT,NTOT1)
      write(6,*) "TR1CTL_DMET7"
      CALL GETMEM('FSQ','ALLO','REAL',LWFSQ,NTOT2)
      write(6,*) "TR1CTL_DMET8"
      CALL GETMEM('DSQ','ALLO','REAL',LWDSQ,NTOT2)
      write(6,*) "TR1CTL_DMET9"
      CALL DCOPY_(NTOT1,HONE,1,WORK(LWFLT),1)
      write(6,*) "TR1CTL_DMET10"
      CALL DCOPY_(NTOT2,0.0D0,0,WORK(LWFSQ),1)
      write(6,*) "TR1CTL_DMET11"
      CALL DCOPY_(NTOT1,0.0D0,0,WORK(LWDLT),1)
      write(6,*) "TR1CTL_DMET12"
      CALL DCOPY_(NTOT2,0.0D0,0,WORK(LWDSQ),1)
      write(6,*) "TR1CTL_DMET13"
      ECOR=0.0D0
      CALL FCIN(WORK(LWFLT),NTOT1,WORK(LWDLT),
     &          WORK(LWFSQ),WORK(LWDSQ),ECOR,CMO)
      write(6,*) "TR1CTL_DMET14"
      CALL GETMEM('DSQ','FREE','REAL',LWDSQ,NTOT2)
      write(6,*) "TR1CTL_DMET15"
      CALL GETMEM('FSQ','FREE','REAL',LWFSQ,NTOT2)
      write(6,*) "TR1CTL_DMET16"
      CALL GETMEM('DLT','FREE','REAL',LWDLT,NTOT1)
      write(6,*) "TR1CTL_DMET17"
*
      ECOR=POTNUC+ECOR
*VB      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
         WRITE(6,'(6X,A,E20.10)') 'TOTAL CORE ENERGY:',ECOR
*VB      ENDIF
*
*     Transform one-electron Fock matrix
*
      CALL GETMEM('FMO','ALLO','REAL',LWFMO,NORBTT)
      write(6,*) "TR1CTL_DMET18"
      CALL GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
      write(6,*) "TR1CTL_DMET19"
      CALL DCOPY_(NORBTT,0.0D0,0,WORK(LWFMO),1)
      write(6,*) "TR1CTL_DMET20"
      CALL DCOPY_(2*N2MAX,0.0D0,0,WORK(LWTMP),1)
      write(6,*) "TR1CTL_DMET21"
      CALL TRAONE_MOTRA(WORK(LWFLT),WORK(LWFMO),WORK(LWTMP),CMO)
      write(6,*) "TR1CTL_DMET22"
*VB      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
            write(6,*) "TR1CTL_DMET23"
        WRITE(6,'(6X,A)') 'Fock matrix in MO basis'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' symmetry species:',ISYM
            write(6,*) "TR1CTL_DMET24"
            CALL TRIPRT(' ',' ',WORK(LWFMO+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
*VB      END IF

      if(iDoNECI) then
c write one-electron integrals into FCIDMP file
        iorboff = 0
        korboff = 0
        lorboff = 0
        ioff    = 0
        infroff = 0
        do ISYM=1,NSYM
          infroff = infroff + nfro(isym)
          IF ( NORB(ISYM).GT.0 ) THEN
            do iorb = iorboff+1,iorboff+norb(isym)
             do jorb = iorboff+1, iorb
c              if(abs(WORK(LWFMO+ioff)).ge.1.0d-10)
       write(LuFCI,'(1X,G20.12,4I5)') WORK(LWFMO+ioff),
     & iorb,jorb,korboff,lorboff
              ioff = ioff + 1
             enddo
            enddo
          END IF
          iorboff = iorb - 1
        end do
c read orbital energies from INPORB file
            write(6,*) "TR1CTL_DMET25"
        call f_Inquire (FnInpOrb,okay)
            write(6,*) "TR1CTL_DMET26"
        If ( okay ) Then
          itotnbas = 0
          do i = 1, nsym
            itotnbas= itotnbas + nbas(i)
          end do
          call getmem('EORB','Allo','Real',ipEOrb,itotnbas)
          Call RdVec(FnInpOrb,LuInpOrb,'E',nSym,nBas,nBas,
     &          Dummy, Dummy, work(ipEOrb), iDummy,
     &          VecTit, 0, iErr)
        Else
          Write (6,*) 'RdCMO: Error finding MO file'
          Call QTrace()
          Call Abend()
        End If
c write orbital energies into FCIDUMP file
        ioff   = 0
        icount = 0
        do ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            do i = 1,norb(isym)
             write(LuFCI,'(1X,G20.12,4I5)') WORK(ipEOrb+ioff+nfro(isym)+
     &                                   i-1),i+icount,0,0,0
            enddo
          END IF
          ioff   = ioff + nbas(isym)
          icount = icount + norb(isym)
        end do
        call getmem('EORB','Free','Real',ipEOrb,itotnbas)
c write core energy into FCIDMP file
        write(LuFCI,'(1X,G20.12,4I5)') ECOR, 0, 0, 0, 0
      end if

      TCONEMO(2)=IDISK
      CALL dDAFILE(LUONEMO,1,WORK(LWFMO),NORBTT,IDISK)
      CALL GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
      CALL GETMEM('FMO','FREE','REAL',LWFMO,NORBTT)
      CALL GETMEM('FLT','FREE','REAL',LWFLT,NTOT1)
*
*     Transform kinetic energy matrix
*
      IRC=0
      CALL GETMEM('KAO','ALLO','REAL',LWKAO,NTOT1)
      CALL GETMEM('KMO','ALLO','REAL',LWKMO,NORBTT)
      CALL GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
      CALL DCOPY_(NORBTT,0.0D0,0,WORK(LWKMO),1)
      CALL DCOPY_(2*N2MAX,0.0D0,0,WORK(LWTMP),1)
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
     &                   NBAS,NORB,NFRO,NDEL,8,BSLBL,2*4*mxOrb)
      CALL DACLOS(LUONEMO)
*
      Call qExit('Tr1Ctl')
      RETURN
      END
