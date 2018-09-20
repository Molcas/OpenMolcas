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
      SUBROUTINE TR1CTL_rasscf(Ovlp,HOne,Kine,CMO,DIAF,ITER)
*
*     Objective: Control section for transformation of one-electron
*                integrals (effective one electron Hamiltonian and
*                kinetic energy )
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "fciqmc_global.fh"
#include "trafo_fciqmc.fh"
#include "files_fciqmc.fh"
#include "fciqmc.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
*
      Real*8 Ovlp(*), HOne(*), Kine(*), CMO(*), DIAF(*)
      logical okay, dbg
      integer ITER
*
      Call qEnter('Tr1Ctl_rasscf')
*
*     Initialize LUONEMO
*
cglm      CALL DANAME(LUONEMO,FNONEMO)
      dbg = .false.
      IDISK=0
* Provisionally initialize ECOR to prevent false alarms from
* automatic detection of uninitialized variables.
      ECOR=0.0D0
*
*     Generate Fock-matrix for inactive orbitals
*     and compute the total core energy
*
      CALL GETMEM('FLT','ALLO','REAL',LWFLT,NTOT1)
      CALL GETMEM('DLT','ALLO','REAL',LWDLT,NTOT1)
      CALL GETMEM('FSQ','ALLO','REAL',LWFSQ,NTOT2)
      CALL GETMEM('DSQ','ALLO','REAL',LWDSQ,NTOT2)
      CALL DCOPY_(NTOT1,HONE,1,WORK(LWFLT),1)
      CALL DCOPY_(NTOT2,0.0D0,0,WORK(LWFSQ),1)
      CALL DCOPY_(NTOT1,0.0D0,0,WORK(LWDLT),1)
      CALL DCOPY_(NTOT2,0.0D0,0,WORK(LWDSQ),1)

      WRITE(6,'(A)') 'going to FCIN_rasscf'
      call xflush(6)
      CALL FCIN_rasscf(WORK(LWFLT),NTOT1,WORK(LWDLT),
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
      CALL DCOPY_(NORBTT,0.0D0,0,WORK(LWFMO),1)
      CALL DCOPY_(2*N2MAX,0.0D0,0,WORK(LWTMP),1)
      If (dbg) then
         Write(6,*)
         Write(6,*) ' CMO in tr1ctl_rasscf'
         Write(6,*) ' ---------------------'
         Write(6,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*)(CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
      End If

      If (dbg) then
         WRITE(6,'(6X,A)') 'HONE in tr1ctl_rasscf'
         do iloop=0,NTOT1-1
           WRITE(6,'(2X,D20.8)') (WORK(LWFLT+iloop))
         end do
      end if

      CALL TRAONE_fciqmc(WORK(LWFLT),WORK(LWFMO),WORK(LWTMP),CMO)
c      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,A)') 'Fock matrix in MO basis in tr1ctl_rassc'
        ISTLT=0
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' symmetry species:',ISYM
            do iloop=0,(NORB(isym))*(NORB(isym)+1)/2 -1
             WRITE(6,'(2X,D20.8)') (WORK(LWFMO+ISTLT+iloop))
            end do
c            CALL TRIPRT(' ',' ',WORK(LWFMO+ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
c      END IF
      if(iDoNECI) then
c write one-electron integrals into FCIDMP file
        iorboff = 0
        korboff = 0
        lorboff = 0
        ioff    = 0
        infroff = 0
        do ISYM=1,NSYM
          infroff = infroff + nfro(isym)
          write(6,*) 'iSym, iorboff, norb', iSym, iorboff, norb(iSym)
          IF ( NORB(ISYM).GT.0 ) THEN
            do iorb = iorboff+1,iorboff+norb(isym)
             do jorb = iorboff+1, iorb
c              if(abs(WORK(LWFMO+ioff)).ge.1.0d-11) then
c       write(LuFCI,'(1X,G20.11,4I5)'),WORK(LWFMO+ioff),
c     & iorb,jorb,korboff,lorboff
       write(6,'(1X,G20.11,4I5)') WORK(LWFMO+ioff),
     & iorb,jorb,korboff,lorboff
c              end if
              ioff = ioff + 1
             enddo
            enddo
            iorboff = iorb - 1
          END IF
        end do
c Orbital energies section ....
        IF(ITER.eq.1) then
c read orbital energies from INPORB file
         call f_Inquire (FnInpOrb,okay)
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
c             write(LuFCI,'(1X,G20.11,4I5)'),WORK(ipEOrb+ioff+nfro(isym)+
c     &                                    i-1),i+icount,0,0,0
             write(6,'(1X,G20.11,4I5)') WORK(ipEOrb+ioff+nfro(isym)+
     &                                    i-1),i+icount,0,0,0
c              end if
            enddo
           END IF
           ioff   = ioff + nbas(isym)
           icount = icount + norb(isym)
         end do
         call getmem('EORB','Free','Real',ipEOrb,itotnbas)
        ELSE
c Here orbital energies are updated from the latest Fock Matrix
c We assume that the Fock Matrix is diagonal and that diagonal elements represent
c orbital energies
         ioff   = 0
         icount = 0
         do ISYM=1,NSYM
           IF ( NORB(ISYM).GT.0 ) THEN
            do i = 1,norb(isym)
c             write(LuFCI,'(1X,G20.11,4I5)'),DIAF(ioff+nfro(isym)+
c     &                                   i),i+icount,0,0,0
             write(6,'(1X,G20.11,4I5)') DIAF(ioff+nfro(isym)+
     &                                   i),i+icount,0,0,0
            enddo
           END IF
           ioff   = ioff + nbas(isym)
           icount = icount + norb(isym)
         end do
        END IF
c write core energy into FCIDMP file
c             if(abs(ECOR).ge.1.0d-11) then
c              write(LuFCI,'(1X,G20.11,4I5)'),ECOR, 0, 0, 0, 0
              write(6,'(1X,G20.11,4I5)') ECOR, 0, 0, 0, 0
c             end if
      end if

cglm      TCONEMO(2)=IDISK
cglm      CALL dDAFILE(LUONEMO,1,WORK(LWFMO),NORBTT,IDISK)
      CALL GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
      CALL GETMEM('FMO','FREE','REAL',LWFMO,NORBTT)
      CALL GETMEM('FLT','FREE','REAL',LWFLT,NTOT1)
*
*     Transform kinetic energy matrix
*
cglm      IRC=0
cglm      CALL GETMEM('KAO','ALLO','REAL',LWKAO,NTOT1)
cglm      CALL GETMEM('KMO','ALLO','REAL',LWKMO,NORBTT)
cglm      CALL GETMEM('TMP','ALLO','REAL',LWTMP,2*N2MAX)
cglm      CALL DCOPY_(NORBTT,0.0D0,0,WORK(LWKMO),1)
cglm      CALL DCOPY_(2*N2MAX,0.0D0,0,WORK(LWTMP),1)
cglm      CALL DCOPY_(NTOT1,KINE,1,WORK(LWKAO),1)
cglm      CALL TRAONE_fciqmc(WORK(LWKAO),WORK(LWKMO),WORK(LWTMP),CMO)
cglm      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
cglm        WRITE(6,'(6X,A)') 'Kinetic integrals in MO basis'
cglm        ISTLT=0
cglm        DO ISYM=1,NSYM
cglm          IF ( NORB(ISYM).GT.0 ) THEN
cglm            WRITE(6,'(6X,A,I2)')' symmetry species:',ISYM
cglm            CALL TRIPRT(' ',' ',WORK(LWKMO+ISTLT),NORB(ISYM))
cglm            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
cglm          END IF
cglm        END DO
cglm      END IF
cglm      TCONEMO(3)=IDISK
cglm      CALL dDAFILE(LUONEMO,1,WORK(LWKMO),NORBTT,IDISK)
cglm      CALL GETMEM('TMP','FREE','REAL',LWTMP,2*N2MAX)
cglm      CALL GETMEM('KMO','FREE','REAL',LWKMO,NORBTT)
cglm      CALL GETMEM('KAO','FREE','REAL',LWKAO,NTOT1)
*
*     Copy overlap matrix to luonem
*
cglm      CALL GETMEM('OVP','ALLO','REAL',LWOVP,NTOT1)
cglm      CALL DCOPY_(NTOT1,OVLP,1,WORK(LWOVP),1)
cglm      TCONEMO(4)=IDISK
cglm      CALL dDAFILE(LUONEMO,1,WORK(LWOVP),NORBTT,IDISK)
cglm      CALL GETMEM('OVP','FREE','REAL',LWOVP,NTOT1)
*
*     Create CIDATA and molecular orbital sections on LUONEM
*
cglm      TCONEMO(5)=IDISK
cglm      IDISK=0
cglm      CALL WR_MOTRA_Info(LUONEMO,1,iDisk,
cglm     &                   TCONEMO,64,ECOR,NSYM,
cglm     &                   NBAS,NORB,NFRO,NDEL,8,BSLBL,LENIN8*mxOrb)
cglm      CALL DACLOS(LUONEMO)
*
      Call qExit('Tr1Ctl_rasscf')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_real_array(Ovlp)
         CALL Unused_integer_array(Kine)
      END IF
      END
