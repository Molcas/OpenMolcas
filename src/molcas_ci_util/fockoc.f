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
!#define _DEBUGPRINT_
      SUBROUTINE FOCKOC(FOCC,F,CMO)
!
!     RASSCF program version IBM-3090: SX section
!
!     PURPOSE: Construct a fock matrix for the occupied orbitals and
!              write it on the job interphase for later use in
!              the RASSCF gradient programs.
!              The fock matrix is in MO basis with the frozen orbitals
!              excluded. It is symmetry blocked in contrast to earlier
!              versions.
!
!     called from FOCK if IFINAL=1.
!
!          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
!

      use Definitions, only: LF => u6
      use wadr, only: FockOcc
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

      Real*8 FOCC(*),F(*),CMO(*)
      Real*8, Allocatable:: SCR1(:), SCR2(:)

!
      ISTFCK=0
      IPQ=0
! Loop over symmetry:
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NOO=NIO+NAO
       NO=NORB(ISYM)
       IF(NOO.NE.0) THEN
        DO NP=1,NOO
         DO NQ=1,NOO
          IPQ=IPQ+1
          FOCC(IPQ)=F(ISTFCK+NO*(NQ-1)+NP)
         END DO
        END DO
       ENDIF
       ISTFCK=ISTFCK+NO**2
      END DO
!
      IAD15=IADR15(5)
      CALL DDAFILE(JOBIPH,1,FOCC,IPQ,IAD15)
!
! fock matrices added -- R L 921008.
!
!-----Start processing the all occupied indices Fock matrix
!
      Call mma_allocate(Scr1,no2m,Label='Scr1')
      Call mma_allocate(Scr2,no2m,Label='Scr2')
      FockOcc(:)=0.0D0
#ifdef _DEBUGPRINT_
      Write(LF,*) 'nTot1=',nTot1
#endif
!
!-----Construct the occupied part of the Fock matrix in SO/AO basis.
!
      ISTFCK= 1
      jFock = 1
      iCMo  = 1
! A long loop over symmetry:
      Do ISYM=1,NSYM
!        Hmm maybe you should check so I use correct nbas/norb
         NO = nOrb(isym)
         IF (NO.NE.0) THEN
!
!-----------Transform to SO/AO basis.
!
#ifdef _DEBUGPRINT_
            Write(LF,*) 'iSym=',iSym
            Call RecPrt('F(iStFck)',' ',F(iStFck),nOrb(iSym),
     &                                            nOrb(iSym))
            Call RecPrt('CMO(iCMO)',' ',CMO(iCMO),nBas(iSym),
     &                                            nBas(iSym))
#endif
            Call DGEMM_('N','N',
     &                  nBas(iSym),nOrb(isym),nOrb(isym),
     &                  1.0d0,CMO(iCMo),nBas(iSym),
     &                  F(iStFck),nOrb(isym),
     &                  0.0d0,Scr1,nBas(iSym))
            Call DGEMM_('N','T',
     &                  nBas(iSym),nBas(iSym),nOrb(isym),
     &                  1.0d0,Scr1,nBas(iSym),
     &                  CMO(iCMo),nBas(iSym),
     &                  0.0d0,Scr2,nBas(iSym))
            ij = jFock
            Do iBas = 1, nBas(iSym)
               Do jBas = 1, iBas-1
                  kl = nBas(iSym)*(jBas-1) + iBas
                  lk = nBas(iSym)*(iBas-1) + jBas
                  FockOcc(ij) = Scr2(kl) + Scr2(lk)
                  ij = ij + 1
               END DO
               kl = nBas(iSym)*(iBas-1) + iBas
               If (ij-jFock+1.gt.nTot1) Then
                  Write(LF,*) ij,jFock,nTot1
                  Call Abend()
               End If
               FockOcc(ij) = Scr2(kl)
               ij = ij + 1
            END DO
         END IF
#ifdef _DEBUGPRINT_
         Call TriPrt('FAO',' ',FockOcc(jFock),nBas(iSym))
#endif
         jFock = jFock + nBas(iSym)*(nBas(iSym)+1)/2
         iCMo  = iCMo  + nBas(iSym)**2
         ISTFCK=ISTFCK+NO**2
! End of long loop over symmetry
      END DO
      Call mma_deallocate(Scr2)
      Call mma_deallocate(Scr1)
!
! End of addition, R L 921008.
      END SUBROUTINE FOCKOC
