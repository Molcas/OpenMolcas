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
      SUBROUTINE FOCKOC(FOCC,F,CMO)
C
C     RASSCF program version IBM-3090: SX section
C
C     PURPOSE: Construct a fock matrix for the occupied orbitals and
C              write it on the job interphase for later use in
C              the RASSCF gradient programs.
C              The fock matrix is in MO basis with the frozen orbitals
C              excluded. It is symmetry blocked in contrast to earlier
C              versions.
C
C     called from FOCK if IFINAL=1.
C
C          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
C

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='FOCKOC  ')
#include "wadr.fh"
#include "WrkSpc.fh"

      DIMENSION FOCC(*),F(*),CMO(*)

C     CALL QENTER('FOCKOC')
C
      ISTFCK=0
      IPQ=0
* Loop over symmetry:
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
C
      IAD15=IADR15(5)
      CALL DDAFILE(JOBIPH,1,FOCC,IPQ,IAD15)
C
c fock matrices added -- R L 921008.
*
*-----Start processing the all occupied indices Fock matrix
*
      Call GetMem('Scr1','ALLO','REAL',ipScr1,no2m)
      Call GetMem('Scr2','ALLO','REAL',ipScr2,no2m)
* NOTE: FOCC is nowadays allocated already in RASSCF.
* Disk address ipFocc is in /WADR/
*      Call GetMem('FOcc','ALLO','REAL',ipFocc,nTot1)
      call dcopy_(nTot1,[0.0D0],0,Work(ipFocc),1)
#ifdef _DEBUG_
      Write(LF,*) 'nTot1=',nTot1
#endif
      nFock = 0
      Do iSym = 1, nSym
         nFock = nFock + (nISh(iSym)+nASh(iSym))**2
      End Do
      Call GetMem('Fockoc','ALLO','REAL',ipFock,nFock)
      iAd15 = iAdr15(5)
*-----Read Fock matrix in MO basis from JOBIPH.
      Call DDaFile(Jobiph,2,Work(ipFock),nFock,iAd15)
#ifdef _DEBUG_
      Call RecPrt('Fock(MO)',' ',Work(ipFock),1,nFock)
#endif
*
*-----Construct the occupied part of the Fock matrix in SO/AO basis.
*
      ISTFCK= 1! ipFock
      jFock = ipFocc
      iCMo  = 1
* A long loop over symmetry:
      Do ISYM=1,NSYM
*        Hmm maybe you should check so I use correct nbas/norb
         NO = nOrb(isym)
         IF (NO.NE.0) THEN
*
*-----------Transform to SO/AO basis.
*
#ifdef _DEBUG_
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
     &                  0.0d0,Work(ipScr1),nBas(iSym))
            Call DGEMM_('N','T',
     &                  nBas(iSym),nBas(iSym),nOrb(isym),
     &                  1.0d0,Work(ipScr1),nBas(iSym),
     &                  CMO(iCMo),nBas(iSym),
     &                  0.0d0,Work(ipScr2),nBas(iSym))
            ij = jFock
            Do iBas = 1, nBas(iSym)
               Do jBas = 1, iBas-1
                  kl = ipScr2-1 + nBas(iSym)*(jBas-1) + iBas
                  lk = ipScr2-1 + nBas(iSym)*(iBas-1) + jBas
                  Work(ij) = Work(kl) + Work(lk)
                  ij = ij + 1
               END DO
               kl = ipScr2-1 + nBas(iSym)*(iBas-1) + iBas
               If (ij-jFock+1.gt.nTot1) Then
                  Write(LF,*) ij,jFock,nTot1
                  Call Abend()
               End If
               Work(ij) = Work(kl)
               ij = ij + 1
            END DO
         END IF
#ifdef _DEBUG_
         Call TriPrt('FAO',' ',Work(jFock),nBas(iSym))
#endif
         jFock = jFock + nBas(iSym)*(nBas(iSym)+1)/2
         iCMo  = iCMo  + nBas(iSym)**2
         ISTFCK=ISTFCK+NO**2
* End of long loop over symmetry
      END DO
      Call GetMem('Fockoc','FREE','REAL',ipFock,nFock)
      Call GetMem('Scr2','FREE','REAL',ipScr2,no2m)
      Call GetMem('Scr1','FREE','REAL',ipScr1,no2m)
*
c End of addition, R L 921008.
      RETURN
      END
