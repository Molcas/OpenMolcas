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
      SUBROUTINE SPINORB(D,CMO,OCC,kroot)
C
C     Purpose: diagonalize the spin density matrix (D) to
C     obtain the eigenvectors (EVEC) and the eigenvalues (EVAL).
C     Then the natural spinorbitals (CMONSO) are computed
C     (only active).
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
      Parameter (ROUTINE='SPINORB ')
      DIMENSION D(*),CMO(*),OCC(*)
C
C     CALL QENTER('SPINORB')
C
C Local print level (if any)
      IPRLEV=IPRLOC(6)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      IPDEN=1
      IPCMO=1
      IPOCC=1
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        NF=NFRO(ISYM)
        NI=NISH(ISYM)
        IF ( NB.NE.0 ) THEN
          NA=NASH(ISYM)
          IF ( NA.NE.0 ) THEN
            CALL GETMEM('SPORB1','ALLO','REAL',LW1,NA*NA)
            CALL GETMEM('SPORB2','ALLO','REAL',LW2,NA*NB)
            CALL DCOPY_(NA*NA,[0.0D0],0,WORK(LW1),1)
            CALL DCOPY_(NA,[1.0D0],0,WORK(LW1),NA+1)
            CALL Jacob(D(IPDEN),WORK(LW1),NA,NA)
            IDIAG=0
            DO I=1,NA
              IDIAG=IDIAG+I
              OCC(IPOCC+NF+NI+I-1)=D(IPDEN+IDIAG-1)
            END DO
            CALL DGEMM_('N','N',
     &                  NB,NA,NA,
     &                  1.0d0,CMO(IPCMO+(NF+NI)*NB),NB,
     &                  WORK(LW1),NA,
     &                  0.0d0,WORK(LW2),NB)
            CALL DCOPY_(NA*NB,WORK(LW2),1,CMO(IPCMO+(NF+NI)*NB),1)
            CALL GETMEM('SPORB2','FREE','REAL',LW2,NA*NB)
            CALL GETMEM('SPORB1','FREE','REAL',LW1,NA*NA)
            IPDEN=IPDEN+NA*(NA+1)/2
          END IF
          IPCMO=IPCMO+NB*NB
          IPOCC=IPOCC+NB
        END IF
      END DO
C
C     CALL QEXIT('SPINORB')
C
      RETURN
c Avoid unused argument warnings
#ifdef _WARNING_WORKAROUND_
      IF (.FALSE.) CALL Unused_integer(kroot)
#endif
      END
