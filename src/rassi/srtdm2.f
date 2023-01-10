************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Bruno Tenorio                                    *
************************************************************************
      SUBROUTINE SRTDM2(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                   PSI1,PSI2,IF21,IF12,SRT2M)
      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*),SRT2M(*)
      REAL*8 COEFF,OVERLAP_RASSI,OVLP
      INTEGER IORBTAB(*),NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER FSBOP,IMODE
      INTEGER LFSBANN1,LFSBANN2,LFSBANN3
      INTEGER ISORB,JSORB,LSORB,JLSORB,IJL
      INTEGER LANN1,LANN2,LANN3
      INTEGER ND1,ND2,ND3
      LOGICAL IF21,IF12
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      EXTERNAL OVERLAP_RASSI

C Calculates the 2-electron Dyson matrix between two states with
C N and N-1 electrons, defined as:
C IF12 D = < N-1 | anni_left anni_right anni_right | N >, or
C IF21 D = < N | anni_left anni_left anni_right | N-1 >
C reduced 2-electron tdm in the space of active spin-orbitals

      NASORB=IORBTAB(4)

C IF12 = eliminte one electron to the left: < N-1 | anni_left (PSI1)
C and then eliminate two to the left (PSI2) anni_right anni_right | N >
      IF(IF12) THEN
       DO ISORB=1,NASORB
C Annihilate a single spin orbital from PSI1 (N-1), ISORB:
        IMODE=-1
        LFSBANN1=FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1)
        ND1=IWORK(LFSBANN1+4)
        COEFF=1.0D0
        CALL GETMEM('ANN1','Allo','Real',LANN1,ND1)
        CALL DCOPY_(ND1,[0.0D0],0,WORK(LANN1),1)
        CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IWORK(LFSBANN1),
     &                   IFSBTAB1,COEFF,WORK(LANN1),PSI1)
CTEST       WRITE(*,*)' The ANN1 wave function, with ISORB=',ISORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,WORK(LANN1))
        DO JSORB=1,NASORB
C Annihilate a single spin orbital from PSI2, the spin orbital JSORB:
         IMODE=-1
         LFSBANN2=FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2)
         ND2=IWORK(LFSBANN2+4)
         COEFF=1.0D0
         CALL GETMEM('ANN2','Allo','Real',LANN2,ND2)
         CALL DCOPY_(ND2,[0.0D0],0,WORK(LANN2),1)
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                IFSBTAB2,COEFF,WORK(LANN2),PSI2)
CTEST       WRITE(*,*)' The ANN2 wave function, with JSORB=',JSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN2),PRTHR,WORK(LANN2))
         DO LSORB=1,NASORB
C Pair index J,L:
          JLSORB=(NASORB*(JSORB-1))+LSORB-1
          OVLP=0.0D0
C Annihilate another spin orbital from PSI2, LSORB:
          IMODE=-1
          LFSBANN3=FSBOP(IMODE,LSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2))
          ND3=IWORK(LFSBANN3+4)
          COEFF=1.0D0
          CALL GETMEM('ANN3','Allo','Real',LANN3,ND3)
          CALL DCOPY_(ND3,[0.0D0],0,WORK(LANN3),1)
          IF (JSORB.ne.LSORB) THEN
           CALL PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,IWORK(LFSBANN3),
     &                   IWORK(LFSBANN2),COEFF,WORK(LANN3),WORK(LANN2))
C Compute the spin transition density matrix element:
           OVLP=OVERLAP_RASSI(IWORK(LFSBANN1),
     &                  IWORK(LFSBANN3),WORK(LANN1),WORK(LANN3))
          ELSE
           OVLP=0.0D0
          END IF
          IJL=ISORB+(NASORB*JLSORB)
          SRT2M(IJL)=OVLP
          CALL GETMEM('ANN3','Free','Real',LANN3,ND3)
          CALL KILLOBJ(LFSBANN3)
         END DO
         CALL GETMEM('ANN2','Free','Real',LANN2,ND2)
         CALL KILLOBJ(LFSBANN2)
        END DO
        CALL GETMEM('ANN1','Free','Real',LANN1,ND1)
        CALL KILLOBJ(LFSBANN1)
       END DO
C ################################################################################
C IF12 = Eliminate to the right (state 2)
      ELSE IF(IF21) THEN
      WRITE(6,*) 'Invalid state combination.
     &      Please, give PSI1=(N-1) and PSI2=(N)'
      ELSE
      WRITE(6,*)'Invalid state combination in 2particle DYSON'
      END IF ! IF10 or IF01
C ################################################################################
      RETURN
      END
