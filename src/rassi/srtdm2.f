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
      use stdalloc, only: mma_allocate, mma_deallocate
      use rassi_global_arrays, only: FSBANN1, FSBANN2, FSBANN3
      IMPLICIT NONE
      INTEGER IORBTAB(*), ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      REAL*8 PSI1(*),PSI2(*)
      LOGICAL IF21,IF12
      REAL*8 SRT2M(*)

      REAL*8 COEFF,OVLP
      INTEGER NASORB
      INTEGER IMODE
      INTEGER ISORB,JSORB,LSORB,JLSORB,IJL
      INTEGER ND1,ND2,ND3
      Real*8, EXTERNAL :: OVERLAP_RASSI
      Real*8, Allocatable:: ANN1(:), ANN2(:), ANN3(:)

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
        Call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
        ND1=FSBANN1(5)
        COEFF=1.0D0
        CALL mma_allocate(ANN1,ND1,Label='ANN1')
        ANN1(:)=0.0D0
        CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,
     &                   IFSBTAB1,COEFF,ANN1,PSI1)
CTEST       WRITE(*,*)' The ANN1 wave function, with ISORB=',ISORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,FSBANN1,PRTHR,ANN1)
        DO JSORB=1,NASORB
C Annihilate a single spin orbital from PSI2, the spin orbital JSORB:
         IMODE=-1
         Call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
         ND2=FSBANN2(5)
         COEFF=1.0D0
         CALL mma_allocate(ANN2,ND2,Label='ANN2')
         ANN2(:)=0.0D0
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,
     &                IFSBTAB2,COEFF,ANN2,PSI2)
CTEST       WRITE(*,*)' The ANN2 wave function, with JSORB=',JSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,FSBANN2,PRTHR,ANN2)
         DO LSORB=1,NASORB
C Pair index J,L:
          JLSORB=(NASORB*(JSORB-1))+LSORB-1
          OVLP=0.0D0
C Annihilate another spin orbital from PSI2, LSORB:
          IMODE=-1
          Call FSBOP(IMODE,LSORB,IORBTAB,ISSTAB,FSBANN2,3)
          ND3=FSBANN3(5)
          COEFF=1.0D0
          CALL mma_allocate(ANN3,ND3,Label='ANN3')
          ANN3(:)=0.0D0
          IF (JSORB.ne.LSORB) THEN
           CALL PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,FSBANN3,
     &                  FSBANN2,COEFF,ANN3,ANN2)
C Compute the spin transition density matrix element:
           OVLP=OVERLAP_RASSI(FSBANN1,FSBANN3,ANN1,ANN3)
          ELSE
           OVLP=0.0D0
          END IF
          IJL=ISORB+(NASORB*JLSORB)
          SRT2M(IJL)=OVLP
          CALL mma_deallocate(ANN3)
          CALL mma_deallocate(FSBANN3)
         END DO
         CALL mma_deallocate(ANN2)
         CALL mma_deallocate(FSBANN2)
        END DO
        CALL mma_deallocate(ANN1)
        CALL mma_deallocate(FSBANN1)
       END DO
C ################################################################################
C IF12 = Eliminate to the right (state 2)
      ELSE IF(IF21) THEN
         WRITE(6,*) 'Invalid state combination.
     &         Please, give PSI1=(N-1) and PSI2=(N)'
      ELSE
         WRITE(6,*)'Invalid state combination in 2particle DYSON'
      END IF ! IF10 or IF01
C ################################################################################
      END SUBROUTINE SRTDM2
