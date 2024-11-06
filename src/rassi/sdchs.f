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
* Copyright (C) 2021, Bruno Tenorio                                    *
************************************************************************
      SUBROUTINE SDCHS(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                   PSI1,PSI2,IF20,IF02,SDCHSM)

      use stdalloc, only: mma_allocate, mma_deallocate
      use rassi_global_arrays, only: FSBANN1, FSBANN2
      IMPLICIT NONE

      INTEGER IORBTAB(*), ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      REAL*8 PSI1(*),PSI2(*),SDCHSM(*)
      LOGICAL IF20,IF02

      REAL*8 COEFF,OVLP
      INTEGER NASORB
      INTEGER IMODE
      INTEGER ISORB,JSORB,IJ
      INTEGER ND1,ND2
      Real*8, Allocatable:: ANN1(:), ANN2(:)

      REAL*8, EXTERNAL:: OVERLAP_RASSI

C Calculates the DCH  matrix elements between two states with
C N and N-1 electrons, defined as:
C IF02 D = < N-2 | anni_right anni_right | N >, or
C IF20 D = < N | anni_left anni_left  | N-2 >
C reduced 2-electron tdm in the space of active spin-orbitals
      NASORB=IORBTAB(4)

C IF02 = eliminte one electron to the right: < N-2 | anni_right
C anni_right | N >
      IF(IF02) THEN

       DO ISORB=1,NASORB
C Symmetry properties:
        !ISMLAB=IORBTAB(KOINFO+1+8*(ISORB-1))
        !ISPLAB=IORBTAB(KOINFO+3+8*(ISORB-1))
C Annihilate a single spin orbital from PSI2, the spin orbital ISORB:
        IMODE=-1
        Call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB2,1)
        ND1=FSBANN1(5)
        COEFF=1.0D0
        CALL mma_allocate(ANN1,ND1,Label='ANN1')
        ANN1(:)=0.0D0
        CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,
     &                IFSBTAB2,COEFF,ANN1,PSI2)

        DO JSORB=1,ISORB-1
C Symmetry properties:
         !JSMLAB=IORBTAB(KOINFO+1+8*(JSORB-1))
         !JSPLAB=IORBTAB(KOINFO+3+8*(JSORB-1))
C Pair index J,L:

         OVLP=0.0D0
C Annihilate another spin orbital from PSI2, LSORB:
         IMODE=-1
         Call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN1,2)
         ND2=FSBANN2(5)
         COEFF=1.0D0
         CALL mma_allocate(ANN2,ND2,Label='ANN2')
         ANN2(:)=0.0D0
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,
     &                FSBANN1,COEFF,ANN2,ANN1)

C Compute the spin transition density matrix element:
         OVLP=OVERLAP_RASSI(IFSBTAB1,FSBANN2,PSI1,ANN2)

         IJ=((ISORB-1)*(ISORB-2))/2+JSORB
         SDCHSM(IJ)=SDCHSM(IJ)+OVLP

         CALL mma_deallocate(ANN2)
         CALL mma_deallocate(FSBANN2)
        END DO
        CALL mma_deallocate(ANN1)
        CALL mma_deallocate(FSBANN1)
       END DO

C ################################################################################
C IF02 = Eliminate to the right (state 2)
      ELSE IF(IF20) THEN
       WRITE(6,*) 'Invalid state combination.
     &      Please, give PSI1=(N-2) and PSI2=(N) '
      ELSE
       WRITE(6,*)'Invalid state combination in DCH states'
      END IF

      END SUBROUTINE SDCHS
