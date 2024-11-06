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
* Copyright (C) 2018, Jesper Norell                                    *
************************************************************************
      SUBROUTINE MKDYSORB(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                 PSI1,PSI2,IF10,IF01,DYSAMP,DYSCOF)

      use Constants, only: One, Zero
      use stdalloc, only: mma_allocate, mma_deallocate
      use rassi_global_arrays, only: FSBANN1, FSBANN2

      IMPLICIT NONE
      INTEGER IORBTAB(*)
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      REAL*8 PSI1(*),PSI2(*)
      LOGICAL IF10,IF01
      REAL*8 DYSAMP, DYSCOF(*)

      REAL*8 COEFF,OVLP
      Real*8, EXTERNAL :: OVERLAP_RASSI
      INTEGER NASORB
      INTEGER IMODE,ISORB
      INTEGER NDETS1,NDETS2
      INTEGER JSORB
      Real*8, Allocatable:: ANN1(:), ANN2(:)

! +++ J. Norell 12/7 - 2018
C Calculates the Dyson orbital between two states with
C N and N-1 electrons, defined as:
C D = < N-1 | anni_right | N >, or
C D = < N | anni_left | N-1 >

C Nr of active spin-orbitals
      NASORB= IORBTAB(4)
      DYSAMP=Zero
      DO ISORB=1,NASORB
       DYSCOF(ISORB)=Zero
      END DO

C IF10 = Eliminate to the left (state 1)
      IF(IF10) THEN

C Loop over all spin orbitals ISORB:
       DO ISORB=1,NASORB
        OVLP=Zero

C Annihilate a single orbital:
        COEFF=One
        IMODE=-1
        Call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
        NDETS1=FSBANN1(5)
        CALL mma_allocate(ANN1,NDETS1,Label='ANN1')
        ANN1(:)=0.0D0
        CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,
     &                   IFSBTAB1,COEFF,ANN1,PSI1)

C Compute the coefficient as the overlap between the N-1 electron w.f.s
        OVLP=OVERLAP_RASSI(FSBANN1,
     &                  IFSBTAB2,ANN1,PSI2)
        Call mma_deallocate(ANN1)
        Call mma_deallocate(FSBANN1)
        DYSCOF(ISORB)=OVLP

C Collect the squared norm of the Dyson orbital
        DYSAMP=DYSAMP+OVLP*OVLP

       END DO ! ISORB LOOP

C IF01 = Eliminate to the right (state 2)
      ELSE IF(IF01) THEN

C Loop over all spin orbitals JSORB:
       DO JSORB=1,NASORB
         OVLP=Zero

C Annihilate a single orbital:
         COEFF=One
         IMODE=-1
         Call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
         NDETS2=FSBANN2(5)
C BRN
         CALL mma_allocate(ANN2,NDETS2,Label='ANN2')
         ANN2(:)=0.0D0
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,
     &                   IFSBTAB2,COEFF,ANN2,PSI2)

C Compute the coefficient as the overlap between the N-1 electron w.f.s
         OVLP=OVERLAP_RASSI(IFSBTAB1,
     &                  FSBANN2,PSI1,ANN2)
         Call mma_deallocate(ANN2)
         Call mma_deallocate(FSBANN2)
         DYSCOF(JSORB)=OVLP

C Collect the squared norm of the Dyson orbital
         DYSAMP=DYSAMP+OVLP*OVLP

       END DO ! JSORB LOOP

      ELSE
       WRITE(6,*)'Invalid state combination in MKDYSORB'
       WRITE(6,*)'(No such Dyson orbital can exist!)'

      END IF ! IF10 or IF01

C The eventual PES amplitude is given by the squared norm,
C but for transformation of the D_ij elements we need to remove the
C square for now
      DYSAMP = SQRT(DYSAMP)

      END SUBROUTINE MKDYSORB
