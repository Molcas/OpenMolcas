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
      SUBROUTINE SPIND(ISYOP,MS2OP,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                 PSI1,PSI2,SPD12)
      use stdalloc, only: mma_allocate, mma_deallocate
      use rassi_global_arrays, only: FSBANN1, FSBANN2
      use Symmetry_Info, only: MUL
      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*),SPD12(*)
      REAL*8 COEFF,OVERLAP_RASSI,OVLP
      INTEGER IORBTAB(*),NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER IJSORB,IMODE,ISORB
      INTEGER NDETS1,NDETS2
      INTEGER JSORB
      INTEGER ISYOP,MS2OP,KOINFO
      INTEGER ISMLAB,ISPLAB,JSYM,JMS2,JSMLAB,JSPLAB
      EXTERNAL OVERLAP_RASSI
      Real*8, Allocatable:: ANN1(:), ANN2(:)

CTEST      write(*,*)' Test prints in SPIND.'
C Nr of active spin-orbitals
      NASORB= IORBTAB(4)
C Loop over pairs of spin orbitals ISORB,JSORB:
      KOINFO=19
      DO ISORB=1,NASORB
       ISMLAB=IORBTAB(KOINFO+1+8*(ISORB-1))
CUNUSED       ISOIND=IORBTAB(KOINFO+2+8*(ISORB-1))
       ISPLAB=IORBTAB(KOINFO+3+8*(ISORB-1))
C Annihilate a single orbital:
       COEFF=1.0D0
       IMODE=-1
       Call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
       NDETS1=FSBANN1(5)
       CALL mma_allocate(ANN1,NDETS1,Label='ANN1')
       ANN1(:)=0.0D0
       CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,
     &                   IFSBTAB1,COEFF,ANN1,PSI1)
CTEST       write(*,*)' Prior to call to PRIMSGM.'
CTEST       write(*,*)' FS block structure at FSBANN1:'
CTEST       CALL PRFSBTAB(FSBANN1)
CTEST       write(*,*)' Wave function ANN1 after PRIMSGM:'
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,FSBANN1,PRTHR,ANN1)
C Compute those ANN2 wave functions, which have the correct properties:
       JSYM=MUL(ISMLAB,ISYOP)
       JMS2= MS2OP+ISPLAB
       DO JSORB=1,NASORB
        OVLP=0.0D0
        JSMLAB=IORBTAB(KOINFO+1+8*(JSORB-1))
        IF(JSMLAB.NE.JSYM) GOTO 100
CUNUSED        JSOIND=IORBTAB(KOINFO+2+8*(JSORB-1))
        JSPLAB=IORBTAB(KOINFO+3+8*(JSORB-1))
        IF(JSPLAB.NE.JMS2) GOTO 100
        COEFF=1.0D0
        IMODE=-1
        Call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
        NDETS2=FSBANN2(5)
        CALL mma_allocate(ANN2,NDETS2,Label='ANN2')
        ANN2(:)=0.0D0
        CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,
     &                   IFSBTAB2,COEFF,ANN2,PSI2)

C Compute the spin transition density matrix element:
        OVLP=OVERLAP_RASSI(FSBANN1,FSBANN2,ANN1,ANN2)
CTEST       write(*,*)' Their overlap:',OVLP
        CALL mma_deallocate(ANN2)
        CALL mma_deallocate(FSBANN2)
 100    CONTINUE
        IJSORB=ISORB+NASORB*(JSORB-1)
        SPD12(IJSORB)=OVLP
       END DO
       CALL mma_deallocate(ANN1)
       CALL mma_deallocate(FSBANN1)
      END DO

      END SUBROUTINE SPIND
