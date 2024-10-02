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
      SUBROUTINE SPIND2(ISYOP,MS2OP,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB4,
     &                 PSI1,PSI4,SPD2)
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      INTEGER ISYOP,MS2OP
      INTEGER IORBTAB(*)
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB4(*)
      REAL*8 PSI1(*),PSI4(*),SPD2(*)

      REAL*8 COEFF,OVLP
      INTEGER NASORB
      INTEGER FSBOP,IMODE,ISORB
      INTEGER LFSBANN1,LFSBANN4
      INTEGER ND1,ND2,ND3,ND4
      INTEGER JSORB
      INTEGER KOINFO
      INTEGER ISMLAB,ISPLAB,JSMLAB,JSPLAB
      INTEGER NASGEM,JISORB,LFSBANN2
      INTEGER LFSBANN3,LSORB,KSORB
      INTEGER KLSORB,KLSYM,KLMS2,LSMLAB,LSPLAB,KSMLAB
      INTEGER KSPLAB,IJKL
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      Real*8, EXTERNAL:: OVERLAP_RASSI
      Real*8, Allocatable:: ANN1(:), ANN2(:), ANN3(:), ANN4(:)

C Nr of active spin-orbitals
      NASORB= IORBTAB(4)
C Nr of active spin-orbital pairs:
      NASGEM=(NASORB*(NASORB-1))/2
      KOINFO=19

      DO ISORB=2,NASORB
C Symmetry properties:
       ISMLAB=IORBTAB(KOINFO+1+8*(ISORB-1))
       ISPLAB=IORBTAB(KOINFO+3+8*(ISORB-1))
C Annihilate a single spin orbital, ISORB:
       IMODE=-1
       LFSBANN1=FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1)
       ND1=IWORK(LFSBANN1+4)
       COEFF=1.0D0
       CALL mma_allocate(ANN1,ND1,Label='ANN1')
       ANN1(:)=0.0D0
       CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IWORK(LFSBANN1),
     &                   IFSBTAB1,COEFF,ANN1,PSI1)
CTEST       WRITE(*,*)' The ANN1 wave function, with ISORB=',ISORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,ANN1)
       DO JSORB=1,ISORB-1
C Symmetry properties:
        JSMLAB=IORBTAB(KOINFO+1+8*(JSORB-1))
        JSPLAB=IORBTAB(KOINFO+3+8*(JSORB-1))
C Pair index:
        JISORB=((ISORB-1)*(ISORB-2))/2+JSORB
C Annihilate once more, the spin orbital JSORB:
        IMODE=-1
        LFSBANN2=FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN1))
        ND2=IWORK(LFSBANN2+4)
        CALL mma_allocate(ANN2,ND2,Label='ANN2')
        ANN2(:)=0.0D0
        CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                IWORK(LFSBANN1),COEFF,ANN2,ANN1)
CTEST       WRITE(*,*)' The ANN2 wave function, with JSORB=',JSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN2),PRTHR,ANN2)

        KLSYM=MUL(MUL(ISMLAB,JSMLAB),ISYOP)
        KLMS2=MS2OP+ISPLAB+JSPLAB
        DO LSORB=2,NASORB
C Symmetry properties:
         LSMLAB=IORBTAB(KOINFO+1+8*(LSORB-1))
         LSPLAB=IORBTAB(KOINFO+3+8*(LSORB-1))
C Annihilate a single spin orbital, LSORB:
         IMODE=-1
         LFSBANN4=FSBOP(IMODE,LSORB,IORBTAB,ISSTAB,IFSBTAB4)
         ND4=IWORK(LFSBANN4+4)
         COEFF=1.0D0
         CALL mma_allocate(ANN4,ND4,Label='ANN4')
         ANN4(:)=0.0D0
         CALL PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,IWORK(LFSBANN4),
     &                   IFSBTAB4,COEFF,ANN4,PSI4)
CTEST       WRITE(*,*)' The ANN4 wave function, with LSORB=',LSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN4),PRTHR,ANN4)
         DO KSORB=1,LSORB-1
          OVLP=0.0D0
C Symmetry properties:
          KSMLAB=IORBTAB(KOINFO+1+8*(KSORB-1))
          KSPLAB=IORBTAB(KOINFO+3+8*(KSORB-1))
          IF(MUL(KSMLAB,LSMLAB).NE.KLSYM) cycle
          IF(KSPLAB+LSPLAB.NE.KLMS2)      cycle
C Pair index:
          KLSORB=((LSORB-1)*(LSORB-2))/2+KSORB
C Annihilate once more, the spin orbital KSORB:
          IMODE=-1
          LFSBANN3=FSBOP(IMODE,KSORB,IORBTAB,ISSTAB,IWORK(LFSBANN4))
          ND3=IWORK(LFSBANN3+4)
          CALL mma_allocate(ANN3,ND3,Label='ANN3')
          ANN3(:)=0.0D0
          CALL PRIMSGM(IMODE,KSORB,IORBTAB,ISSTAB,IWORK(LFSBANN3),
     &                IWORK(LFSBANN4),COEFF,ANN3,ANN4)
CTEST       WRITE(*,*)' The ANN3 wave function, with KSORB=',KSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN3),PRTHR,ANN3)
C Compute the spin transition density matrix element:
          OVLP=OVERLAP_RASSI(IWORK(LFSBANN2),
     &                  IWORK(LFSBANN3),ANN2,ANN3)
CTEST       write(*,*)' Their overlap:',OVLP
            IJKL=JISORB+NASGEM*(KLSORB-1)
            SPD2(IJKL)=SPD2(IJKL)+OVLP
            Call mma_deallocate(ANN3)
            CALL KILLOBJ(LFSBANN3)
         END DO
         Call mma_deallocate(ANN4)
         CALL KILLOBJ(LFSBANN4)
        END DO
        Call mma_deallocate(ANN2)
        CALL KILLOBJ(LFSBANN2)
       END DO
       Call mma_deallocate(ANN1)
       CALL KILLOBJ(LFSBANN1)
      END DO

      END SUBROUTINE SPIND2
