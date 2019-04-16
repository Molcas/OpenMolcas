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
      IMPLICIT NONE
      REAL*8 PSI1(*),PSI4(*),SPD2(*)
C     REAL*8 COEFF,OVERLAP,OVLP,PRTHR
      REAL*8 COEFF,OVERLAP_RASSI,OVLP
      INTEGER IORBTAB(*),NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB4(*)
      INTEGER FSBOP,IMODE,ISORB
      INTEGER LFSBANN1,LFSBANN4
      INTEGER ND1,ND2,ND3,ND4
      INTEGER JSORB,LANN1,LANN2
      INTEGER ISYOP,MS2OP,KOINFO
      INTEGER ISMLAB,ISPLAB,JSMLAB,JSPLAB
      INTEGER NASGEM,JISORB,LFSBANN2
      INTEGER LFSBANN3,LANN3,LANN4,LSORB,KSORB
      INTEGER KLSORB,KLSYM,KLMS2,LSMLAB,LSPLAB,KSMLAB
      INTEGER KSPLAB,IJKL
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      EXTERNAL OVERLAP_RASSI

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
       CALL GETMEM('ANN1','Allo','Real',LANN1,ND1)
       CALL DCOPY_(ND1,[0.0D0],0,WORK(LANN1),1)
       CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IWORK(LFSBANN1),
     &                   IFSBTAB1,COEFF,WORK(LANN1),PSI1)
CTEST       WRITE(*,*)' The ANN1 wave function, with ISORB=',ISORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,WORK(LANN1))
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
        CALL GETMEM('ANN2','Allo','Real',LANN2,ND2)
        CALL DCOPY_(ND2,[0.0D0],0,WORK(LANN2),1)
        CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                IWORK(LFSBANN1),COEFF,WORK(LANN2),WORK(LANN1))
CTEST       WRITE(*,*)' The ANN2 wave function, with JSORB=',JSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN2),PRTHR,WORK(LANN2))

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
         CALL GETMEM('ANN4','Allo','Real',LANN4,ND4)
         CALL DCOPY_(ND4,[0.0D0],0,WORK(LANN4),1)
         CALL PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,IWORK(LFSBANN4),
     &                   IFSBTAB4,COEFF,WORK(LANN4),PSI4)
CTEST       WRITE(*,*)' The ANN4 wave function, with LSORB=',LSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN4),PRTHR,WORK(LANN4))
         DO KSORB=1,LSORB-1
          OVLP=0.0D0
C Symmetry properties:
          KSMLAB=IORBTAB(KOINFO+1+8*(KSORB-1))
          KSPLAB=IORBTAB(KOINFO+3+8*(KSORB-1))
          IF(MUL(KSMLAB,LSMLAB).NE.KLSYM) GOTO 100
          IF(KSPLAB+LSPLAB.NE.KLMS2)      GOTO 100
C Pair index:
          KLSORB=((LSORB-1)*(LSORB-2))/2+KSORB
C Annihilate once more, the spin orbital KSORB:
          IMODE=-1
          LFSBANN3=FSBOP(IMODE,KSORB,IORBTAB,ISSTAB,IWORK(LFSBANN4))
          ND3=IWORK(LFSBANN3+4)
          CALL GETMEM('ANN3','Allo','Real',LANN3,ND3)
          CALL DCOPY_(ND3,[0.0D0],0,WORK(LANN3),1)
          CALL PRIMSGM(IMODE,KSORB,IORBTAB,ISSTAB,IWORK(LFSBANN3),
     &                IWORK(LFSBANN4),COEFF,WORK(LANN3),WORK(LANN4))
CTEST       WRITE(*,*)' The ANN3 wave function, with KSORB=',KSORB
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN3),PRTHR,WORK(LANN3))
C Compute the spin transition density matrix element:
          OVLP=OVERLAP_RASSI(IWORK(LFSBANN2),
     &                  IWORK(LFSBANN3),WORK(LANN2),WORK(LANN3))
CTEST       write(*,*)' Their overlap:',OVLP
            IJKL=JISORB+NASGEM*(KLSORB-1)
            SPD2(IJKL)=SPD2(IJKL)+OVLP
            CALL GETMEM('ANN3','Free','Real',LANN3,ND3)
            CALL KILLOBJ(LFSBANN3)
 100      CONTINUE
         END DO
         CALL GETMEM('ANN4','Free','Real',LANN4,ND4)
         CALL KILLOBJ(LFSBANN4)
        END DO
        CALL GETMEM('ANN2','Free','Real',LANN2,ND2)
        CALL KILLOBJ(LFSBANN2)
       END DO
       CALL GETMEM('ANN1','Free','Real',LANN1,ND1)
       CALL KILLOBJ(LFSBANN1)
      END DO
      RETURN
      END
