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

      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*),SPD12(*)
C     REAL*8 COEFF,OVERLAP,OVLP,PRTHR
      REAL*8 COEFF,OVERLAP_RASSI,OVLP
      INTEGER IORBTAB(*),NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER FSBOP,IJSORB,IMODE,ISORB
      INTEGER NDETS1,NDETS2
      INTEGER LFSBANN1,LFSBANN2
      INTEGER JSORB,LANN1,LANN2
      INTEGER ISYOP,MS2OP,KOINFO
      INTEGER ISMLAB,ISPLAB,JSYM,JMS2,JSMLAB,JSPLAB
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      EXTERNAL OVERLAP_RASSI

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
       LFSBANN1=FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1)
       NDETS1=IWORK(LFSBANN1+4)
       CALL GETMEM('ANN1','Allo','Real',LANN1,NDETS1)
       CALL DCOPY_(NDETS1,0.0D0,0,WORK(LANN1),1)
       CALL PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IWORK(LFSBANN1),
     &                   IFSBTAB1,COEFF,WORK(LANN1),PSI1)
CTEST       write(*,*)' Prior to call to PRIMSGM.'
CTEST       write(*,*)' FS block structure at LFSBANN1:'
CTEST       CALL PRFSBTAB(IWORK(LFSBANN1))
CTEST       write(*,*)' Wave function ANN1 after PRIMSGM:'
CTEST       PRTHR=0.01D0
CTEST       CALL PRWVF(IORBTAB,ISSTAB,IWORK(LFSBANN1),PRTHR,WORK(LANN1))
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
        LFSBANN2=FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2)
        NDETS2=IWORK(LFSBANN2+4)
        CALL GETMEM('ANN2','Allo','Real',LANN2,NDETS2)
        CALL DCOPY_(NDETS2,0.0D0,0,WORK(LANN2),1)
        CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                   IFSBTAB2,COEFF,WORK(LANN2),PSI2)

C Compute the spin transition density matrix element:
        OVLP=OVERLAP_RASSI(IWORK(LFSBANN1),
     &                  IWORK(LFSBANN2),WORK(LANN1),WORK(LANN2))
CTEST       write(*,*)' Their overlap:',OVLP
        CALL GETMEM('ANN2','Free','Real',LANN2,NDETS2)
        CALL KILLOBJ(LFSBANN2)
 100    CONTINUE
        IJSORB=ISORB+NASORB*(JSORB-1)
        SPD12(IJSORB)=OVLP
       END DO
       CALL GETMEM('ANN1','Free','Real',LANN1,NDETS1)
       CALL KILLOBJ(LFSBANN1)
      END DO
      RETURN
      END
