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
      SUBROUTINE MKDYORB(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                 PSI1,PSI2,SPD12,IF10,IF01,DYSAMP)

      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*),SPD12(*)
      REAL*8 COEFF,OVERLAP_RASSI,OVLP,DYSAMP
      INTEGER IORBTAB(*),NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER FSBOP,IJSORB,IMODE,ISORB
      INTEGER NDETS1,NDETS2
      INTEGER LFSBANN1,LFSBANN2
      INTEGER JSORB,LANN1,LANN2
      INTEGER KOINFO
      INTEGER ISMLAB,ISPLAB
      LOGICAL IF10,IF01
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      EXTERNAL OVERLAP_RASSI

! +++ J. Norell 12/7 - 2018
C Calculates the Dyson orbital between two states with
C N and N-1 electrons, defined as:
C D = < N-1 | anni_right | N >, or
C D = < N | anni_left | N-1 >
C |D|^2 gives the PES intensity as evaluated within the
C sudden approximation, for more precise treatment
C the Dyson orbital must be numerically overlapped
C with a dipole operator and continuum electron.

C Nr of active spin-orbitals
      NASORB= IORBTAB(4)
      DYSAMP=0.0

C IF10 = Eliminate to the left (state 1)
      IF(IF10) THEN
       KOINFO=19

C Loop over all spin orbitals ISORB:
       DO ISORB=1,NASORB
        OVLP=0.0
        ISMLAB=IORBTAB(KOINFO+1+8*(ISORB-1))
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

C Compute the coefficient as the overlap between the N-1 electron w.f.s
        OVLP=OVERLAP_RASSI(IWORK(LFSBANN1),
     &                  IFSBTAB2,WORK(LANN1),PSI2)
        CALL GETMEM('ANN1','Free','Real',LANN1,NDETS1)
        CALL KILLOBJ(LFSBANN1)
C Collect the squared norm of the Dyson orbital
        DYSAMP=DYSAMP+OVLP*OVLP

       END DO ! ISORB LOOP

C IF01 = Eliminate to the right (state 2)
      ELSE IF(IF01) THEN

C Loop over all spin orbitals JSORB:
       DO JSORB=1,NASORB
         OVLP=0.0

C Annihilate a single orbital:
         COEFF=1.0D0
         IMODE=-1
         LFSBANN2=FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2)
         NDETS2=IWORK(LFSBANN2+4)
         CALL GETMEM('ANN2','Allo','Real',LANN2,NDETS2)
         CALL DCOPY_(NDETS2,0.0D0,0,WORK(LANN2),1)
         CALL PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,IWORK(LFSBANN2),
     &                   IFSBTAB2,COEFF,WORK(LANN2),PSI2)

C Compute the coefficient as the overlap between the N-1 electron w.f.s
         OVLP=OVERLAP_RASSI(IFSBTAB1,
     &                  IWORK(LFSBANN2),PSI1,WORK(LANN2))
        CALL GETMEM('ANN2','Free','Real',LANN2,NDETS2)
        CALL KILLOBJ(LFSBANN2)
C Collect the squared norm of the Dyson orbital
        DYSAMP=DYSAMP+OVLP*OVLP

       END DO ! JSORB LOOP

      ELSE
       WRITE(6,*)'Invalid state combination in MKDYORB'
       WRITE(6,*)'(No such Dyson orbital can exist!)'

      END IF ! IF10 or IF01

C The eventual PES amplitude is given by the squared norm,
C but for transformation of the D_ij elements we need to remove the
C square for now
      DYSAMP = SQRT(DYSAMP)

      RETURN

      END
