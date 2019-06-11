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
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE MKTDAB    IBM-3090 RELEASE 89 01 31
*  PURPOSE: CALCULATE TRANSITION DENSITY MATRIX FOR CI EXPANSIONS IN
*  BIORTHONORMAL ORBITAL BASES A AND B.
*****************************************************************
      SUBROUTINE MKTDAB(OVER,GAMMA1,TDMAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TDMAB(NTDMAB)
      DIMENSION GAMMA1(NASHT,NASHT)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "WrkSpc.fh"
      DIMENSION IOFFA(8)
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO 5 I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
5     CONTINUE
C  INITIALIZE TRANSITION DENSITY MATRIX:
      CALL FZERO(TDMAB,NTDMAB)
C CONTRIBUTION FROM INACTIVE ORBITALS:
      IF(LSYM1.EQ.LSYM2) THEN
       IF(OVER.NE.0.0D0) THEN
        IOFFTD=0
        DO 50 ISY=1,NSYM
          II=0
          DO 40 I=1,NISH(ISY)
            II=II+1
            IPOS=IOFFTD+(II-1)*NOSH(ISY)+II
            TDMAB(IPOS)=2.0D0*OVER
40        CONTINUE
          IOFFTD=IOFFTD+NOSH(ISY)**2
50      CONTINUE
       END IF
      END IF
C THEN ADD CONTRIBUTION FROM ACTIVE SPACE.
      ISY12=MUL(LSYM1,LSYM2)
      IOFFTD=0
      DO 120 ISY1=1,NSYM
        NO1=NOSH(ISY1)
        IF(NO1.EQ.0) GOTO 120
        ISY2=MUL(ISY1,ISY12)
        NO2=NOSH(ISY2)
        IF(NO2.EQ.0) GOTO 120
        NA1=NASH(ISY1)
        IF(NA1.EQ.0) GOTO 110
        NA2=NASH(ISY2)
        IF(NA2.EQ.0) GOTO 110
        NI1=NISH(ISY1)
        NI2=NISH(ISY2)
        DO 100 I=1,NA1
          IA=IOFFA(ISY1)+I
          II=NI1+I
          DO 100 J=1,NA2
            JA=IOFFA(ISY2)+J
            JJ=NI2+J
            IPOS=IOFFTD+II+(JJ-1)*NO1
            TDMAB(IPOS)=GAMMA1(IA,JA)
100     CONTINUE
110     IOFFTD=IOFFTD+NO1*NO2
120   CONTINUE
!      print *, 'TDMAB...'
!      print *, TDMAB(1:NTDMAB)
      RETURN
      END
