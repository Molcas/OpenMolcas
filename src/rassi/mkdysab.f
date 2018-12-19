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
*               2018, Jesper Norell                                    *
************************************************************************

*****************************************************************
* Modified from MKTAB to MKDYSAB by Jesper Norell, 2018
*  SUBROUTINE MKDYSAB
*  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
*  BIORTHONORMAL ORBITAL BASE A,
*  IN ANALOGUE TO MKTDAB FOR TRANSITION DENSITY MATRIX.
*****************************************************************

      SUBROUTINE MKDYSAB(DYSCOF,DYSAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DYSCOF(*),DYSAB(*)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "WrkSpc.fh"

C *** Symmetry is likely not handled correctly, since the effect
C *** of the annihilated electron is not accounted for.

      DIMENSION IOFFA(8)
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      END DO

C CONTRIBUTION FROM INACTIVE ORBITALS:
C (By definition 0 for Dyson orbitals,
C but we need to fill out the full vector for easier
C transformation.)
      IF(LSYM1.EQ.LSYM2) THEN
        IOFFTD=0
        DO 50 ISY=1,NSYM
          II=0
          DO 40 I=1,NISH(ISY)
            II=II+1
            IPOS=IOFFTD+II
            DYSAB(IPOS)=0.0D0
40        CONTINUE
          IOFFTD=IOFFTD+NOSH(ISY)
50      CONTINUE
      END IF

C THEN ADD CONTRIBUTION FROM ACTIVE SPACE.
      ISY12=MUL(LSYM1,LSYM2)
      IOFFTD=0
      ICOFF=1
      DO 120 ISY1=1,NSYM
        NO1=NOSH(ISY1)
        IF(NO1.EQ.0) GOTO 120
        NA1=NASH(ISY1)
        IF(NA1.EQ.0) GOTO 110
        NI1=NISH(ISY1)
        DO 100 I=1,NA1
          IA=IOFFA(ISY1)+I
          II=NI1+I
          IPOS=IOFFTD+II

! Alpha and beta contributions are added up, in analogue to other rassi
! routines.
! Alpha
          DYSAB(IPOS)=DYSCOF(ICOFF) ! Overwrite "old" values
          ICOFF=ICOFF+1
! Beta
          DYSAB(IPOS)=DYSAB(IPOS)+DYSCOF(ICOFF)
          ICOFF=ICOFF+1
100     CONTINUE
110     IOFFTD=IOFFTD+NO1
120   CONTINUE

      RETURN
      END
