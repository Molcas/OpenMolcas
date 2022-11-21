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
*               2020, Bruno Tenorio                                    *
************************************************************************

**********************************************************************
* Modified from MKTAB to MKDYSAB by Jesper Norell, 2018
*  SUBROUTINE MKDYSAB
*  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
*  BIORTHONORMAL ORBITAL BASE A,
*  IN ANALOGUE TO MKTDAB FOR TRANSITION DENSITY MATRIX.
**********************************************************************
*  MODIFIED BY BRUNO TENORIO TO ADDRESS SYMMETRY
*  SEPTEMBER 2020
**********************************************************************

      SUBROUTINE MKDYSAB(DYSCOF,DYSAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DYSCOF(*),DYSAB(*)
      INTEGER :: IOFFA(8)
      REAL*8 GAA,GBB,OVLP
      INTEGER IORB,ISORB
      real*8, Allocatable:: DYSCOF2(:)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
!+++BRN Create a scalar spin summed Dyson coefficients DYSCOF2
!Alpha and beta contributions are added up here
      Call mma_allocate(DYSCOF2,NASHT,Label='DYSCOF2')
      DO IORB=1,NASHT
       ISORB=2*IORB-1
       GAA=DYSCOF(ISORB)
       GBB=DYSCOF(ISORB+1)
       OVLP=GBB+GAA
       !normally GAA gives just zeros...
       DYSCOF2(IORB)=OVLP
      END DO
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      END DO

C CONTRIBUTION FROM INACTIVE ORBITALS:
C (By definition 0 for Dyson orbitals,
C but we need to fill out the full vector for easier
C transformation.)
        IOFFTD=0
        DO 50 ISY=1,NSYM
         IF(NISH(ISY).NE.0) THEN
          II=0
          DO 40 I=1,NISH(ISY)
            II=II+1
            IPOS=IOFFTD+II
            DYSAB(IPOS)=0.0D0
40        CONTINUE
          IOFFTD=IOFFTD+NOSH(ISY)
         END IF
50      CONTINUE
C THEN ADD CONTRIBUTION FROM ACTIVE SPACE.
      IOFFTD=0
      ICOFF=1
      DO 120 ISY1=1,NSYM
        NO1=NOSH(ISY1)
        IF(NO1.EQ.0) GOTO 120
        NA1=NASH(ISY1)
        IF(NA1.EQ.0) GOTO 110
        NI1=NISH(ISY1)
        DO 100 I=1,NA1
          II=NI1+I
          IPOS=IOFFTD+II
          DYSAB(IPOS)=DYSCOF2(ICOFF)
          ICOFF=ICOFF+1
100     CONTINUE
110     IOFFTD=IOFFTD+NO1
120   CONTINUE
      Call mma_deallocate(DYSCOF2)
      RETURN
      END
