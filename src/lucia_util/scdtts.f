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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE SCDTTS(  BLOCKS,  IBLOCK,  NBLOCK,   NSMST,
     &                      NSASO,   NSBSO,     IDC,    IWAY,
     &                     IPRNT)
*
* Scale batch of
* blocks between determinant and combination form
*
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* The blocks are assumed to be in packed form !!
*
*. Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
*
C?    LOGICAL DIAGBL
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from SCDTTS '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTTTS(   BLOCKS,   IBLOCK,   NBLOCK,    NSMST,
     &                 NSASO,    NSBSO,        2)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IOFFP= IBLOCK(6, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*. Is this block diagonal in packed form
        IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA   = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
        IF(IPACK .EQ. 1 ) THEN
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
*Ms combinations
        IF(IDC.EQ.2) THEN
          IF(IWAY.EQ.1) THEN
            FACTOR = SQ2
          ELSE
            FACTOR = SQ2I
          END IF
          CALL SCALVE(BLOCKS(IOFFP),FACTOR,NELMNT)
          IF(IPACK.EQ.1 ) THEN
            FACTOR = 1.0D0/FACTOR
            CALL SCLDIA(BLOCKS(IOFFP),FACTOR,NIA,1)
          END IF
        END IF
*
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output vector '
        CALL WRTTTS(   BLOCKS,   IBLOCK,   NBLOCK,    NSMST,
     &                 NSASO,    NSBSO,        2)
      END IF
*
      RETURN
      END
