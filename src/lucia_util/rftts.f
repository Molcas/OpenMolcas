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
      SUBROUTINE RFTTS( BLOCKSI, BLOCKSO,  IBLOCK,  NBLOCK,   ICOPY,
     &                    NSMST,    NSASO,   NSBSO,
     &                      IDC,      PS,    IWAY,   IPRNT)
*
* Reformat between determinant and combination form of
* matrices. No scaling is performed .
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* Combination storage mode is defined BY IDC
*
*. Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKSI(*),BLOCKSO(*)
      INTEGER IBLOCK(8,NBLOCK)
*

      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
      LENGTH = 0
      IF(IWAY.EQ.1) THEN
        ISCI = 1
        ISCO = 2
      ELSE
        ISCI = 2
        ISCO = 1
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTTTS(  BLOCKSI,   IBLOCK,   NBLOCK,    NSMST,
     &                NSASO,    NSBSO,     ISCI)
      END IF
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF(IWAY.EQ.1) THEN
          IOFFI = IBLOCK(5,JBLOCK)
          IOFFO = IBLOCK(6,JBLOCK)
        ELSE
          IOFFO = IBLOCK(5,JBLOCK)
          IOFFI = IBLOCK(6,JBLOCK)
        END IF
*. Is this block diagonal in packed form
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
*. Number of elements in output block
        IF(IPACK .EQ. 1 .AND. ISCO.EQ.2 ) THEN
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
C?     WRITE(6,*) ' JBLOCK, NELMNT = ', JBLOCK,NELMNT
C?     write(6,*)
C?   & ' RFTTS : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
C?     WRITE(6,*)
C?   & ' RFTTS : NIA NIB IOFFI,IOFFO',NIA,NIB,IOFFI,IOFFO
*
        IF(IPACK.EQ.0) THEN
*. Just copy
          CALL COPVEC(BLOCKSI(IOFFI),BLOCKSO(IOFFO),NELMNT)
        ELSE
          IF(IWAY.EQ.1) THEN
*. unpacked => packed
C TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
            CALL TRIPK3(BLOCKSI(IOFFI),BLOCKSO(IOFFO),1,NIA,NIA,
     &                        PS)
          ELSE
*. Packed => unpacked
            CALL TRIPK3(BLOCKSO(IOFFO),BLOCKSI(IOFFI),2,NIA,NIA,
     &                        PS)
          END IF
        END IF
        LENGTH = LENGTH + NELMNT
        END IF
      END DO
*
      IF(ICOPY.NE.0) THEN
        CALL COPVEC(BLOCKSO,BLOCKSI,LENGTH)
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Output vector '
        CALL WRTTTS(  BLOCKSO,   IBLOCK,   NBLOCK,    NSMST,
     &                NSASO,    NSBSO,     ISCO)
      END IF
*
      RETURN
      END
