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
      SUBROUTINE WRTTTS(  BLOCKS,  IBLOCK,  NBLOCK,   NSMST,
     &                    NSASO,   NSBSO,     ISC)
*
* Print a batch of TTS blocks as given by IBLOCK
*
*
* ISC = 1 : In slater determinant form
* ISC = 2 : In Combination        form
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
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF (ISC.EQ.1 ) THEN
          IOFF = IBLOCK(5,JBLOCK)
        ELSE
          IOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(IPACK.EQ.1) THEN
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL PRSM2(BLOCKS(IOFF) ,NIA)
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL WRTMAT(BLOCKS(IOFF) ,NIA,NIB,NIA,NIB)
          END IF
        END IF
*
        END IF
      END DO
*
      RETURN
      END
