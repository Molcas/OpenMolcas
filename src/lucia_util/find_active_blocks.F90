!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
!
!. Find the active (nonvanishing blocks) on LUIN
!. Non vanishing block is flagged by a 1.0 ( note : real)
!  in BLK_A
!
      use lucia_data, only: IDISK
      IMPLICIT None
      INTEGER LUIN,LBLK
!. Output
      REAL*8 BLK_A(*)
!. Scratch
      REAL*8 SEGMNT(*)

      INTEGER LBL(1),IDUMMY(1)
      INTEGER IBLK,NBLK_A,KBLK,NO_ZEROING,IMZERO,NBLK,NTEST,IAMPACK
!
      IDISK(LUIN)=0
!
      IBLK = 0
      NBLK_A = 0
!. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL(1) = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          CALL IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
          CALL IDAFILE(LUIN,2,IDUMMY,1,IDISK(LUIN))
        END IF
        IF( LBL(1) .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL(1)
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(  SEGMNT,  LBL(1),    KBLK,    LUIN,  IMZERO,    &
     &                  IAMPACK,NO_ZEROING)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
      IF( LBL(1) .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
      NBLK =  IBLK-1
!
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*)                                                      &
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks '
        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
      END IF
!
      END SUBROUTINE FIND_ACTIVE_BLOCKS
!
! Obtain property integrals with LABEL LABEL from LU91,
! LUCIA format
!
! Jeppe Olsen, Feb.98

!
!
!
