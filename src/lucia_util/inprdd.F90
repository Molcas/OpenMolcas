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
      REAL*8 FUNCTION INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
!
! DISC VERSION OF INPROD
!
! LBLK DEFINES STRUCTURE OF FILE
!
      use Constants, only: Zero
      use lucia_data, only: IDISK
      IMPLICIT NONE
      REAL*8 VEC1(*),VEC2(*)
      INTEGER LU1,LU2,IREW,LBLK

      REAL*8 INPROD,X
      LOGICAL DIFVEC
      INTEGER IDUMMY(1),NBL1,KBLK,IAMPACK,IMZERO
!
      X = Zero
      IF( LU1 .NE. LU2 ) THEN
        DIFVEC = .TRUE.
      ELSE
        DIFVEC =  .FALSE.
      END IF
!
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          IDISK(LU1)=0
          IF(DIFVEC) IDISK(LU2)=0
         ELSE
          IDISK(LU1)=0
          IF(DIFVEC) IDISK(LU2)=0
         END IF
      END IF
!
! LOOP OVER BLOCKS OF VECTORS
!
 1000 CONTINUE
!
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          CALL IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
          NBL1=IDUMMY(1)
          IF( DIFVEC) THEN
            CALL IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
          END IF
        ELSE IF ( LBLK .LT. 0 ) THEN
          CALL IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
          NBL1=IDUMMY(1)
          CALL IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
          IF( DIFVEC) THEN
            CALL IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
            CALL IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
          END IF
        END IF
!
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC(     VEC1,     NBL1,     KBLK,      LU1,   IMZERO,&
     &                  IAMPACK)
          IF( DIFVEC) THEN
            CALL FRMDSC(    VEC2,    NBL1,    KBLK,     LU2,  IMZERO,   &
     &                   IAMPACK)
            IF(NBL1 .GT. 0 )                                            &
     &      X = X + INPROD(VEC1,VEC2,NBL1)
!?          write(6,*) ' vec1 and vec2 in INPRDD '
!?         CALL WRTMAT(VEC1,1,NBL1,1,NBL1)
!?         CALL WRTMAT(VEC2,1,NBL1,1,NBL1)
          ELSE
          IF(NBL1 .GT. 0 )                                              &
     &    X = X + INPROD(VEC1,VEC1,NBL1)
        END IF
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
!
      INPRDD = X
!
      END FUNCTION INPRDD
