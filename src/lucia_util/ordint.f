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
      SUBROUTINE ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
*
* ORDER A STRING OF INTEGERS TO ASCENDING ORDER
*
* IINST : INPUT STRING
* IOUTST : OUTPUT STRING
* NELMNT : NUMBER OF INTEGERS
* INO : Mapping array from new to old order
*
* THIS CODE CONTAINS THE OLD ORDER CODE OF JOE GOLAB
* ( HE IS HEREBY AKNOWLEDGED , AND I AM EXCUSED )
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IINST(NELMNT),IOUTST(NELMNT),INO(NELMNT)
*
      IF(NELMNT.EQ.0) GOTO 1001
      CALL ICOPVE(IINST,IOUTST,NELMNT)
      DO  5 I = 1, NELMNT
        INO(I) = I
    5 CONTINUE
C
C       BEGIN TO ORDER
C
        JOE = 1
  10    I = JOE
  20    CONTINUE
        IF(I.EQ.NELMNT) GO TO 50
        IF(IOUTST(I).LE.IOUTST(I+1)) GO TO 40
        JOE = I + 1
  30    ISWAP = IOUTST(I)
        IOUTST(I) = IOUTST(I+1)
        IOUTST(I+1) = ISWAP
        ISWAP = INO(I)
        INO(I) = INO(I+1)
        INO(I+1) = ISWAP
        IF(I.EQ.1) GO TO 10
        I = I - 1
        IF(IOUTST(I).GT.IOUTST(I+1)) GO TO 30
        GO TO 10
 40     I = I + 1
      GO TO 20
C
C     END ORDER
C
 50   CONTINUE
*
 1001 CONTINUE
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GE.200) THEN
        WRITE(6,*) ' Result from ORDINT '
        WRITE(6,*)
        WRITE(6,*)  ' Input string '
        CALL IWRTMA(IINST,1,NELMNT,1,NELMNT)
        WRITE(6,*)  ' Ordered string '
        CALL IWRTMA(IOUTST,1,NELMNT,1,NELMNT)
        WRITE(6,*) ' New to old order '
        CALL IWRTMA(INO,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
