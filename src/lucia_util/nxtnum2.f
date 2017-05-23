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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An set of numbers INUM(I),I=1,NELMNT is
* given. Find next compund number.
* Digit I must be in the range MINVAL,MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
C08-08-01   INUM(JPLACE) = 1
            INUM(JPLACE) = MINVAL
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
