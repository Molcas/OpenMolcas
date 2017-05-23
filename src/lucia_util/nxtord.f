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
* Copyright (C) 1989, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE NXTORD(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An ordered set of numbers INUM(I),I=1,NELMNT is
* given in strictly ascending order. Values of INUM(*) is
* restricted to the interval MINVAL,MAXVAL .
*
* Find next higher number.
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen May 1989
*
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTORD '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF( IPLACE .LT. NELMNT .AND.
     &      INUM(IPLACE)+1 .LT. INUM(IPLACE+1)
     &  .OR.IPLACE.EQ. NELMNT .AND.
     &      INUM(IPLACE)+1.LE.MAXVAL) THEN
              INUM(IPLACE) = INUM(IPLACE) + 1
              NONEW = 0
              GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
              IF(IPLACE .EQ. 1 ) THEN
                INUM(IPLACE) = MINVAL
              ELSE
                INUM(IPLACE) = INUM(IPLACE-1) + 1
              END IF
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
