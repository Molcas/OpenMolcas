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
* Copyright (C) Jeppe Olsen                                            *
************************************************************************
      SUBROUTINE BLKCLS(   IBLKS,   NBLKS, IBLKCLS, ISPSPCL,    NCLS,
     &                      LCLS,  NOCTPA,  NOCTPB)
*
* Class of each block, and dimension of each class
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLKS(8,NBLKS)
      INTEGER ISPSPCL(NOCTPA,NOCTPB)
*. Output
      INTEGER IBLKCLS(NBLKS),LCLS(NCLS)
*
C?    WRITE(6,*) ' ISPSPCL'
C?    CALL IWRTMA(ISPSPCL,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      IZERO = 0
      CALL ISETVC(LCLS,IZERO,NCLS)
      DO JBLK = 1, NBLKS
        IICLS = ISPSPCL(IBLKS(1,JBLK),IBLKS(2,JBLK))
        IBLKCLS(JBLK) = IICLS
        LCLS(IICLS) = LCLS(IICLS) + IBLKS(8,JBLK)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' BLKCLS Speaking '
        WRITE(6,*) ' ==============='
        WRITE(6,*)
        WRITE(6,*) ' Dimension of each class '
        CALL IWRTMA(LCLS,1,NCLS,1,NCLS)
        WRITE(6,*)
        WRITE(6,*) ' Class of each block : '
        CALL IWRTMA(IBLKCLS,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
