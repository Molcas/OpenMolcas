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
      SUBROUTINE CHO_MEM(TEXT,JOB,TYPE,KPOINT,LENGTH)
C
C     Purpose: memory management using getmem. Only differences are that
C     'flush' flushes back to and including KPOINT, and 'max ' allocates
C     max. memory available. To get max. memory available without
C     allocating, specify JOB='GETM'.
C
      IMPLICIT NONE
      CHARACTER*(*) TEXT,JOB,TYPE
      INTEGER       KPOINT, LENGTH
      INTEGER       LJOB, I, LENDUM
      CHARACTER*4   MYJOB

      LJOB = LEN(JOB)
      IF (LJOB .LT. 4) THEN
         DO I = 1,LJOB
            MYJOB(I:I) = JOB(I:I)
         END DO
         DO I = LJOB+1,4
            MYJOB(I:I) = ' '
         END DO
      ELSE
         DO I = 1,4
            MYJOB(I:I) = JOB(I:I)
         END DO
      END IF
      CALL UPCASE(MYJOB)

      IF (MYJOB(1:4) .EQ. 'MAX ') THEN
         CALL GETMEM(TEXT,'MAX ',TYPE,KPOINT,LENGTH)
         CALL GETMEM(TEXT,'ALLO',TYPE,KPOINT,LENGTH)
      ELSE IF (MYJOB(1:4) .EQ. 'FLUS') THEN
         Write (6,*) 'CHO_MEM: keyword  FLUSH is now disabled'
         Call Abend()
         LENDUM = -1
         CALL GETMEM(TEXT,'FLUSH',TYPE,KPOINT,LENDUM)
         CALL GETMEM(TEXT,'FREE',TYPE,KPOINT,LENGTH)
      ELSE IF (MYJOB(1:4) .EQ. 'GETM') THEN
         CALL GETMEM(TEXT,'MAX ',TYPE,KPOINT,LENGTH)
      ELSE
         CALL GETMEM(TEXT,JOB,TYPE,KPOINT,LENGTH)
      END IF

      END
