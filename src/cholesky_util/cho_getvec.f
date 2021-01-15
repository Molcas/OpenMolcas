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
      SUBROUTINE CHO_GETVEC(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                      SCR,LSCR)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set. The algorithm used for
C              reading is taken from input (via cholesky.fh header
C              file).
C
      use ChoArr, only: iScr
#include "implicit.fh"
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      DIMENSION SCR(LSCR)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      external ddot_

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (N2 = INFVEC_N2)

      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

C     Return if no vectors requested.
C     -------------------------------

      IF (NUMVEC .LT. 1) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
            WRITE(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
         END IF
         RETURN
      END IF

C     Check vector dimension: should be identical to current reduced
C     set. Check also symmetry and vector index.
C     --------------------------------------------------------------

      IF (LOCDBG) THEN
         IFAIL = 0
         IF (LENVEC.NE.NNBSTR(ISYM,2) .OR. LENVEC.LT.1) THEN
            WRITE(LUPRI,*) SECNAM,': illegal vector dimension:'
            WRITE(LUPRI,*) SECNAM,': LENVEC = ',LENVEC
            IFAIL = IFAIL + 1
         END IF
         IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
            WRITE(LUPRI,*) SECNAM,': illegal symmetry input:'
            WRITE(LUPRI,*) SECNAM,': ISYM = ',ISYM
            IFAIL = IFAIL + 1
         END IF
         IVEC2 = IVEC1 + NUMVEC - 1
         IF (IVEC1.LT.1 .OR. IVEC1.GT.MAXVEC .OR.
     &       IVEC2.LT.1 .OR. IVEC2.GT.MAXVEC) THEN
            WRITE(LUPRI,*) SECNAM,': illegal vector index:'
            WRITE(LUPRI,*) SECNAM,': IVEC1,IVEC2 = ',IVEC1,IVEC2
            IFAIL = IFAIL + 1
         ELSE
            IF (INFVEC(IVEC1,3,ISYM) .LT. 0) THEN
               WRITE(LUPRI,*) SECNAM,': illegal first vector address:'
               WRITE(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC1,3,ISYM)
               IFAIL = IFAIL + 1
            END IF
            IF (INFVEC(IVEC2,3,ISYM) .LT. 0) THEN
               WRITE(LUPRI,*) SECNAM,': illegal last vector address:'
               WRITE(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC2,3,ISYM)
               IFAIL = IFAIL + 1
            END IF
         END IF
         IF (CHO_IOVEC.EQ.1 .OR. CHO_IOVEC.EQ.2 .OR. CHO_IOVEC.EQ.3 .OR.
     &       CHO_IOVEC.EQ.4) THEN
            IF (SIZE(ISCR) .LT. NNBSTR(ISYM,2)) THEN
               WRITE(LUPRI,*) SECNAM,': insufficient iscratch:'
               WRITE(LUPRI,*) SECNAM,': SIZE(ISCR) = ',SIZE(ISCR)
               IFAIL = IFAIL + 1
            END IF
         END IF
         IF (IFAIL .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': unable to continue!'
            CALL CHO_QUIT('Error in '//SECNAM,104)
         END IF
      END IF

C     Call reading routine.
C     ---------------------

      IF (CHO_IOVEC .EQ. 1) THEN
         CALL CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
      ELSE IF (CHO_IOVEC.EQ.2 .OR. CHO_IOVEC.EQ.3 .OR. CHO_IOVEC.EQ.4)
     & THEN
         CALL CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
      ELSE
         CALL CHO_GETVEC0(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
      END IF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         CALL CHO_FLUSH(LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,
     &                  ' of symmetry ',ISYM,' read from unit ',
     &                  LUCHO(ISYM)
         WRITE(LUPRI,*) 'Vector dimension: ',LENVEC,
     &                  ' (current reduced set)'
         WRITE(LUPRI,*) 'Algorithm: ',CHO_IOVEC
         DO IVEC = 1,NUMVEC
            JVEC = IVEC1 + IVEC - 1
            JADR = INFVEC(JVEC,3,ISYM)
            XNRM = SQRT(DDOT_(LENVEC,CHOVEC(1,IVEC),1,CHOVEC(1,IVEC),1))
            WRITE(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',
     &                     XNRM
         END DO
         CALL CHO_FLUSH(LUPRI)
      END IF

      END
