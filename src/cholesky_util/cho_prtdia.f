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
      SUBROUTINE CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
C
C     Purpose: print requested symmetry block(s) of diagonal in
C              first (IRED=1) or current (IRED=2) reduced set.
C
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh
#include "implicit.fh"
      DIMENSION DIAG(*)
      INTEGER   ISYLST(NSYLST)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_PRTDIA')

      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)

C     Check dimension of symmetry list.
C     ---------------------------------

      IF (NSYLST .LT. 1) THEN
         RETURN
      ELSE IF (NSYLST .GT. NSYM) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': NSYLST <= NSYM required!'
         WRITE(LUPRI,'(1X,A,I10)')    'NSYLST = ',NSYLST
         WRITE(LUPRI,'(1X,A,I10,/)')  'NSYM   = ',NSYM
         CALL CHO_QUIT('[0] Symmetry error in '//SECNAM,102)
      END IF

C     Code for first or current reduced set.
C     --------------------------------------

      IF (IRED .EQ. 1) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal in Original Reduced Set',
     &                 '=',80,LUPRI)
         DO ILST = 1,NSYLST
            ISYM = ISYLST(ILST)
            IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
               WRITE(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,
     &                        ' of list ISYLST is out of bounds!'
               CALL CHO_QUIT('ISYLST input error in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
               WRITE(LUPRI,'(/,A,/,A)')
     &         '  Element Shell-Pair  SP Index         Diagonal',
     &         '-----------------------------------------------'
               DO ISHLAB = 1,NNSHL
                  IAB1 = IIBSTR(ISYM,IRED) + IIBSTRSH(ISYM,ISHLAB,IRED)
     &                 + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,IRED) - 1
                  DO IAB = IAB1,IAB2
                     IF (INDRSH(IAB) .NE. ISP2F(ISHLAB)) THEN
                        WRITE(LUPRI,*)
     &                  'Shell Pair error: INDRSH,ISP2F,ISHLAB',
     &                  INDRSH(IAB),ISP2F(ISHLAB),ISHLAB
                        CALL CHO_QUIT('Shell-Pair error in '//SECNAM,
     &                                104)
                     ELSE
                        JAB = IAB
                        WRITE(LUPRI,'(I9,2X,I9,1X,I9,1X,1P,D16.8)')
     &                  JAB,ISP2F(ISHLAB),INDRED(IAB,IRED),DIAG(JAB)
                     END IF
                  END DO
               END DO
               WRITE(LUPRI,'(A)')
     &         '-----------------------------------------------'
            END IF
         END DO
      ELSE IF (IRED .EQ. 2) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal in Current Reduced Set',
     &                 '=',80,LUPRI)
         DO ILST = 1,NSYLST
            ISYM = ISYLST(ILST)
            IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
               WRITE(LUPRI,*) SECNAM,': element ',ILST,': ',ISYM,
     &                        ' of list ISYLST is out of bounds!'
               CALL CHO_QUIT('ISYLST input error in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,'(/,A,I2)') 'Symmetry block:',ISYM
               WRITE(LUPRI,'(/,A,A,/,A,A)')
     &         '  Element  RedSet 1 Shell-Pair  SP Index',
     &         '         Diagonal',
     &         '----------------------------------------',
     &         '-----------------'
               DO ISHLAB = 1,NNSHL
                  IAB1 = IIBSTR(ISYM,IRED) + IIBSTRSH(ISYM,ISHLAB,IRED)
     &                 + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,IRED) - 1
                  DO IAB = IAB1,IAB2
                     JAB = INDRED(IAB,IRED)
                     IF (INDRSH(JAB) .NE. ISP2F(ISHLAB)) THEN
                        WRITE(LUPRI,*)
     &                  'Shell Pair error: INDRSH,ISP2F,ISHLAB',
     &                  INDRSH(JAB),ISP2F(ISHLAB),ISHLAB
                        CALL CHO_QUIT('Shell-Pair error in '//SECNAM,
     &                                104)
                     ELSE
                       WRITE(LUPRI,'(I9,1X,I9,2X,I9,1X,I9,1X,1P,D16.8)')
     &                  IAB,JAB,ISP2F(ISHLAB),INDRED(JAB,1),DIAG(JAB)
                     END IF
                  END DO
               END DO
               WRITE(LUPRI,'(A,A)')
     &         '----------------------------------------',
     &         '-----------------'
            END IF
         END DO
      END IF

      END
