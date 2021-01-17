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
      SUBROUTINE CHO_MCA_CALCINT_3(XINT,LINT,ISHLAB)
C
C     Purpose: calculate qualified integral columns from
C              shell pair distribution (**|ISHLA ISHLB).
C
C     Version 3: avoid storage of full shell quadruple in interface to
C                seward; get qualified directly as in Version 2!
C                Changes from Version 2:
C                - addressing of qualified columns
C                - integrals returned in core (no I/O)
C
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh
#include "implicit.fh"
      DIMENSION XINT(LINT)
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT_3')

      INTEGER NAB(8)
      INTEGER CHO_ISUMELM

      LOGICAL   DOINTS, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (INFINT = INF_INT, INFIN2 = INF_IN2)

C     Initializations.
C     ----------------

      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)

      XXSHL = DBLE(NNSHL)
      XSKIP = 0.0D0

      IF (IPRINT .GE. INFINT) WRITE(LUPRI,*)

C     Set mapping from shell pair AB to qualified columns.
C     ----------------------------------------------------

      IRC  = 0
      ILOC = 2
      CALL CHO_SETSHP2Q_2(IRC,ILOC,ISHLAB,NAB)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2Q_2 returned ',IRC
         CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
      END IF

C     Print.
C     ------

      IF (IPRINT .GE. INF_IN2) THEN
         NCOLAB = CHO_ISUMELM(NAB,NSYM)
         WRITE(LUPRI,'(/,A,I5,1X,I5,A,I9,A)')
     &   'Calculating shell pair (**|',ISHLA,ISHLB,
     &   '):',NCOLAB,' columns have been qualified'
         WRITE(LUPRI,'(80A)') ('=',i=1,77)
      END IF

C     Loop over shell quadruples.
C     ---------------------------

      DO ISHLCD = 1,NNSHL

C        Set left shell pair index.
C        --------------------------

         CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)

C        Find out if this shell pair (CD) contributes to
C        current reduced set.
C        -----------------------------------------------

         ISYM   = 1
         DOINTS = (NAB(ISYM).GT.0) .AND.
     &            (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         DO WHILE ((ISYM.LT.NSYM) .AND. (.NOT.DOINTS))
            ISYM   = ISYM + 1
            DOINTS = (NAB(ISYM).GT.0) .AND.
     &               (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         END DO

         IF (DOINTS) THEN

C           Print message.
C           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'Invoking Seward for shell quadruple (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

C           Set mapping from shell pair CD to reduced set.
C           ----------------------------------------------

            IRC  = 0
            ILOC = 2
            CALL CHO_SETSHP2RS(IRC,ILOC,ISHLCD,NAB)
            IF (IRC .NE. 0) THEN
               WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2RS returned ',IRC
               CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
            END IF

C           Calculate integrals.
C           --------------------

            CALL CHO_TIMER(C1,W1)
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,
     &                         XINT,LINT,
     &                         LOCDBG.OR.(IPRINT.GE.100))
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

         ELSE

C           Update skip counter.
C           --------------------

            XSKIP = XSKIP + 1.0D0

C           Print message.
C           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

         END IF

      END DO

C     Print skip statistics.
C     ----------------------

      IF (IPRINT .GE. INFIN2) THEN
         PCT = 1.0D2*XSKIP/XXSHL
         WRITE(LUPRI,'(A,F7.2,A)')
     &   'Skipped',PCT,'% of rows (shell pairs) in this distribution'
         CALL CHO_FLUSH(LUPRI)
      END IF

      END
