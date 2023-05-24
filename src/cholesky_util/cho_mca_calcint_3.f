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
      SUBROUTINE CHO_MCA_CALCINT_3(XINT,LINT,ISHLAB)
!
!     Purpose: calculate qualified integral columns from
!              shell pair distribution (**|ISHLA ISHLB).
!
!     Version 3: avoid storage of full shell quadruple in interface to
!                seward; get qualified directly as in Version 2!
!                Changes from Version 2:
!                - addressing of qualified columns
!                - integrals returned in core (no I/O)
!
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh
      Implicit Real*8 (a-h,o-z)
      REAL*8 XINT(LINT)
#include "cholesky.fh"
#include "choprint.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT_3')

      INTEGER NAB(8)
      INTEGER CHO_ISUMELM

      LOGICAL   DOINTS, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (INFINT = INF_INT, INFIN2 = INF_IN2)

!     Initializations.
!     ----------------

      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)

      XXSHL = DBLE(NNSHL)
      XSKIP = 0.0D0

      IF (IPRINT .GE. INFINT) WRITE(LUPRI,*)

!     Set mapping from shell pair AB to qualified columns.
!     ----------------------------------------------------

      IRC  = 0
      ILOC = 2
      CALL CHO_SETSHP2Q_2(IRC,ILOC,ISHLAB,NAB)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2Q_2 returned ',IRC
         CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
      END IF

!     Print.
!     ------

      IF (IPRINT .GE. INF_IN2) THEN
         NCOLAB = CHO_ISUMELM(NAB,NSYM)
         WRITE(LUPRI,'(/,A,I5,1X,I5,A,I9,A)')
     &   'Calculating shell pair (**|',ISHLA,ISHLB,
     &   '):',NCOLAB,' columns have been qualified'
         WRITE(LUPRI,'(80A)') ('=',i=1,77)
      END IF

!     Loop over shell quadruples.
!     ---------------------------

      DO ISHLCD = 1,NNSHL

!        Set left shell pair index.
!        --------------------------

         CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)

!        Find out if this shell pair (CD) contributes to
!        current reduced set.
!        -----------------------------------------------

         ISYM   = 1
         DOINTS = (NAB(ISYM).GT.0) .AND.
     &            (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         DO WHILE ((ISYM.LT.NSYM) .AND. (.NOT.DOINTS))
            ISYM   = ISYM + 1
            DOINTS = (NAB(ISYM).GT.0) .AND.
     &               (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         END DO

         IF (DOINTS) THEN

!           Print message.
!           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'Invoking Seward for shell quadruple (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

!           Set mapping from shell pair CD to reduced set.
!           ----------------------------------------------

            IRC  = 0
            ILOC = 2
            CALL CHO_SETSHP2RS(IRC,ILOC,ISHLCD,NAB)
            IF (IRC .NE. 0) THEN
               WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2RS returned ',IRC
               CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
            END IF

!           Calculate integrals.
!           --------------------

            CALL CHO_TIMER(C1,W1)
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,
     &                         XINT,LINT,
     &                         LOCDBG.OR.(IPRINT.GE.100))
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

         ELSE

!           Update skip counter.
!           --------------------

            XSKIP = XSKIP + 1.0D0

!           Print message.
!           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

         END IF

      END DO

!     Print skip statistics.
!     ----------------------

      IF (IPRINT .GE. INFIN2) THEN
         PCT = 1.0D2*XSKIP/XXSHL
         WRITE(LUPRI,'(A,F7.2,A)')
     &   'Skipped',PCT,'% of rows (shell pairs) in this distribution'
         CALL CHO_FLUSH(LUPRI)
      END IF

      END
