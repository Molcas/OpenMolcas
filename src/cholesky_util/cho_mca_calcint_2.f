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
      SUBROUTINE CHO_MCA_CALCINT_2(ISHLAB)
C
C     Purpose: calculate qualified integral columns from
C              shell pair distribution (**|ISHLA ISHLB).
C
C     Version 2: avoid storage of full shell quadruple in interface to
C                seward; get qualified directly!
C
      use ChoArr, only: iSP2F, MySP
      use ChoSwp, only: nnBstRSh
      use Constants
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "choprint.fh"
#include "stdalloc.fh"

      CHARACTER(LEN=17), PARAMETER:: SECNAM = 'CHO_MCA_CALCINT_2'

      LOGICAL   DOINTS
      LOGICAL, PARAMETER:: LOCDBG = .FALSE.
      Integer, PARAMETER:: INFINT = INF_INT, INFIN2 = INF_IN2

      INTEGER NAB(8)

      Real*8, Allocatable:: IntCol(:)

#if defined (_DEBUGPRINT_)
      Call mma_maxDBLE(MEM_START)
#endif

C     Initializations.
C     ----------------

      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)

      NAB(1) = NQUAL(1) - IOFFQ(1)
      DO ISYM = 2,NSYM
         NAB(ISYM) = NQUAL(ISYM) - IOFFQ(ISYM)
      END DO

      IOFF_COL(1) = 0
      LCOL    = NNBSTR(1,2)*NAB(1)
      DO ISYM = 2,NSYM
         IOFF_COL(ISYM) = LCOL
         LCOL = LCOL + NNBSTR(ISYM,2)*NAB(ISYM)
      END DO

      XXSHL = DBLE(NNSHL)
      XSKIP = Zero

      IF (IPRINT .GE. INFINT) WRITE(LUPRI,*)

C     Allocate memory and initialize:
C     qualified columns in reduced set,
C     max. shell quadruple.
C     ---------------------------------

      Call mma_allocate(IntCol,LCOL,Label='IntCol')
      IntCol(:)=Zero

C     Set mapping from shell pair AB to qualified columns.
C     ----------------------------------------------------

      IRC  = 0
      ILOC = 2
      CALL CHO_P_SETSHP2Q(IRC,ILOC,ISHLAB,NAB)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2Q returned ',IRC
         CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
      END IF

C     Set memory used by seward.
C     --------------------------

      Call mma_maxDBLE(LINT)
      CALL XSETMEM_INTS(LINT)

C     Loop over shell quadruples.
C     ---------------------------

      DO ISHLCD = 1,NNSHL

C        Set left shell pair index.
C        --------------------------

         ISCD = MYSP(ISHLCD)
         CALL CHO_INVPCK(ISP2F(ISCD),ISHLC,ISHLD,.TRUE.)

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
            CALL CHO_MCA_INT_1(ISCD,ISHLAB,
     &                         IntCol,LCOL,
     &                         LOCDBG.OR.(IPRINT.GE.100))
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

         ELSE

C           Update skip counter.
C           --------------------

            XSKIP = XSKIP + One

C           Print message.
C           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

         END IF

      END DO

C     Write the columns to disk.
C     --------------------------

      CALL CHO_TIMER(C1,W1)
      DO ISYM = 1,NSYM
         LTOT = NNBSTR(ISYM,2)*NAB(ISYM)
         IF (LTOT .GT. 0) THEN
            IOPT = 1
            KOFF = 1 + IOFF_COL(ISYM)
            IADR = NNBSTR(ISYM,2)*IOFFQ(ISYM)
            CALL DDAFILE(LUSEL(ISYM),IOPT,IntCol(KOFF),LTOT,IADR)
         END IF
      END DO
      CALL CHO_TIMER(C2,W2)
      TINTEG(1,2) = TINTEG(1,2) + C2 - C1
      TINTEG(2,2) = TINTEG(2,2) + W2 - W1

C     Free memory: both memory used by seward and used here.
C     ------------------------------------------------------

      CALL XRLSMEM_INTS()
      Call mma_deallocate(IntCol)

C     Print skip statistics.
C     ----------------------

      IF (IPRINT .GE. INFIN2) THEN
         PCT = 1.0D2*XSKIP/XXSHL
         WRITE(LUPRI,'(A,F7.2,A)')
     &   'Skipped',PCT,'% of rows (shell pairs) in this distribution'
      END IF

#if defined (_DEBUGPRINT_)
      Call mma_maxDBLE(MEM_END)
      LEAK = MEM_END - MEM_START
      IF (LEAK .NE. 0) THEN
         WRITE(LUPRI,'(//,A,A,I9)')
     &   SECNAM,': Memory leak:',LEAK
         CALL CHO_QUIT('Memory leak detected in '//SECNAM,104)
      END IF
#endif

      END
