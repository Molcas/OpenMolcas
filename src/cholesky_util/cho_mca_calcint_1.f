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
      SUBROUTINE CHO_MCA_CALCINT_1(ISHLAB)
C
C     Purpose: calculate qualified integral columns from
C              shell pair distribution (**|ISHLA ISHLB).
C
C     Version 1: store full shell quadruple.
C
      use ChoArr, only: nBstSh, iSP2F
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh
#include "implicit.fh"
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT_1')

      LOGICAL   DOINTS, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (INFINT = INF_INT, INFIN2 = INF_IN2)

      INTEGER NAB(8)

      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)

#if defined (_DEBUGPRINT_)
      CALL GETMEM('INT.LEAK','MAX ','REAL',KLEAK,LLEAK)
      MEM_START = LLEAK
#endif

C     Initializations.
C     ----------------

      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)

      NAB(1) = NQUAL(1) - IOFFQ(1)
      DO ISYM = 2,NSYM
         NAB(ISYM) = NQUAL(ISYM) - IOFFQ(ISYM)
      END DO

      IOFF_COL(1) = 0
      LCOL = NNBSTR(1,2)*NAB(1)
      DO ISYM = 2,NSYM
         IOFF_COL(ISYM) = LCOL
         LCOL = LCOL + NNBSTR(ISYM,2)*NAB(ISYM)
      END DO

      IF (ISHLA .EQ. ISHLB) THEN
         NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLA) + 1)/2
      ELSE
         NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
      END IF
      MAXCD = 0
      DO ISHLCD = 1,NNSHL
         ISYM   = 1
         DOINTS = (NAB(ISYM).GT.0) .AND.
     &            (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         DO WHILE ((ISYM.LT.NSYM) .AND. (.NOT.DOINTS))
            ISYM   = ISYM + 1
            DOINTS = (NAB(ISYM).GT.0) .AND.
     &               (NNBSTRSH(ISYM,ISHLCD,2).GT.0)
         END DO
         IF (DOINTS) THEN
            CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
            IF (ISHLC .EQ. ISHLD) THEN
               NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
            ELSE
               NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
            END IF
            MAXCD = MAX(MAXCD,NUMCD)
         END IF
      END DO
      L4SHMX = MAXCD*NUMAB

      XXSHL = DBLE(NNSHL)
      XSKIP = 0.0D0

      IF (IPRINT .GE. INFINT) WRITE(LUPRI,*)

C     Allocate memory and initialize:
C     qualified columns in reduced set,
C     max. shell quadruple.
C     ---------------------------------

      CALL GETMEM('INT.4sh','ALLO','REAL',K4SH,L4SHMX)
      CALL GETMEM('INT.col','ALLO','REAL',KCOL,LCOL)
      CALL CHO_DZERO(WORK(KCOL),LCOL)

C     Set memory used by seward.
C     --------------------------

      CALL GETMEM('INT.MAX','MAX ','REAL',KINT,LINT)
      CALL XSETMEM_INTS(LINT)

C     Loop over shell quadruples.
C     ---------------------------

      DO ISHLCD = 1,NNSHL

C        Set left shell pair index.
C        --------------------------

         CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
         IF (ISHLC .EQ. ISHLD) THEN
            NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
         ELSE
            NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
         END IF

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

C           Calculate integrals.
C           --------------------

            CALL CHO_TIMER(C1,W1)
            L4SH = NUMCD*NUMAB
            CALL CHO_DZERO(WORK(K4SH),L4SH)
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,
     &                         WORK(K4SH),L4SH,
     &                         LOCDBG.OR.(IPRINT.GE.100))
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

C           Extract columns in reduced set.
C           IAB: index AB within full shell pair.
C           JAB: index AB within current reduced set.
C           KAB: index AB within qualifieds.
C           -----------------------------------------

            DO ISYM = 1,NSYM
               DO KAB = 1,NAB(ISYM)

                  JAB = IQUAB(IOFFQ(ISYM)+KAB,ISYM)
                  IAB = INDRED(INDRED(JAB,2),1)

                  DO JCD0 = 1,NNBSTRSH(ISYM,ISHLCD,2)

                     JCDS = IIBSTRSH(ISYM,ISHLCD,2) + JCD0
                     JCD  = IIBSTR(ISYM,2) + JCDS
                     ICD  = INDRED(INDRED(JCD,2),1)

                     KOFF1 = KCOL + IOFF_COL(ISYM)
     &                     + NNBSTR(ISYM,2)*(KAB - 1) + JCDS - 1
                     KOFF2 = K4SH + NUMCD*(IAB - 1) + ICD - 1

                     WORK(KOFF1) = WORK(KOFF2)

                  END DO

               END DO
            END DO

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

C     Write the columns to disk.
C     --------------------------

      CALL CHO_TIMER(C1,W1)
      DO ISYM = 1,NSYM
         LTOT = NNBSTR(ISYM,2)*NAB(ISYM)
         IF (LTOT .GT. 0) THEN
            IOPT = 1
            KOFF = KCOL + IOFF_COL(ISYM)
            IADR = NNBSTR(ISYM,2)*IOFFQ(ISYM)
            CALL DDAFILE(LUSEL(ISYM),IOPT,WORK(KOFF),LTOT,IADR)
         END IF
      END DO
      CALL CHO_TIMER(C2,W2)
      TINTEG(1,2) = TINTEG(1,2) + C2 - C1
      TINTEG(2,2) = TINTEG(2,2) + W2 - W1

C     Free memory: both memory used by seward and used here.
C     ------------------------------------------------------

      CALL XRLSMEM_INTS
      CALL GETMEM('INT.col','FREE','REAL',KCOL,LCOL)
      CALL GETMEM('INT.4sh','FREE','REAL',K4SH,L4SHMX)

C     Print skip statistics.
C     ----------------------

      IF (IPRINT .GE. INFIN2) THEN
         PCT = 1.0D2*XSKIP/XXSHL
         WRITE(LUPRI,'(A,F7.2,A)')
     &   'Skipped',PCT,'% of rows (shell pairs) in this distribution'
      END IF

#if defined (_DEBUGPRINT_)
      CALL GETMEM('INT.LEAK','MAX ','REAL',KLEAK,LLEAK)
      MEM_END = LLEAK
      LEAK = MEM_END - MEM_START
      IF (LEAK .NE. 0) THEN
         WRITE(LUPRI,'(//,A,A,I9)')
     &   SECNAM,': Memory leak:',LEAK
         CALL CHO_QUIT('Memory leak detected in '//SECNAM,104)
      END IF
#endif

      END
