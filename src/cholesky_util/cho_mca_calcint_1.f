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
      SUBROUTINE CHO_MCA_CALCINT_1(ISHLAB)
!
!     Purpose: calculate qualified integral columns from
!              shell pair distribution (**|ISHLA ISHLB).
!
!     Version 1: store full shell quadruple.
!
      use ChoArr, only: nBstSh, iSP2F
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRed
      use Constants
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "choprint.fh"
#include "stdalloc.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT_1')

      LOGICAL   DOINTS, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (INFINT = INF_INT, INFIN2 = INF_IN2)

      INTEGER NAB(8)

      Real*8, Allocatable:: Int4Sh(:), IntCol(:)

#if defined (_DEBUGPRINT_)
      Call mma_maxDBLE(LLEAK)
      MEM_START = LLEAK
#endif

!     Initializations.
!     ----------------

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
      XSKIP = Zero

      IF (IPRINT .GE. INFINT) WRITE(LUPRI,*)

!     Allocate memory and initialize:
!     qualified columns in reduced set,
!     max. shell quadruple.
!     ---------------------------------

      Call mma_allocate(Int4Sh,L4SHMX,Label='Int4Sh')
      Call mma_allocate(IntCol,LCOL,Label='IntCol')
      IntCol(:)=Zero

!     Set memory used by seward.
!     --------------------------

      Call mma_maxDBLE(LINT)
      CALL XSETMEM_INTS(LINT)

!     Loop over shell quadruples.
!     ---------------------------

      DO ISHLCD = 1,NNSHL

!        Set left shell pair index.
!        --------------------------

         CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
         IF (ISHLC .EQ. ISHLD) THEN
            NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
         ELSE
            NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
         END IF

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

!           Calculate integrals.
!           --------------------

            CALL CHO_TIMER(C1,W1)
            L4SH = NUMCD*NUMAB
            Int4Sh(1:L4SH)=Zero
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,
     &                         Int4SH,L4SH,
     &                         LOCDBG.OR.(IPRINT.GE.100))
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

!           Extract columns in reduced set.
!           IAB: index AB within full shell pair.
!           JAB: index AB within current reduced set.
!           KAB: index AB within qualifieds.
!           -----------------------------------------

            DO ISYM = 1,NSYM
               DO KAB = 1,NAB(ISYM)

                  JAB = IQUAB(IOFFQ(ISYM)+KAB,ISYM)
                  IAB = INDRED(INDRED(JAB,2),1)

                  DO JCD0 = 1,NNBSTRSH(ISYM,ISHLCD,2)

                     JCDS = IIBSTRSH(ISYM,ISHLCD,2) + JCD0
                     JCD  = IIBSTR(ISYM,2) + JCDS
                     ICD  = INDRED(INDRED(JCD,2),1)

                     KOFF1 = IOFF_COL(ISYM)
     &                     + NNBSTR(ISYM,2)*(KAB - 1) + JCDS
                     KOFF2 = NUMCD*(IAB - 1) + ICD

                     IntCol(KOFF1) = Int4Sh(KOFF2)

                  END DO

               END DO
            END DO

         ELSE

!           Update skip counter.
!           --------------------

            XSKIP = XSKIP + One

!           Print message.
!           --------------

            IF (IPRINT .GE. INFINT) THEN
                WRITE(LUPRI,'(A,I5,1X,I5,A,I5,1X,I5,A)')
     &          'NOTICE: skipping shell quadruple    (',ISHLC,ISHLD,
     &          '|',ISHLA,ISHLB,')'
            END IF

         END IF

      END DO

!     Write the columns to disk.
!     --------------------------

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

!     Free memory: both memory used by seward and used here.
!     ------------------------------------------------------

      CALL XRLSMEM_INTS()
      Call mma_deallocate(IntCol)
      Call mma_deallocate(Int4Sh)

!     Print skip statistics.
!     ----------------------

      IF (IPRINT .GE. INFIN2) THEN
         PCT = 1.0D2*XSKIP/XXSHL
         WRITE(LUPRI,'(A,F7.2,A)')
     &   'Skipped',PCT,'% of rows (shell pairs) in this distribution'
      END IF

#if defined (_DEBUGPRINT_)
      Call mma_maxDBLE(LLEAK)
      MEM_END = LLEAK
      LEAK = MEM_END - MEM_START
      IF (LEAK .NE. 0) THEN
         WRITE(LUPRI,'(//,A,A,I9)')
     &   SECNAM,': Memory leak:',LEAK
         CALL CHO_QUIT('Memory leak detected in '//SECNAM,104)
      END IF
#endif

      END
