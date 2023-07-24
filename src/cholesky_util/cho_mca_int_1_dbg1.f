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
      SUBROUTINE CHO_MCA_INT_1_DBG1(DIAG,IRED)
!
!     Purpose: test diagonal, reduced set IRED. Note that the
!              diagonal *must* be the original diagonal stored
!              in reduced set 1.
!
      use ChoArr, only: nBstSh, iSP2F
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, IndRed
      use Constants
      use stdalloc
      Implicit Real*8 (a-h,o-z)
      Integer IRED
      Real*8 DIAG(*)
#include "cholesky.fh"

      CHARACTER(LEN=18), PARAMETER:: SECNAM = 'CHO_MCA_INT_1_DBG1'

      LOGICAL, PARAMETER:: PRTINT = .FALSE.

      Real*8, Allocatable:: xINT(:)

      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) SECNAM,': testing diagonal, reduced set ',IRED
      WRITE(LUPRI,*)

!     Force computation of full shell quadruple.
!     ------------------------------------------

      IF (IFCSEW .NE. 1) THEN
         WRITE(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',
     &                  IFCSEW,' to 1.'
         IFCSEW = 1
      END IF

      LINT1 = MX2SH*MX2SH
      Call mma_allocate(xINT,LINT1,Label='INT')
      Call mma_maxDBLE(LSEW)
      CALL XSETMEM_INTS(LSEW)

      NERR = 0
      NTST = 0
      DO ISHLAB = 1,NNSHL

!        Allocate memory for shell quadruple (AB|AB).
!        --------------------------------------------

         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         IF (ISHLB .EQ. ISHLA) THEN
            NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLB) + 1)/2
         ELSE
            NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
         END IF
         LINT = NUMAB*NUMAB

!        Calculate integrals.
!        --------------------

         xINT(1:LINT)=Zero
         CALL CHO_MCA_INT_1(ISHLAB,ISHLAB,xINT,LINT,PRTINT)

!        Look up all diagonal elements in DIAG and compare to
!        values just calculated.
!        ----------------------------------------------------

         IERR = 0
         IF (IRED .EQ. 1) THEN

            DO ISYM = 1,NSYM

               JAB1 = IIBSTR(ISYM,1) + IIBSTRSH(ISYM,ISHLAB,1)
     &              + 1
               JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,1) - 1

               DO JAB = JAB1,JAB2   ! loop over elements in diagonal

                  IF ((JAB.LT.1) .OR. (JAB.GT.NNBSTRT(1))) THEN
                     WRITE(LUPRI,*) SECNAM,': JAB = ',JAB
                     WRITE(LUPRI,*) SECNAM,
     &                              ': should be between 1 and ',
     &                              NNBSTRT(1)
                     CALL CHO_QUIT(SECNAM//': index error (IRED=1)',
     &                             103)
                  END IF

                  JSHLAB = INDRSH(JAB)
                  IF (JSHLAB .NE. ISP2F(ISHLAB)) THEN
                     WRITE(LUPRI,*) SECNAM,': test is meaningless!'
                     WRITE(LUPRI,*) SECNAM,': JSHLAB must equal ',
     &                              'ISP2F(ISHLAB)'
                     WRITE(LUPRI,*) SECNAM,': JSHLAB,ISP2F(ISHLAB): ',
     &                              JSHLAB,ISP2F(ISHLAB)
                     CALL CHO_QUIT(SECNAM//': shell bug (IRED=1)',
     &                             103)
                  END IF

                  IAB = INDRED(JAB,1)
                  IF ((IAB.LT.1) .OR. (IAB.GT.NUMAB)) THEN
                     WRITE(LUPRI,*) SECNAM,': IAB = ',IAB
                     WRITE(LUPRI,*) SECNAM,
     &                             ': should be between 1 and ',NUMAB
                     CALL CHO_QUIT(SECNAM//': index error (IRED=1)',
     &                             103)
                  END IF

                  KABAB = NUMAB*(IAB - 1) + IAB
                  DIFF  = DIAG(JAB) - xINT(KABAB)
                  IF (ABS(DIFF) .GT. 1.0D-14) THEN
                     WRITE(LUPRI,*)
     &               SECNAM,': ISHLA,ISHLB,JAB,IAB,DIFF: ',
     &               ISHLA,ISHLB,JAB,IAB,DIFF
                     IERR = IERR + 1
                  END IF

                  NTST = NTST + 1

               END DO

            END DO

            WRITE(LUPRI,*)
     &      SECNAM,': ISHLA,ISHLB,#errors: ',ISHLA,ISHLB,IERR

            NERR = NERR + IERR

         ELSE IF ((IRED.EQ.2) .OR. (IRED.EQ.3)) THEN

            DO ISYM = 1,NSYM

               JAB1 = IIBSTR(ISYM,IRED) + IIBSTRSH(ISYM,ISHLAB,IRED)
     &              + 1
               JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,IRED) - 1

               DO JAB = JAB1,JAB2   ! loop over elements in diagonal

                  IF ((JAB.LT.1) .OR. (JAB.GT.NNBSTRT(IRED))) THEN
                     WRITE(LUPRI,*) SECNAM,': JAB = ',JAB
                     WRITE(LUPRI,*) SECNAM,
     &                              ': should be between 1 and ',
     &                              NNBSTRT(IRED)
                     CALL CHO_QUIT(SECNAM//': index error (IRED>1)',
     &                             103)
                  END IF

                  JSHLAB = INDRSH(INDRED(JAB,IRED))
                  IF (JSHLAB .NE. ISP2F(ISHLAB)) THEN
                     WRITE(LUPRI,*) SECNAM,': test is meaningless!'
                     WRITE(LUPRI,*) SECNAM,': JSHLAB must equal ',
     &                              'ISP2F(ISHLAB)'
                     WRITE(LUPRI,*) SECNAM,': JSHLAB,ISP2F(ISHLAB): ',
     &                              JSHLAB,ISP2F(ISHLAB)
                     CALL CHO_QUIT(SECNAM//': shell bug (IRED>1)',
     &                             103)
                  END IF

                  KAB = INDRED(JAB,IRED)  ! index in red. set 1
                  IF ((KAB.LT.1) .OR. (KAB.GT.NNBSTRT(1))) THEN
                     WRITE(LUPRI,*) SECNAM,': KAB = ',KAB
                     WRITE(LUPRI,*) SECNAM,
     &                              ': should be between 1 and ',
     &                              NNBSTRT(1)
                     CALL CHO_QUIT(SECNAM//': index error (IRED>1)',
     &                             103)
                  END IF

                  IAB = INDRED(KAB,1)
                  IF ((IAB.LT.1) .OR. (IAB.GT.NUMAB)) THEN
                     WRITE(LUPRI,*) SECNAM,': IAB = ',IAB
                     WRITE(LUPRI,*) SECNAM,
     &                             ': should be between 1 and ',NUMAB
                     CALL CHO_QUIT(SECNAM//': index error (IRED>1)',
     &                             103)
                  END IF

                  KABAB = NUMAB*(IAB - 1) + IAB
                  DIFF  = DIAG(KAB) - xINT(KABAB)
                  IF (ABS(DIFF) .GT. 1.0D-14) THEN
                     WRITE(LUPRI,*)
     &               SECNAM,': ISHLA,ISHLB,JAB,IAB,DIFF: ',
     &               ISHLA,ISHLB,JAB,IAB,DIFF
                     IERR = IERR + 1
                  END IF

                  NTST = NTST + 1

               END DO

            END DO

            WRITE(LUPRI,*)
     &      SECNAM,': ISHLA,ISHLB,#errors: ',ISHLA,ISHLB,IERR

            NERR = NERR + IERR

         ELSE

            CALL CHO_QUIT(SECNAM//': IRED out of bounds!',104)

         END IF

      END DO

      CALL XRLSMEM_INTS()
      Call mma_deallocate(xINT)

      WRITE(LUPRI,*) '***END OF ',SECNAM,': #tests: ',NTST,
     &               ' #errors: ',NERR

      END
