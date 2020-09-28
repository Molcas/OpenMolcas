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
      SUBROUTINE CHO_SETRED(DIAG,IIBSTRSH,NNBSTRSH,INDRED,
     &                      MSYM,LMMBSTRT,MMSHL)
C
C     Purpose: set next reduced set. A copy of the previous set
C              is stored in location 3.
C
#include "implicit.fh"
      DIMENSION DIAG(*)
      INTEGER   INDRED(LMMBSTRT,3)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_SETRED')

      ISP2F(I)=IWORK(ip_iSP2F-1+I)
      IATOMSHL(I)=IWORK(ip_IATOMSHL-1+I)

#if defined (_DEBUG_)
#endif

C     Debug print.
C     ------------

      IF (CHO_TRCNEG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,
     &                  ': tracing of negative diagonals activated.'
         WRITE(LUPRI,*) SECNAM,': flag SCDIAG     is ',SCDIAG
         WRITE(LUPRI,*) SECNAM,': flag CHO_USEABS is ',CHO_USEABS
         IF (SCDIAG) THEN
            WRITE(LUPRI,*) SECNAM,': MODE_SCREEN     is ',MODE_SCREEN
         END IF
         WRITE(LUPRI,*) SECNAM,': checking for negative diagonals in ',
     &                  'first reduced set:'
         NNEG = 0
         DO ISYM = 1,NSYM
            JAB1 = IIBSTR(ISYM,1) + 1
            JAB2 = JAB1 + NNBSTR(ISYM,1) - 1
            INEG = 0
            DO JAB = JAB1,JAB2
               IF (DIAG(JAB) .LT. 0.0D0) THEN
                  INEG = INEG + 1
               END IF
            END DO
            NNEG = NNEG + INEG
            WRITE(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,
     &                     ': ',INEG
         END DO
         WRITE(LUPRI,*) SECNAM,': total #negative: ',NNEG
      END IF

C     Copy index arrays from location 2 to location 3.
C     ------------------------------------------------

      CALL CHO_RSCOPY(IIBSTRSH,NNBSTRSH,INDRED,2,3,MSYM,MMSHL,LMMBSTRT,
     &                3)

C     Re-initialize index arrays at location 2.
C     -----------------------------------------

      CALL CHO_IZERO(INDRED(1,2),LMMBSTRT)
      CALL CHO_IZERO(IIBSTRSH(1,1,2),MSYM*MMSHL)
      CALL CHO_IZERO(NNBSTRSH(1,1,2),MSYM*MMSHL)
      CALL CHO_IZERO(IIBSTR(1,2),MSYM)
      CALL CHO_IZERO(NNBSTR(1,2),MSYM)
      NNBSTRT(2) = 0

C     Set new reduced set: mapping and SP counter.
C     --------------------------------------------

      IF (SCDIAG) THEN  ! do screening

         IF (MODE_SCREEN .EQ. 1) THEN ! both conv. and unconv. included

            IF (CHO_USEABS) THEN ! neg. diag. might be included

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     IAB = INDRED(JAB1,3)
                     XM  = ABS(DIAG(IAB))
                     DO JAB = JAB1+1,JAB2
                        IAB = INDRED(JAB,3)
                        XM  = MAX(XM,ABS(DIAG(IAB)))
                     END DO

                     IF (XM .GT. THRCOM) THEN  ! only if not converged
                        DO ISHLAB = 1,NNSHL
                           JAB1 = IIBSTR(ISYM,3)
     &                          + IIBSTRSH(ISYM,ISHLAB,3) + 1
                           JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                           DO JAB = JAB1,JAB2
                              IAB = INDRED(JAB,3)
                              TST = SQRT(ABS(DIAG(IAB))*XM)*DAMP(2)
                              IF (TST .GT. THRCOM) THEN
                                 KAB = KAB + 1
                                 INDRED(KAB,2) = IAB
                                 NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                              END IF
                           END DO
                        END DO
                     END IF

                  END IF
               END DO

            ELSE ! neg. diag. excluded

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     JAB1 = IIBSTR(ISYM,3) + 1
                     JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                     IAB = INDRED(JAB1,3)
                     XM  = ABS(DIAG(IAB))
                     DO JAB = JAB1+1,JAB2
                        IAB = INDRED(JAB,3)
                        XM  = MAX(XM,ABS(DIAG(IAB)))
                     END DO

                     IF (XM .GT. THRCOM) THEN  ! only if not converged
                        DO ISHLAB = 1,NNSHL
                           JAB1 = IIBSTR(ISYM,3)
     &                          + IIBSTRSH(ISYM,ISHLAB,3) + 1
                           JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                           DO JAB = JAB1,JAB2
                              IAB = INDRED(JAB,3)
                              IF (DIAG(IAB) .GT. 0.0D0) THEN ! neg=>conv
                                 TST = SQRT(DIAG(IAB)*XM)*DAMP(2)
                                 IF (TST .GT. THRCOM) THEN
                                    KAB = KAB + 1
                                    INDRED(KAB,2) = IAB
                                    NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                                 END IF
                              END IF
                           END DO
                        END DO
                     END IF

                  END IF
               END DO

            END IF

         ELSE IF (MODE_SCREEN .EQ. 2) THEN ! only unconv. included

            IF (CHO_USEABS) THEN ! neg. diag. might be included

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     DO ISHLAB = 1,NNSHL
                        JAB1 = IIBSTR(ISYM,3)
     &                       + IIBSTRSH(ISYM,ISHLAB,3) + 1
                        JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           IF (ABS(DIAG(IAB)) .GT. THRCOM) THEN
                              KAB = KAB + 1
                              INDRED(KAB,2) = IAB
                              NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                           END IF
                        END DO
                     END DO

                  END IF
               END DO

            ELSE ! neg. diag. excluded

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     DO ISHLAB = 1,NNSHL
                        JAB1 = IIBSTR(ISYM,3)
     &                       + IIBSTRSH(ISYM,ISHLAB,3) + 1
                        JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           IF (DIAG(IAB) .GT. THRCOM) THEN
                              KAB = KAB + 1
                              INDRED(KAB,2) = IAB
                              NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                           END IF
                        END DO
                     END DO

                  END IF
               END DO

            END IF

         ELSE IF (MODE_SCREEN .EQ. 3) THEN ! only 1-center unconv. incl.

            IF (CHO_USEABS) THEN ! neg. diag. might be included

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     DO ISHLAB = 1,NNSHL
                        CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,
     &                                  .TRUE.)
                        IF (IATOMSHL(ISHLA) .EQ. IATOMSHL(ISHLB)) THEN
                           JAB1 = IIBSTR(ISYM,3)
     &                          + IIBSTRSH(ISYM,ISHLAB,3) + 1
                           JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                           DO JAB = JAB1,JAB2
                              IAB = INDRED(JAB,3)
                              IF (ABS(DIAG(IAB)) .GT. THRCOM) THEN
                                 KAB = KAB + 1
                                 INDRED(KAB,2) = IAB
                                 NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                              END IF
                           END DO
                        END IF
                     END DO

                  END IF
               END DO

            ELSE ! neg. diag. excluded

               KAB = 0
               DO ISYM = 1,NSYM
                  IF (NNBSTR(ISYM,3) .GT. 0) THEN

                     DO ISHLAB = 1,NNSHL
                        CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,
     &                                  .TRUE.)
                        IF (IATOMSHL(ISHLA) .EQ. IATOMSHL(ISHLB)) THEN
                           JAB1 = IIBSTR(ISYM,3)
     &                          + IIBSTRSH(ISYM,ISHLAB,3) + 1
                           JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                           DO JAB = JAB1,JAB2
                              IAB = INDRED(JAB,3)
                              IF (DIAG(IAB) .GT. THRCOM) THEN
                                 KAB = KAB + 1
                                 INDRED(KAB,2) = IAB
                                 NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                              END IF
                           END DO
                        END IF
                     END DO

                  END IF
               END DO

            END IF

         ELSE ! MODE_SCREEN out of bounds

            WRITE(LUPRI,*) SECNAM,': MODE_SCREEN = ',MODE_SCREEN
            CALL CHO_QUIT('MODE_SCREEN out of bounds in '//SECNAM,103)

         END IF

      ELSE ! no screening; remove zero diagonals and check convergence

         IF (CHO_USEABS) THEN ! neg diag might be incl.

            KAB = 0
            DO ISYM = 1,NSYM
               IF (NNBSTR(ISYM,3) .GT. 0) THEN

                  JAB1 = IIBSTR(ISYM,3) + 1
                  JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                  IAB = INDRED(JAB1,3)
                  XM  = ABS(DIAG(IAB))
                  DO JAB = JAB1+1,JAB2
                     IAB = INDRED(JAB,3)
                     XM  = MAX(XM,ABS(DIAG(IAB)))
                  END DO

                  IF (XM .GT. THRCOM) THEN  ! only if not converged
                     DO ISHLAB = 1,NNSHL
                        JAB1 = IIBSTR(ISYM,3)
     &                       + IIBSTRSH(ISYM,ISHLAB,3) + 1
                        JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           IF (ABS(DIAG(IAB)) .GT. 0.0D0) THEN
                              KAB = KAB + 1
                              INDRED(KAB,2) = IAB
                              NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                           END IF
                        END DO
                     END DO
                  END IF

               END IF
            END DO

         ELSE ! neg diag excl.

            KAB = 0
            DO ISYM = 1,NSYM
               IF (NNBSTR(ISYM,3) .GT. 0) THEN

                  JAB1 = IIBSTR(ISYM,3) + 1
                  JAB2 = JAB1 + NNBSTR(ISYM,3) - 1

                  IAB = INDRED(JAB1,3)
                  XM  = ABS(DIAG(IAB))
                  DO JAB = JAB1+1,JAB2
                     IAB = INDRED(JAB,3)
                     XM  = MAX(XM,ABS(DIAG(IAB)))
                  END DO

                  IF (XM .GT. THRCOM) THEN  ! only if not converged
                     DO ISHLAB = 1,NNSHL
                        JAB1 = IIBSTR(ISYM,3)
     &                       + IIBSTRSH(ISYM,ISHLAB,3) + 1
                        JAB2 = JAB1 + NNBSTRSH(ISYM,ISHLAB,3) - 1
                        DO JAB = JAB1,JAB2
                           IAB = INDRED(JAB,3)
                           IF (DIAG(IAB) .GT. 0.0D0) THEN
                              KAB = KAB + 1
                              INDRED(KAB,2) = IAB
                              NNBSTRSH(ISYM,ISHLAB,2) =
     &                                       NNBSTRSH(ISYM,ISHLAB,2) + 1
                           END IF
                        END DO
                     END DO
                  END IF

               END IF
            END DO

         END IF

      END IF

C     Set remaining index arrays.
C     ---------------------------

      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,MSYM,MMSHL,2)

C     Debug print.
C     ------------

      IF (CHO_TRCNEG) THEN
         WRITE(LUPRI,*) SECNAM,': checking for negative diagonals ',
     &                  'in next reduced set:'
         NNEG = 0
         DO ISYM = 1,NSYM
            JAB1 = IIBSTR(ISYM,2) + 1
            JAB2 = JAB1 + NNBSTR(ISYM,2) - 1
            INEG = 0
            DO JAB = JAB1,JAB2
               IAB = INDRED(JAB,2)
               IF (DIAG(IAB) .LT. 0.0D0) THEN
                  INEG = INEG + 1
               END IF
            END DO
            NNEG = NNEG + INEG
            WRITE(LUPRI,*) SECNAM,': #negative in symmetry ',ISYM,
     &                     ': ',INEG
         END DO
         WRITE(LUPRI,*) SECNAM,': total #negative: ',NNEG
      END IF

#if defined (_DEBUG_)
#endif

      END
