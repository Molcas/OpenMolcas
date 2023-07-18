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
      SUBROUTINE CHO_DECOM(DIAG,WRK,LWRK,IPASS,NUM)
!
!     Purpose: calculate Cholesky vectors from qualified integral
!              columns (from disk).
!
      use ChoSwp, only: iQuAB, IndRed
      use ChoVecBuf, only: nVec_in_Buf
      Implicit Real*8 (a-h,o-z)
      Real*8 Diag(*), WRK(LWRK)
#include "cholesky.fh"
#include "choprint.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_DECOM')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      LOGICAL LAST

      INTEGER NUMCHO_OLD(8)

      LENLIN = 0  ! to avoid compiler warnings...
      IF (IPRINT .GE. INF_PROGRESS) THEN
         CALL CHO_HEAD(SECNAM//                                         &
     &                 ': Decomposition of Qualified Diagonals','=',    &
     &                 80,LUPRI)
         WRITE(LUPRI,'(/,A,I5,A,I4,A)')                                 &
     &   'Integral pass number',IPASS,' (',NUM,                         &
     &   ' shell pair distributions calculated)'
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#Cholesky vec.: ',(NUMCHO(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#vec. in buff.: ',(NVEC_IN_BUF(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#qualified    : ',(NQUAL(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   'Current  dim. : ',(NNBSTR(ISYM,2),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   'Original dim. : ',(NNBSTR(ISYM,1),ISYM=1,NSYM)
         WRITE(LUPRI,'(/,A,/,A,A)')                                     &
     &   '           #Vectors             Treated Diagonal',            &
     &   'Sym.     Sym.     Total     Index     Before      After',     &
     &   '   Conv. Neg.   New Max'
         LENLIN = 79
         WRITE(LUPRI,'(80A)') ('-',I=1,LENLIN)
         CALL CHO_FLUSH(LUPRI)
         CALL ICOPY(NSYM,NUMCHO,1,NUMCHO_OLD,1)
      ELSE IF (IPRINT .GE. INF_PASS) THEN
         WRITE(LUPRI,'(/,A,I4)')                                        &
     &   'Number of shell pair distributions calculated:',NUM
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#Cholesky vec.: ',(NUMCHO(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#vec. in buff.: ',(NVEC_IN_BUF(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#qualified    : ',(NQUAL(ISYM),ISYM=1,NSYM)
         CALL CHO_FLUSH(LUPRI)
         CALL ICOPY(NSYM,NUMCHO,1,NUMCHO_OLD,1)
      END IF

!     Decompose each symmetry block.
!     ------------------------------

      DO ISYM = 1,NSYM

!        Cycle loop if nothing to do in this symmetry.
!        ---------------------------------------------

         IF ((NQUAL(ISYM).LT.1) .OR. (NNBSTR(ISYM,2).LT.1)) GO TO 100

!        Reserve space for qualified integral columns.
!        ---------------------------------------------

         LINT1 = NNBSTR(ISYM,2)*NQUAL(ISYM) ! integrals

         KINT1 = 1
         KEND0 = KINT1 + LINT1
         LWRK0 = LWRK  - KEND0 + 1
         IF (LWRK0 .LE. 0) THEN
            CALL CHO_QUIT('[0] Insufficient memory in '//SECNAM,101)
         END IF

!        Determine size of Cholesky vector (output) buffer.
!        --------------------------------------------------

         NUMBUF = MIN(LWRK0/NNBSTR(ISYM,2),NQUAL(ISYM))
         IF (NUMBUF .LT. 1) THEN
            CALL CHO_QUIT('[1] Insufficient memory in '//SECNAM,101)
         END IF
         KCHO1 = KEND0
         KEND1 = KCHO1 + NNBSTR(ISYM,2)*NUMBUF
         LWRK1 = LWRK  - KEND1 + 1
         IF (LWRK1 .LT. 0) THEN  ! should be redundant...
            CALL CHO_QUIT('Buffer allocation error in '//SECNAM,101)
         END IF

!        Read qualified integral columns.
!        --------------------------------

         CALL CHO_TIMER(C1,W1)
         IOPT = 2
         LTOT = NNBSTR(ISYM,2)*NQUAL(ISYM)
         IADR = 0
         CALL DDAFILE(LUSEL(ISYM),IOPT,WRK(KINT1),LTOT,IADR)
         CALL CHO_TIMER(C2,W2)
         TDECOM(1,1) = TDECOM(1,1) + C2 - C1
         TDECOM(2,1) = TDECOM(2,1) + W2 - W1

!        Subtract contributions from previous vectors.
!        ---------------------------------------------

         CALL CHO_SUBTR(WRK(KINT1),WRK(KEND0),LWRK0,ISYM)

!        Debug: check diagonal elements in updated integrals.
!        ----------------------------------------------------

         IF (CHO_DIACHK .OR. LOCDBG) THEN
            TOL  = TOL_DIACHK
            NERR = 0
            CALL CHO_CHKINT(WRK(KINT1),DIAG,ISYM,NERR,TOL,.TRUE.)
            IF (NERR .NE. 0) THEN
               WRITE(LUPRI,*) SECNAM,': ',NERR,' diagonal errors found!'
               WRITE(LUPRI,*) '          #tests: ',NQUAL(ISYM)
!              WRITE(LUPRI,*) '          Printing integrals:'
!              CALL CHO_OUTPUT(WRK(KINT1),
!    &                         1,NNBSTR(ISYM,2),1,NQUAL(ISYM),
!    &                         NNBSTR(ISYM,2),NQUAL(ISYM),1,LUPRI)
               CALL CHO_QUIT('Diagonal errors in '//SECNAM,104)
            ELSE
               WRITE(LUPRI,*) SECNAM,': comparison of qual. integrals ',&
     &                     'and current diagonal: no errors !'
            END IF
         END IF

!        Decompose in loop over qualified columns.
!        -----------------------------------------

         IVEC  = NUMCHO(ISYM)
         IVECT = NUMCHT
         IDUMP = 0
         DO ICHO = 1,NQUAL(ISYM)

!           Find max. diagonal among qualified.
!           -----------------------------------

            IAB  = 1
            IABG = INDRED(IQUAB(IAB,ISYM),2)
            XC   = DIAG(IABG)
            DO I = 2,NQUAL(ISYM)
               KAB = INDRED(IQUAB(I,ISYM),2)
               IF (DIAG(KAB) .GT. XC) THEN
                  IAB  = I
                  IABG = KAB
                  XC   = DIAG(KAB)
               END IF
            END DO

!           Decompose if max. diagonal is still qualified.
!           ----------------------------------------------

            LAST = (XC.LT.DIAMIN(ISYM)) .OR. (XC.LT.THRCOM)
            IF (.NOT. LAST) THEN

!              Offset to max. diagonal column.
!              -------------------------------

               KOFF0 = KINT1 + NNBSTR(ISYM,2)*(IAB - 1) - 1

!              Scale column corresponding to max. diagonal to obtain
!              the Cholesky vector.
!              -----------------------------------------------------

               FAC  = ONE/SQRT(XC)
               KOFF = KOFF0 + 1
               CALL DSCAL_(NNBSTR(ISYM,2),FAC,WRK(KOFF),1)

!              Zero entries in Cholesky vector corresponding to zero
!              diagonals.
!              -----------------------------------------------------

               DO I = 1,NNBSTR(ISYM,2)
                  II = IIBSTR(ISYM,2) + I
                  JJ = INDRED(II,2)
                  IF (DIAG(JJ) .EQ. ZERO) THEN
                     KOFF = KOFF0 + I
                     WRK(KOFF) = ZERO
                  END IF
               END DO

!              Update diagonal.
!              ----------------

               DO I = 1,NNBSTR(ISYM,2)
                  II   = IIBSTR(ISYM,2) + I
                  JJ   = INDRED(II,2)
                  KOFF = KOFF0 + I
                  DIAG(JJ) = DIAG(JJ) - WRK(KOFF)*WRK(KOFF)
               END DO

!              Zero treated diagonal element and analyze updated diagonal.
!              -----------------------------------------------------------

               OLDIAG     = DIAG(IABG)
               DIAG(IABG) = ZERO
               CALL CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XM,NNEGT,NNEG,NCONV)

!              Update total number of zeroed negative diagonals.
!              -------------------------------------------------

               NNZTOT = NNZTOT + NNEG

!              Update DIAMIN from max. abs. diagonal element XM.
!              CHO_1CENTER: update from max. diagonal element among
!                           qualified.
!              CHO_SIMP   : "simulate parallel algorithm" = do not
!                           update DIAMIN.
!              ----------------------------------------------------

               IF (.NOT. CHO_SIMP) THEN
                  IF (CHO_1CENTER) THEN
                     YM = DIAG(INDRED(IQUAB(1,ISYM),2))
                     DO I = 2,NQUAL(ISYM)
                        YM = MAX(YM,DIAG(INDRED(IQUAB(I,ISYM),2)))
                     END DO
                  ELSE
                     YM = XM
                  END IF
                  DIAMIN(ISYM) = MAX(YM*SPAN,THRCOM)
               END IF

!              Subtract this Cholesky vector from integrals. If
!              the corresponding diagonal element is zero, the
!              column will no longer be qualified and subtraction
!              can safely be skipped.
!              --------------------------------------------------

               KOFF1 = KOFF0 + 1
               DO I = 1,NQUAL(ISYM)
                  II = IQUAB(I,ISYM)
                  JJ = INDRED(II,2)
                  IF (DIAG(JJ) .NE. ZERO) THEN
                     KOFF2 = KINT1 + NNBSTR(ISYM,2)*(I - 1)
                     KOFF3 = KOFF0 + II - IIBSTR(ISYM,2)
                     FAC   = -WRK(KOFF3)
                     CALL DAXPY_(NNBSTR(ISYM,2),FAC,WRK(KOFF1),1,       &
     &                                             WRK(KOFF2),1)
                  END IF
               END DO

!              Store Cholesky vector in buffer.
!              --------------------------------

               IDUMP = IDUMP + 1

               KOFF1 = KOFF0 + 1
               KOFF2 = KCHO1 + NNBSTR(ISYM,2)*(IDUMP - 1)
               CALL DCOPY_(NNBSTR(ISYM,2),WRK(KOFF1),1,WRK(KOFF2),1)

!              Update Cholesky vector counters.
!              --------------------------------

               IVEC  = IVEC  + 1
               IVECT = IVECT + 1

!              Set info for this vector.
!              -------------------------

               CALL CHO_SETVECINF(IVEC,ISYM,IABG,IPASS,2)

!              Print progress report.
!              ----------------------

               IF (IPRINT .GE. INF_PROGRESS) THEN
              WRITE(LUPRI,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)')&
     &            ISYM,IVEC,IVECT,IABG,XC,OLDIAG,NCONV,NNEG,XM
               END IF

            END IF

!           Dump vectors to disk when there is no more to be done, or
!           when the buffer is full.
!           ---------------------------------------------------------

            IF (LAST .OR. (IDUMP.EQ.NUMBUF)) THEN
               CALL CHO_TIMER(C1,W1)
               IVEC1 = NUMCHO(ISYM) + 1
               CALL CHO_PUTVEC(WRK(KCHO1),NNBSTR(ISYM,2),IDUMP,IVEC1,   &
     &                         ISYM)
               CALL CHO_VECBUF_COPY(WRK(KCHO1),IDUMP,ISYM)
               NUMCHO(ISYM) = NUMCHO(ISYM) + IDUMP
               NUMCHT       = NUMCHT       + IDUMP
               CALL CHO_TIMER(C2,W2)
               TDECOM(1,2) = TDECOM(1,2) + C2 - C1
               TDECOM(2,2) = TDECOM(2,2) + W2 - W1
               IF (LAST) THEN
                  GO TO 100   ! cycle symmetry loop
               ELSE
                  IVEC  = NUMCHO(ISYM)
                  IVECT = NUMCHT
                  IDUMP = 0
               END IF
            END IF

         END DO

!        Cycle point: go to next symmetry.
!        ---------------------------------

  100    CONTINUE
         IF (IPRINT .GE. INF_PROGRESS) CALL CHO_FLUSH(LUPRI)

      END DO

      IF (IPRINT .GE. INF_PROGRESS) THEN
         DO ISYM = 1,NSYM
            NUMCHO_OLD(ISYM) = NUMCHO(ISYM) - NUMCHO_OLD(ISYM)
         END DO
         WRITE(LUPRI,'(80A)') ('-',I=1,LENLIN)
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#vec. gener.  : ',(NUMCHO_OLD(ISYM),ISYM=1,NSYM)
      ELSE IF (IPRINT .GE. INF_PASS) THEN
         DO ISYM = 1,NSYM
            NUMCHO_OLD(ISYM) = NUMCHO(ISYM) - NUMCHO_OLD(ISYM)
         END DO
         WRITE(LUPRI,'(A,8I8)')                                         &
     &   '#vec. gener.  : ',(NUMCHO_OLD(ISYM),ISYM=1,NSYM)
      END IF

      END
