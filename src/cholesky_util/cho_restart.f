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
      SUBROUTINE CHO_RESTART(DIAG,WRK,LWRK,DSKDIA,LCONV)
C
C     Purpose: update and analyze diagonal for restart
C              (or check of decomposition). Index arrays
C              for first reduced set must be set up before
C              this routine is called. Reduced set 2, on the
C              other hand, is set up here.
C
      use ChoArr, only: iSP2F
#include "implicit.fh"
      DIMENSION DIAG(*), WRK(LWRK)
      LOGICAL   DSKDIA, LCONV
#include "cholesky.fh"
#include "choptr.fh"
#include "choptr2.fh"
#include "choprint.fh"
#include "chosimri.fh"
#include "WrkSpc.fh"

      external ddot_

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_RESTART')

      INTEGER  CHO_F2SP
      EXTERNAL CHO_F2SP

      LOGICAL SYNC, SCDIAG_SAVE

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0)

      INDRSH(I)=IWORK(ip_INDRSH-1+I)
      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)
      IIBSTRSH(I,J,K)=IWORK(ip_IIBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      IATOMSHL(I)=IWORK(ip_IATOMSHL-1+I)
      MYSP(I)=IWORK(ip_MYSP-1+I)
      ISIMRI(I)=IWORK(ip_ISIMRI-1+I)


C     Read diagonal (in reduced set 1).
C     ---------------------------------

      IF (DSKDIA) THEN
         IOPT = 2
         CALL CHO_IODIAG(DIAG,IOPT)
         CALL CHO_P_SYNCDIAG(DIAG,1)
      END IF

      IF (IPRINT .GE. INF_PASS) THEN
         WRITE(LUPRI,'(/,A,I10,/)')
     &   'Number of diagonal elements (1st reduced set): ',NNBSTRT(1)
      END IF

C     Analyze diagonal before update.
C     -------------------------------

      IF (IPRINT .GE. INF_PASS) THEN
         BIN1 = 1.0D2
         STEP = 1.0D-1
         NBIN = 18
         SYNC = .FALSE.
         CALL CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.TRUE.)
      END IF

C     Copy reduced set 1 to 2.
C     ------------------------

      CALL CHO_RSCOPY(IWORK(ip_IIBSTRSH),IWORK(ip_NNBSTRSH),
     &                IWORK(ip_INDRED),1,2,NSYM,NNSHL,NNBSTRT(1),3)

      IMXAB  = 0
      IMNAB  = 0
      NCONVT = 0
      DO ISYM = 1,NSYM

         NDIM = NNBSTR(ISYM,2)
         NVEC = NUMCHO(ISYM)

         IF (IPRINT .GE. INF_PASS) THEN
            WRITE(LUPRI,'(//,A,I2)')
     &      'Check information, symmetry',ISYM
            WRITE(LUPRI,'(/,A,6X,I12)')
     &      'Dimension, 1st reduced set: ',NDIM
            WRITE(LUPRI,'(A,6X,I12)')
     &      'Number of Cholesky vectors: ',NVEC
         END IF

         IF ((NVEC.GT.0) .AND. (NDIM.GT.0)) THEN

            KDIAG = 1
            KEND1 = KDIAG + NDIM
            LWRK1 = LWRK  - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               CALL CHO_QUIT('Insufficient memory in '//SECNAM,101)
            END IF

C           Save a copy of the original diagonal.
C           -------------------------------------

            CALL DCOPY_(NDIM,DIAG(IIBSTR(ISYM,1)+1),1,WRK(KDIAG),1)

C           Calculate Cholesky diagonal.
C           ----------------------------

            CALL CHO_DIACHO(DIAG,ISYM,WRK(KEND1),LWRK1)

C           Find min. and max. error and save original value.
C           -------------------------------------------------

            ERRMX = -1.0D10
            ERRMN =  1.0D10
            EXAMX = ZERO
            EXAMN = ZERO
            DO JAB = 1,NDIM
               IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
               SAV = WRK(KDIAG-1+IAB-IIBSTR(ISYM,1))
               ERR = ABS(DIAG(IAB))
               IF (ERR .GT. ERRMX) THEN
                  ERRMX = ERR
                  EXAMX = SAV
                  IMXAB = IAB
               END IF
               IF (ERR .LT. ERRMN) THEN
                  ERRMN = ERR
                  EXAMN = SAV
                  IMNAB = IAB
               END IF
            END DO

C           Find min. and max. diagonals, zero too negative diagonals,
C           and screen (if requested).
C           ----------------------------------------------------------

            IF (CHO_DECALG .EQ. 4) THEN
               SCDIAG_SAVE = SCDIAG
               SCDIAG = .FALSE. ! do NOT screen (no zeroing of diags)
               DMX = 1.0d0
               CALL CHO_CHKDIA_A4(DIAG,DMX,ISYM,NNEG,NNEGT,NSCR,XMAX,
     &                            XMIN,XAMAX)
               SCDIAG = SCDIAG_SAVE
            ELSE
               CALL CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XAMAX,NNEGT,NNEG,
     &                         NSCR)
            END IF

C           Count converged diagonals.
C           --------------------------

            NCONV = 0
            DO JAB = 1,NDIM
               IAB = INDRED(IIBSTR(ISYM,2)+JAB,2)
               IF (ABS(DIAG(IAB)) .LE. THRCOM) NCONV = NCONV + 1
            END DO

C           Calculate average and RMS error.
C           --------------------------------

            KOFF   = IIBSTR(ISYM,1) + 1
            XDIM   = DBLE(NDIM)
            RMSERR = SQRT(DDOT_(NDIM,DIAG(KOFF),1,DIAG(KOFF),1)/XDIM)
            AVEERR = CHO_DSUMELM(DIAG(KOFF),NDIM)/XDIM

C           Print.
C           ------

            IF (IPRINT .GE. INF_PASS) THEN
               WRITE(LUPRI,'(A,1P,D18.8)')
     &         'Minimum diagonal          : ',XMIN
               WRITE(LUPRI,'(A,1P,D18.8)')
     &         'Maximum diagonal          : ',XMAX
               WRITE(LUPRI,'(A,1P,D18.8,1X,D18.8)')
     &         'Minimum absolute error    : ',ERRMN,EXAMN
               WRITE(LUPRI,'(A,1P,D18.8,1X,D18.8)')
     &         'Maximum absolute error    : ',ERRMX,EXAMX
               WRITE(LUPRI,'(A,1P,D18.8)')
     &         'Average error             : ',AVEERR
               WRITE(LUPRI,'(A,1P,D18.8)')
     &         'Root-mean-square error    : ',RMSERR
               WRITE(LUPRI,'(A,6X,I12)')
     &         'Converged diagonals       : ',NCONV
               WRITE(LUPRI,'(A,6X,I12)')
     &         'Unconverged diagonals     : ',NDIM-NCONV
               WRITE(LUPRI,'(A,6X,I12)')
     &         'Zeroed negative diagonals : ',NNEG
               IF (CHO_DECALG .NE. 4) THEN ! NSCR is useless here
                  IF (SCDIAG) THEN
                     WRITE(LUPRI,'(A,6X,I12)')
     &               'Screened diagonals        : ',NSCR
                  ELSE
                     WRITE(LUPRI,'(A,6X,I12,A)')
     &               'Screenable diagonals      : ',NSCR,
     &               ' (not screened)'
                  END IF
               END IF
            END IF
            NCONVT = NCONVT + NCONV

         END IF

         CALL CHO_FLUSH(LUPRI)

      END DO

C     Sync diagonal and set reduced set.
C     ----------------------------------

      SYNC = .TRUE.
      CALL CHO_P_SETRED(DIAG,SYNC)
      KRED = XNPASS + 1
      CALL CHO_SETRSDIM(IWORK(ip_NDIMRS),NSYM,MAXRED,KRED,2)

C     Sync and analyze (histogram) updated diagonal.
C     ----------------------------------------------

      CALL CHO_P_SYNCDIAG(DIAG,2)
      IF (IPRINT .GE. INF_PASS) THEN
         BIN1 = 1.0D2
         STEP = 1.0D-1
         NBIN = 18
         SYNC = .FALSE.
         CALL CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.FALSE.)
      END IF

C     Set data for minimal integral checking.
C     ---------------------------------------

      IF (CHO_MINCHK) THEN
         IF (IMXAB .GT. 0) THEN
            ISHLAB = CHO_F2SP(INDRSH(IMXAB))
            IF (ISHLAB .GT. 0) THEN
               CALL CHO_INTCHK_REG('MAX DIAG',ISHLAB,ISHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
         END IF
         IF (IMNAB .GT. 0) THEN
            ISHLAB = CHO_F2SP(INDRSH(IMNAB))
            IF (ISHLAB .GT. 0) THEN
               CALL CHO_INTCHK_REG('MIN DIAG',ISHLAB,ISHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
         END IF
         IF (IMXAB.GT.0 .AND. IMNAB.GT.0) THEN
            ISHLAB = CHO_F2SP(INDRSH(IMXAB))
            JSHLAB = CHO_F2SP(INDRSH(IMNAB))
            IF (ISHLAB.GT.0 .AND. JSHLAB.GT.0) THEN
               CALL CHO_INTCHK_REG('MAX|MIN ',ISHLAB,JSHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
         END IF
         JMXAB = 0
         I0AB  = 0
         XX    = 0.0D0
         DO ISHLAB = 1,NNSHL
            NTOT = 0
            DO ISYM = 1,NSYM
               IF (NNBSTRSH(ISYM,ISHLAB,1) .GT. 0) THEN
                  IAB1 = IIBSTR(ISYM,1) + IIBSTRSH(ISYM,ISHLAB,1) + 1
                  IAB2 = IAB1 + NNBSTRSH(ISYM,ISHLAB,1) - 1
                  DO IAB = IAB1,IAB2
                     IF (DIAG(IAB) .LT. XX) THEN
                        XX    = DIAG(IAB)
                        JMXAB = IAB
                     END IF
                  END DO
                  NTOT = NTOT + 1
               END IF
            END DO
            IF (NTOT .EQ. 0) THEN
               I0AB = ISHLAB
            END IF
         END DO
         IF (JMXAB.GT.0 .AND. JMXAB.NE.IMNAB) THEN
            JSHLAB = CHO_F2SP(INDRSH(JMXAB))
            IF (JSHLAB .GT. 0) THEN
               CALL CHO_INTCHK_REG('NEG DIAG',JSHLAB,JSHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
            IF (IMXAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMXAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MAX|NEG ',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
            IF (IMNAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMNAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MIN|NEG ',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
         END IF
         IF (I0AB .GT. 0) THEN
            JSHLAB = CHO_F2SP(INDRSH(I0AB))
            IF (JSHLAB .GT. 0) THEN
               CALL CHO_INTCHK_REG('EXCL RS1',JSHLAB,JSHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
            IF (IMXAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMXAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MAX|XRS1',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
            IF (IMNAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMNAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MIN|XRS1',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
         END IF
         IF (IABMNZ.GT.0 .AND. IABMNZ.LE.NNSHL) THEN
            JSHLAB = CHO_F2SP(INDRSH(IABMNZ))
            IF (JSHLAB .GT. 0) THEN
               CALL CHO_INTCHK_REG('NEG->ZER',JSHLAB,JSHLAB)
            ELSE
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF
            IF (IMXAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMXAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MAX|NEGZ',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
            IF (IMNAB .GT. 0) THEN
               ISHLAB = CHO_F2SP(INDRSH(IMNAB))
               IF (ISHLAB .GT. 0) THEN
                  CALL CHO_INTCHK_REG('MIN|NEGZ',ISHLAB,JSHLAB)
               ELSE
                  CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
               END IF
            END IF
         END IF
      END IF

C     Set convergence flag.
C     ---------------------

      LCONV = NCONVT.EQ.NNBSTRT(1)

C     For 1-center decomposition, convergence is determined only by
C     judging the 1-center diagonals.
C     For RI simulations, diagonals that were initially zeroed may not
C     be converged, but that is OK as long as the diagonal is smaller
C     than the threshold for deletion, THR_SIMRI.
C     ----------------------------------------------------------------

      IF (CHO_1CENTER .AND. .NOT.LCONV) THEN
#if defined (_DEBUGPRINT_)
         IF (NSYM .NE. 1) THEN
            CALL CHO_QUIT(SECNAM//': CHO_1CENTER on, but NSYM != 1',103)
         END IF
         IF (l_IATOMSHL .LT. NSHELL) THEN
            CALL CHO_QUIT(SECNAM//': iAtomShl not allocated correctly!',
     &                    103)
         END IF
#endif
         LCONV  = .TRUE.
         ISHLAB = 0
         DO WHILE (ISHLAB.LT.NNSHL .AND. LCONV)
            ISHLAB = ISHLAB + 1
            CALL CHO_INVPCK(ISP2F(MYSP(ISHLAB)),ISHLA,ISHLB,.TRUE.)
            IF (IATOMSHL(ISHLA) .EQ. IATOMSHL(ISHLB)) THEN
               IAB1 = IIBSTRSH(1,ISHLAB,1) + 1
               IAB2 = IAB1 + NNBSTRSH(1,ISHLAB,1) - 1
               NCONV = 0
               DO IAB = IAB1,IAB2
                  IF (ABS(DIAG(IAB)) .LE. THRCOM) THEN
                     NCONV = NCONV + 1
                  ELSE IF (CHO_SIMRI) THEN
                     IF (ISIMRI(IAB).EQ.1 .AND.
     &                   ABS(DIAG(IAB)) .LE. THR_SIMRI) THEN
                        NCONV = NCONV + 1
                     END IF
                  END IF
               END DO
               LCONV = NCONV.EQ.NNBSTRSH(1,ISHLAB,1)
            END IF
         END DO
      END IF


      END
