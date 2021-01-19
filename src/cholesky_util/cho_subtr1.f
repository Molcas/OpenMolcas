************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      SUBROUTINE CHO_SUBTR1(XINT,WRK,LWRK,ISYM,FXDMEM)
C
C     Purpose: subtract contributions from previous vectors
C              from the qualified integrals (in XINT).
C              This version is I/O-driven.
C
C     Screening in subtraction introduced Jan. 2006, TBP.
C
      use ChoArr, only: iSP2F, iScr
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh, IndRSh, InfRed,
     &                  InfVec, IndRed
#include "implicit.fh"
      DIMENSION XINT(*), WRK(LWRK)
      LOGICAL   FXDMEM
#include "cholesky.fh"
#include "chovecbuf.fh"
#include "choprint.fh"
#include "chosubscr.fh"
#include "cholq.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_SUBTR1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1), IVSTAT(2,2)

      DIMENSION TIMLOC(2,3)

      PARAMETER (INFO = INF_SUBTR1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER  CHO_X_NUMRD
      EXTERNAL CHO_X_NUMRD

      DSUBSCR(I)=WORK(ip_DSUBSCR-1+I)
      DSPNM(I)=WORK(ip_DSPNM-1+I)

C     Return if nothing to do.
C     ------------------------

      IF (NUMCHO(ISYM) .LT. 1) RETURN

      NVEC_TO_READ = NUMCHO(ISYM) - NVEC_IN_BUF(ISYM)
      IF (NVEC_TO_READ .EQ. 0) RETURN
      IF (NVEC_TO_READ .LT. 0) THEN
         CALL CHO_QUIT('Vector buffer error in '//SECNAM,104)
      END IF

C     Allocate "junk yard".
C     ---------------------

      KJUNK = 1
      KEND0 = KJUNK + 1
      LWRK0 = LWRK  - KEND0 + 1

      MUST = NNBSTR(ISYM,1) + NNBSTR(ISYM,2) + NQUAL(ISYM)
      IF (LWRK0 .LT. MUST) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': insufficient memory:'
         WRITE(LUPRI,*) 'Need at least: ',MUST+KEND0-1
         WRITE(LUPRI,*) 'Available    : ',LWRK
         WRITE(LUPRI,*) '(A significant increase of memory is ',
     &                  'needed for efficient execution.)'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Memory available in ',
     &                  SECNAM,' may also be increased by reducing:'
         WRITE(LUPRI,*)
     &   '1) max. #qualified per symmetry, currently: ',MAXQUAL,
     &   '   to less than than ',NQUAL(ISYM)
         WRITE(LUPRI,*)
     &   '2) max. memory fraction used by qualified, ',
     &   ' currently: ',N1_QUAL,'/',N2_QUAL
         CALL CHO_QUIT('Insufficient memory in '//SECNAM//' [0]',101)
      END IF

      WRK(KJUNK) = ZERO
      IOFF(0)    = KJUNK

C     Split memory for subtraction of previous vectors:
C     {Read buffer},{L(cd,#J),L({ab},#J)}
C     Initially, the fraction N1_VECRD/N2_VECRD of total memory is
C     reserved for reading. Then, the split aims at
C              1<=#J<=MAX(NQUAL(ISYM),N_SUBTR)
C     I.e., the number of vectors in the read buffer is at least #J.
C     Obviously, if MAX(NQUAL(ISYM),N_SUBTR) > NVEC_TO_READ then
C     1<=#J<=NVEC_TO_READ. (Keeping #J within bounds is the role of
C     NUMSMN).
C     N1_VECRD, N2_VECRD, and N_SUBTR can be user-defined (input)
C     and should have been checked as part of configuration check
C     (CHO_CHKCONF) during initialization.
C     --------------------------------------------------------------

      KREAD = KEND0 ! pointer to read buffer
      IREDC = -1    ! id of red. set index array at location 3

      X1 = DBLE(N1_VECRD)
      X2 = DBLE(N2_VECRD)
      XM = DBLE(LWRK0)
      XREAD  = XM*X1/X2  ! initial guess for buffer size (from input)
      LREAD  = INT(XREAD)
      LEFT   = LWRK0 - LREAD
      MINLFT = NNBSTR(ISYM,2) + NQUAL(ISYM)

      IF (FXDMEM) THEN ! fixed-size read buffer

         IF (LEFT .LT. MINLFT) THEN
            LEFT  = MINLFT
            LREAD = LWRK0 - LEFT
         END IF

      ELSE ! attempt to optimize the split

         NUMSMN = MIN(MAX(NQUAL(ISYM),N_SUBTR),NVEC_TO_READ)
         NUMSUB = MAX(MIN(NUMSMN,LEFT/MINLFT),1) ! 1<=NUMSUB<=NUMSMN
         LREAD  = LWRK0 - NUMSUB*MINLFT

         NUMRD = CHO_X_NUMRD(1,ISYM,IREDC,LREAD) ! # that can be read
         IF (NUMRD .LT. 0) CALL CHO_QUIT('NUMRD error in '//SECNAM,104)
         DO WHILE (NUMRD .LT. NUMSUB) ! reduce NUMSUB until NUMRD=NUMSUB
            NUMSUB = NUMSUB - 1
            IF (NUMSUB .LT. 1) THEN ! should never occur (checked above)
               CALL CHO_QUIT('Insufficient memory for split in '
     &                       //SECNAM,101)
            END IF
            LREAD = LWRK0 - NUMSUB*MINLFT
            NUMRD = CHO_X_NUMRD(1,ISYM,IREDC,LREAD)
            IF (NUMRD .LT. 0) THEN
               CALL CHO_QUIT('NUMRD error in '//SECNAM,104)
            END IF
         END DO

      END IF

C     Initializations.
C     ----------------

      NUMRD  = 0
      NUMBAT = 0
      XTOT = 0.0D0
      XDON = 0.0D0
      CALL CHO_IZERO(IVSTAT,4)
      CALL CHO_DZERO(TIMLOC,6)

C     Start buffer batch loop.
C     ------------------------

      IVEC1 = NVEC_IN_BUF(ISYM) + 1
      IMAPC = -1
      DO WHILE (IVEC1 .LE. NUMCHO(ISYM))

C        Read as many vectors as possible into buffer.
C        ---------------------------------------------

         CALL CHO_TIMER(C1,W1)
         NVRD  = 0
         MUSED = 0
         CALL CHO_VECRD(WRK(KREAD),LREAD,IVEC1,NUMCHO(ISYM),ISYM,
     &                  NVRD,IREDC,MUSED)
         NUMRD = NUMRD + 1
         CALL CHO_TIMER(C2,W2)
         TIMLOC(1,1) = TIMLOC(1,1) + C2 - C1
         TIMLOC(2,1) = TIMLOC(2,1) + W2 - W1

C        Quit if no vectors were read.
C        -----------------------------

         IF (NVRD .LT. 1) THEN
            CALL CHO_QUIT('Insufficient scratch space for read in '
     &                    //SECNAM,101)
         END IF

C        Compute memory available for subtraction batching.
C        --------------------------------------------------

         KEND1 = KREAD + MUSED
         LWRK1 = LWRK  - KEND1 + 1

         IF (LWRK1 .LT. 1) THEN
            CALL CHO_QUIT('Insufficient memory in '//SECNAM//' [1]',101)
         END IF

C        Set up batch.
C        -------------

         MMEM = NNBSTR(ISYM,2) + NQUAL(ISYM)
         IF (MMEM .LT. 1) THEN
            CALL CHO_QUIT('Batch setup corrupted in '//SECNAM,104)
            NVEC = -999999
         ELSE
            NVEC = MIN(LWRK1/MMEM,NVRD)
         END IF
         IF (NVEC .LT. 1) THEN
            CALL CHO_QUIT('Batch failure in '//SECNAM,101)
            NBATCH = -999999
         ELSE
            NBATCH = (NVRD - 1)/NVEC + 1
         END IF

C        Set local statistics info.
C        --------------------------

         NUMBAT = NUMBAT + NBATCH
         IF (NUMRD .EQ. 1) THEN
            DO I = 1,2
               IVSTAT(I,1) = NVRD
               IVSTAT(I,2) = NVEC
            END DO
         ELSE
            IVSTAT(1,1) = MIN(IVSTAT(1,1),NVRD)
            IVSTAT(2,1) = MAX(IVSTAT(2,1),NVRD)
            IVSTAT(1,2) = MIN(IVSTAT(1,2),NVEC)
            IVSTAT(2,2) = MAX(IVSTAT(2,2),NVEC)
         END IF

C        Start batch loop.
C        -----------------

         IOFF(1) = KREAD - 1
         DO IBATCH = 1,NBATCH

            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NVRD - NVEC*(NBATCH - 1)
            ELSE
               NUMV = NVEC
            END IF
            IVEC1_1 = IVEC1 + NVEC*(IBATCH-1)

C           Set memory pointers for this batch.
C           -----------------------------------

            KCHO1 = KEND1
            KCHO2 = KCHO1 + NNBSTR(ISYM,2)*NUMV
            KEND2 = KCHO2 + NQUAL(ISYM)*NUMV
            LWRK2 = LWRK  - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               CALL CHO_QUIT('Batch error in '//SECNAM,104)
            END IF

C           Get the next NUMV vectors sorted according to current
C           reduced set (originally, this section was part of the
C           I/O; hence, it is timed as if it was I/O for
C           compatibility).
C           -----------------------------------------------------

            CALL CHO_TIMER(C1,W1)

            JVEC1   = NVEC*(IBATCH - 1) + 1
            JVEC2   = JVEC1 + NUMV - 1
            JRED1   = INFVEC(IVEC1+JVEC1-1,2,ISYM)
            JRED2   = INFVEC(IVEC1+JVEC2-1,2,ISYM)
            LVEC1   = JVEC1
            KVEC1   = 1
            DO JRED = JRED1,JRED2

               LNUM = 0
               LVEC = LVEC1 - 1
               DO WHILE (LVEC.LT.JVEC2)
                  LVEC = LVEC + 1
                  LRED = INFVEC(IVEC1+LVEC-1,2,ISYM)
                  IF (LRED .EQ. JRED) THEN
                     LNUM = LNUM + 1
                  ELSE
                     LVEC = JVEC2
                  END IF
               END DO

               IF (LNUM .GT. 0) THEN

                  IF (JRED .NE. IREDC) THEN
                     ILOC = 3
                     CALL CHO_GETRED(INFRED,nnBstRSh(:,:,ILOC),
     &                               IndRed(:,ILOC),INDRSH,iSP2F,
     &                               MAXRED,NSYM,NNSHL,MMBSTRT,JRED,
     &                               .FALSE.)
                     CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,
     &                                  ILOC)
                     IREDC = JRED
                  END IF

                  IF (JRED .NE. IMAPC) THEN
                     CALL CHO_RS2RS(ISCR,SIZE(ISCR),2,3,JRED,ISYM)
                     IMAPC = JRED
                  END IF

                  DO LVEC = 1,LNUM
                     KVEC = KVEC1 + LVEC - 1
                     DO IAB = 1,NNBSTR(ISYM,2)
                        KOFF1 = KCHO1 + NNBSTR(ISYM,2)*(KVEC - 1) + IAB
     &                        - 1
                        KOFF2 = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                        WRK(KOFF1) = WRK(KOFF2)
                     END DO
                     IOFF(1) = IOFF(1) + NNBSTR(ISYM,3)
                  END DO

                  LVEC1 = LVEC1 + LNUM
                  KVEC1 = KVEC1 + LNUM

               END IF

            END DO

            CALL CHO_TIMER(C2,W2)
            TIMLOC(1,2) = TIMLOC(1,2) + C2 - C1
            TIMLOC(2,2) = TIMLOC(2,2) + W2 - W1

C           Screened or unscreened subtraction section.
C           The screened version uses level 2 blas, while the unscreened
C           one employs level 3 blas.
C           ------------------------------------------------------------

            CALL CHO_TIMER(C1,W1)

            IF (CHO_SSCREEN) THEN ! screened subtraction

C              Copy out sub-blocks corresponding to qualified diagonals:
C              L(#J,{ab})
C              ---------------------------------------------------------

               KOFB0 = KCHO1 - 1 - IIBSTR(ISYM,2)
               DO J = 1,NUMV
                  KOFFA = KCHO2 + J - 1
                  KOFFB = KOFB0 + NNBSTR(ISYM,2)*(J - 1)
                  DO IAB = 1,NQUAL(ISYM)
                     WRK(KOFFA+NUMV*(IAB-1))=WRK(KOFFB+IQUAB(IAB,ISYM))
                  END DO
               END DO

C              Subtract:
C              (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L(#J,{ab})
C              for each ab in {ab}.
C              ----------------------------------------------------

               CALL CHO_SUBSCR_DIA(WRK(KCHO1),NUMV,ISYM,2,SSNORM)
               DO IAB = 1,NQUAL(ISYM)
                  DO ISHGD = 1,NNSHL
                     NGD = NNBSTRSH(ISYM,ISHGD,2)
                     IF (NGD .GT. 0) THEN
                        XTOT = XTOT + 1.0D0
                        JAB = IQUAB(IAB,ISYM) - IIBSTR(ISYM,2)
                        TST = SQRT(DSPNM(ISHGD)*DSUBSCR(JAB))
                        IF (TST .GT. SSTAU) THEN
                           XDON  = XDON + 1.0D0
                           KOFF1 = KCHO1 + IIBSTRSH(ISYM,ISHGD,2)
                           KOFF2 = KCHO2 + NUMV*(IAB-1)
                           KOFF3 = NNBSTR(ISYM,2)*(IAB-1)
     &                           + IIBSTRSH(ISYM,ISHGD,2) + 1
                           CALL DGEMV_('N',NGD,NUMV,
     &                                XMONE,WRK(KOFF1),NNBSTR(ISYM,2),
     &                                WRK(KOFF2),1,ONE,XINT(KOFF3),1)
                        END IF
                     END IF
                  END DO
               END DO

            ELSE ! unscreened subtraction

               IF (L_LQ_SYM(ISYM) .GT. 0) THEN

C                 If the qualified block, L({ab},#J), is already in
C                 core, use this block.
C                 -------------------------------------------------

                  LOFF = IP_LQ_SYM(ISYM) + LDLQ(ISYM)*(IVEC1_1-1)

                  CALL DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,
     &                       XMONE,WRK(KCHO1),NNBSTR(ISYM,2),
     &                             WORK(LOFF),LDLQ(ISYM),
     &                       ONE,XINT,NNBSTR(ISYM,2))

               ELSE

C                 Copy out sub-blocks corresponding to qualified
C                 diagonals: L({ab},#J)
C                 ----------------------------------------------

                  KOFB0 = KCHO1 - 1 - IIBSTR(ISYM,2)
                  DO J = 1,NUMV
                     KOFFA = KCHO2 + NQUAL(ISYM)*(J - 1) - 1
                     KOFFB = KOFB0 + NNBSTR(ISYM,2)*(J - 1)
                     DO IAB = 1,NQUAL(ISYM)
                        WRK(KOFFA+IAB) = WRK(KOFFB+IQUAB(IAB,ISYM))
                     END DO
                  END DO

C                 Subtract:
C                 (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
C                 ----------------------------------------------------

                  CALL DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,
     &                       XMONE,WRK(KCHO1),NNBSTR(ISYM,2),
     &                       WRK(KCHO2),NQUAL(ISYM),
     &                       ONE,XINT,NNBSTR(ISYM,2))

               END IF

            END IF

            CALL CHO_TIMER(C2,W2)
            TIMLOC(1,3) = TIMLOC(1,3) + C2 - C1
            TIMLOC(2,3) = TIMLOC(2,3) + W2 - W1

         END DO

C        Update counter.
C        ---------------

         IVEC1 = IVEC1 + NVRD

      END DO

C     Update global statistics info.
C     ------------------------------

      NSYS_CALL   = NSYS_CALL + NUMRD
      NDGM_CALL   = NDGM_CALL + NUMBAT
      TDECOM(1,2) = TDECOM(1,2) + TIMLOC(1,1) + TIMLOC(1,2)
      TDECOM(2,2) = TDECOM(2,2) + TIMLOC(2,1) + TIMLOC(2,2)
      TDECOM(1,3) = TDECOM(1,3) + TIMLOC(1,3)
      TDECOM(2,3) = TDECOM(2,3) + TIMLOC(2,3)
      IF (CHO_SSCREEN) THEN
         SUBSCRSTAT(1) = SUBSCRSTAT(1) + XTOT
         SUBSCRSTAT(2) = SUBSCRSTAT(2) + XDON
      END IF

C     Print statistics.
C     -----------------

      IF (LOCDBG .OR. IPRINT.GE.INFO) THEN
         IF (NUMRD .EQ. 0) THEN
            XAVERD = -9.99999D5
         ELSE
            XAVERD = DBLE(NUMCHO(ISYM))/DBLE(NUMRD)
         END IF
         IF (NUMBAT .EQ. 0) THEN
            XAVEVC = -9.99999D5
         ELSE
            XAVEVC = DBLE(NUMCHO(ISYM))/DBLE(NUMBAT)
         END IF
         WRITE(LUPRI,'(A)') '*****'
         WRITE(LUPRI,'(A,A,I2,A)')
     &   SECNAM,' statistics, symmetry',ISYM,':'
         WRITE(LUPRI,'(A,I12)')
     &   'Number of previous vectors                           : ',
     &   NUMCHO(ISYM)
         WRITE(LUPRI,'(A,I12)')
     &   'Number of vectors in buffer                          : ',
     &   NVEC_IN_BUF(ISYM)
         WRITE(LUPRI,'(A,I12)')
     &   'Memory available for subtraction of previous vectors : ',
     &   LWRK
         WRITE(LUPRI,'(A,I12)')
     &   'Memory reserved for buffered vector read             : ',
     &   LREAD
         WRITE(LUPRI,'(A,I12)')
     &   'Number of batches needed for reading vectors         : ',
     &   NUMRD
         IF (CHO_SSCREEN) THEN
            WRITE(LUPRI,'(A,F12.2)')
     &      'Number of calls to DGEMV                             : ',
     &      XDON
            WRITE(LUPRI,'(A,1P,D12.2)')
     &      'Screening threshold                                  : ',
     &      SSTAU
            IF (XTOT .GT. 0.0D0) THEN
               SCRPCT = 1.0D2*(XTOT-XDON)/XTOT
            ELSE
               SCRPCT = 1.0D15
            END IF
            WRITE(LUPRI,'(A,F12.2,A)')
     &      'Screening percent                                    : ',
     &      SCRPCT,'%'
         ELSE
            WRITE(LUPRI,'(A,I12)')
     &      'Number of calls to DGEMM                             : ',
     &      NUMBAT
         END IF
         WRITE(LUPRI,'(A,I12,I12,F12.2)')
     &   'Minimum, maximum, and average #vectors read          : ',
     &   IVSTAT(1,1),IVSTAT(2,1),XAVERD
         WRITE(LUPRI,'(A,I12,I12,F12.2)')
     &   'Minimum, maximum, and average #vecs per call to BLAS : ',
     &   IVSTAT(1,2),IVSTAT(2,2),XAVEVC
         WRITE(LUPRI,'(A,2F12.2)')
     &   'Time for reading vectors into buffer (CPU/Wall; sec.): ',
     &   TIMLOC(1,1),TIMLOC(2,1)
         WRITE(LUPRI,'(A,2F12.2)')
     &   'Time for reduced set vector reorder  (CPU/Wall; sec.): ',
     &   TIMLOC(1,2),TIMLOC(2,2)
         WRITE(LUPRI,'(A,2F12.2)')
     &   'Time for qual. copy + subtraction    (CPU/Wall; sec.): ',
     &   TIMLOC(1,3),TIMLOC(2,3)
         WRITE(LUPRI,'(A)') '*****'
      END IF

      END
