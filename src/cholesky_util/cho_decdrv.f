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
      SUBROUTINE CHO_DECDRV(DIAG)
C
C     Purpose: driver for the decomposition of the two-electron integral
C              matrix based on the reduced diagonal.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfRed
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"
#include "choprint.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      INTEGER ISYLST(8)
      REAL*8  DIAMAX_SIMP(8)

      CHARACTER(LEN=10), PARAMETER:: SECNAM = 'CHO_DECDRV'

      CHARACTER*7 FILSEL

      CHARACTER*20 STRING

      LOGICAL CONV, SYNC
      LOGICAL, PARAMETER:: LOCDBG = .FALSE.

      INTEGER, EXTERNAL:: CHO_P_GETMPASS

      Integer, Allocatable:: LSTQSP(:)

C     Start timing.
C     -------------

      CALL CHO_TIMER(TCPU1,TWALL1)

C     Initializations and static settings.
C     IRED=2: points to current reduced set in index arrays.
C     ------------------------------------------------------

      IRED  = 2
      CONV  = .FALSE.

C     Initialize Cholesky vector buffer.
C     ----------------------------------

      CALL CHO_VECBUF_INIT(FRAC_CHVBUF,NNBSTR(1,1))
      IF (LOCDBG .OR. IPRINT.GE.INF_VECBUF) THEN
         CALL CHO_VECBUF_PRINT(LUPRI,NSYM)
      END IF

C     Allocate memory for shell pair based diagonal.
C     It is important that the DIASH allocation is first, as memory is
C     released by flushing back to and including this allocation.
C     ----------------------------------------------------------------

      CALL CHO_P_GETGSP(NGSP)
      CALL GETMEM('DIASH','ALLO','REAL',KDIASH,NGSP)
      CALL GETMEM('ISYSH','ALLO','INTE',KISYSH,NGSP)

C     Set first integral pass.
C     ------------------------

      SYNC = .FALSE.
      NPOTSH = 0
      CALL CHO_P_SETPASS(DIAG,SYNC,WORK(KDIASH),IWORK(KISYSH),IRED,CONV,
     &                   NPOTSH)
      IF (NPOTSH .GT. 0) THEN
         IF (CONV) THEN
            CALL CHO_QUIT('Logical error [0.1] in '//SECNAM,103)
         END IF
      ELSE
         IF (.NOT. CONV) THEN
            CALL CHO_QUIT('Logical error [0.2] in '//SECNAM,103)
         END IF
      END IF

C     Allocate shell pair list.
C     -------------------------

      Call mma_allocate(LSTQSP,MAX(NPOTSH,1),Label='LSTQSP')

C     Loop over integral passes. Continue until convergence or
C     until the max. number of integral passes has been reached.
C     To each integral pass there is associated a reduced set,
C     so the IPASS counter is also used as identifier of reduced
C     set during I/O.
C     ----------------------------------------------------------

      IPASS = XNPASS
      JPASS = 0
      MPASS = CHO_P_GETMPASS(IRED)
      DO WHILE ((.NOT.CONV) .AND. (JPASS.LT.MPASS))

C        Update integral pass counter.
C        -----------------------------

         JPASS = JPASS + 1
         IPASS = IPASS + 1

C        Print.
C        ------

         IF (IPRINT .GE. INF_PASS) THEN
             CALL CHO_TIMER(TLTOT1,WLTOT1)
             WRITE(STRING,'(A13,I7)') 'Integral Pass',IPASS
             CALL CHO_HEAD(STRING,'*',80,LUPRI)
         END IF

C        Update idle proc info.
C        ----------------------

         If (Trace_Idle) Then
            nDim_Now=nnBstR(1,2)
            Do iSym=2,nSym
               nDim_Now=nDim_Now+nnBstR(iSym,2)
            End Do
            Call Cho_TrcIdl_Update(nDim_Now.lt.1)
         End If

C        Debug: print diagonal.
C        ----------------------

         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': debug: diagonal before pass ',
     &                     IPASS
            DO ISYM = 1,NSYM
               ISYLST(ISYM) = ISYM
            END DO
            SYNC = .FALSE.
            CALL CHO_P_PRTDIA(DIAG,SYNC,ISYLST,NSYM,IRED)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) SECNAM,': INFRED before pass ',IPASS
            WRITE(LUPRI,'(10I8)') (INFRED(I),I=1,MIN(IPASS,MAXRED))
         END IF

C        Write index arrays for reduced set to disk
C        and update disk address.
C        ------------------------------------------

         CALL CHO_P_PUTRED(IPASS,IRED)

C        Maintain Cholesky vector buffer.
C        The logicals request that statistics informations are updated
C        in the maintainance routine.
C        -------------------------------------------------------------

         IRC = 0
         IPASS_PREV = IPASS - 1
         CALL CHO_VECBUF_MAINTAIN(IRC,IPASS_PREV,.TRUE.,.TRUE.)
         IF (IRC .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': CHO_VECBUF_MAINTAIN returned ',IRC
            CALL CHO_QUIT('Error detected in '//SECNAM,IRC)
         END IF

C        Open scratch files for qualified integral columns.
C        --------------------------------------------------

         DO ISYM = 1,NSYM
            IF (NNBSTR(ISYM,2) .GT. 0) THEN
               LUSEL(ISYM) = 7
               WRITE(FILSEL,'(A6,I1)') 'CHOSEL',ISYM
               CALL DANAME_WA(LUSEL(ISYM),FILSEL)
            ELSE
               LUSEL(ISYM) = -1
            END IF
         END DO

C        Get integral columns on disk stored in current reduced set.
C        -----------------------------------------------------------

         IF (IPRINT .GE. INF_PASS) CALL CHO_TIMER(TLINT1,WLINT1)
         NUM = 0
         CALL CHO_GETINT(DIAG,WORK(KDIASH),IWORK(KISYSH),
     &                   LSTQSP,NPOTSH,NUM)
         CALL CHO_FLUSH(LUPRI)
         IF (IPRINT .GE. INF_PASS) CALL CHO_TIMER(TLINT2,WLINT2)

C        Decompose the qualified integral columns.
C        -----------------------------------------

         IF (IPRINT .GE. INF_PASS) CALL CHO_TIMER(TLDEC1,WLDEC1)
         IF (CHO_DECALG.EQ.4 .OR. CHO_DECALG.EQ.5 .OR. CHO_DECALG.EQ.6)
     &   THEN
            CALL CHO_DECOM_A4(DIAG,LSTQSP,NUM,IPASS)
         ELSE
            IF (CHO_SIMP) THEN
               CALL CHO_MAXDX(DIAG,DIAMAX_SIMP)
               DO ISYM = 1,NSYM
                  DIAMIN(ISYM) = MAX(THRCOM,DIAMAX_SIMP(ISYM)*SPAN)
               END DO
            END IF
            CALL CHO_MEM('MaxMem','MAX ','REAL',KWRK,LWRK)
            CALL CHO_DECOM(DIAG,WORK(KWRK),LWRK,IPASS,NUM)
            CALL CHO_MEM('MaxMem','FREE','REAL',KWRK,LWRK)
         END IF
         CALL CHO_FLUSH(LUPRI)
         IF (IPRINT .GE. INF_PASS) CALL CHO_TIMER(TLDEC2,WLDEC2)

C        Sync global vector counter.
C        ---------------------------

         CALL CHO_P_SYNCNUMCHO(NUMCHO,NSYM)

C        Write restart info to disk.
C        ---------------------------

         CALL CHO_P_WRRSTC(IPASS)

C        Close scratch files for qualified integral columns.
C        ---------------------------------------------------

         DO ISYM = 1,NSYM
            IF (LUSEL(ISYM) .GT. 0) THEN
               CALL DACLOS(LUSEL(ISYM))
            END IF
         END DO

C        Sync diagonal.
C        --------------

         CALL CHO_P_SYNCDIAG(DIAG,2)

C        Analyze diagonal.
C        -----------------

         IF (IPRINT .GE. INF_PASS) THEN
            BIN1 = 1.0D2
            STEP = 1.0D-1
            NBIN = 18
            SYNC = .FALSE.
            CALL CHO_P_ANADIA(DIAG,SYNC,BIN1,STEP,NBIN,.FALSE.)
         END IF

C        Get next reduced set.
C        ---------------------

         SYNC = .FALSE.
         CALL CHO_P_SETRED(DIAG,SYNC)
         KRED = IPASS + 1
         CALL CHO_SETRSDIM(NDIMRS,NSYM,MAXRED,KRED,IRED)
         IF (IPRINT .GE. INF_PASS) THEN
            CALL CHO_P_PRTRED(2)
            CALL CHO_FLUSH(LUPRI)
         END IF

C        Check convergence and, if not converged, set next integral
C        pass.
C        ----------------------------------------------------------

         SYNC = .FALSE.
         NPOTSH = 0
         CALL CHO_P_SETPASS(DIAG,SYNC,WORK(KDIASH),IWORK(KISYSH),IRED,
     &                      CONV,NPOTSH)
         IF (NPOTSH .GT. 0) THEN
            IF (CONV) THEN
               CALL CHO_QUIT('Logical error [1.1] in '//SECNAM,103)
            END IF
         ELSE
            IF (.NOT. CONV) THEN
               CALL CHO_QUIT('Logical error [1.2] in '//SECNAM,103)
            END IF
         END IF

C        Update bookmarks: store largest diagonal (integral accuracy)
C        and number of Cholesky vectors.
C        ------------------------------------------------------------

         Call Cho_P_UpdateBookmarks(iPass)

C        Print idle report.
C        ------------------

         If (Trace_Idle) Then
            Call Cho_TrcIdl_Report()
         End If

C        Print timing for this pass.
C        ---------------------------

         IF (IPRINT .GE. INF_PASS) THEN
            TLINT = TLINT2 - TLINT1
            WLINT = WLINT2 - WLINT1
            TLDEC = TLDEC2 - TLDEC1
            WLDEC = WLDEC2 - WLDEC1
            CALL CHO_TIMER(TLTOT2,WLTOT2)
            TLTOT = TLTOT2 - TLTOT1
            WLTOT = WLTOT2 - WLTOT1
            WRITE(LUPRI,'(/,A,I7,A)')
     &      'Overall timings for integral pass',IPASS,
     &      ' (CPU/Wall in seconds):'
            WRITE(LUPRI,'(A,F12.2,1X,F12.2)')
     &      'Integrals (incl. qualified I/O etc.): ',TLINT,WLINT
            WRITE(LUPRI,'(A,F12.2,1X,F12.2)')
     &      'Decomposition of qualified columns  : ',TLDEC,WLDEC
            WRITE(LUPRI,'(A,F12.2,1X,F12.2)')
     &      'Total (incl. restart info I/O etc.) : ',TLTOT,WLTOT
         END IF

      END DO

C     Free memory for shell pair based diagonal.
C     ------------------------------------------

      Call mma_deallocate(LSTQSP)
      CALL GETMEM('ISYSH','FREE','INTE',KISYSH,NGSP)
      CALL GETMEM('DIASH','FREE','REAL',KDIASH,NNSHL)

C     Shut down the Cholesky vector buffer.
C     -------------------------------------

      CALL CHO_VECBUF_FINAL()

C     Set stuff for statistics.
C     -------------------------

      DID_DECDRV = .TRUE.
      XNPASS     = IPASS

C     Timing.
C     -------

      CALL CHO_TIMER(TCPU2,TWALL2)
      TDECDRV(1) = TCPU2  - TCPU1
      TDECDRV(2) = TWALL2 - TWALL1


      END
