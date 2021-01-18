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
      SUBROUTINE CHO_INIT(SKIP_PRESCREEN,ALLOCATE_BOOKMARKS)
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iAtomShl,
     &                  iShlSO, IntMap
C
C     Purpose: initializations.
C
C              IF (SKIP_PRESCREEN): skip prescreening of diagonal.
C              In this case, NNSHL and array iSP2F must be set
C              externally (the allocation is checked here).
C
C              IF (ALLOCATE_BOOKMARKS): allocate arrays needed to
C              record bookmarks during Cholesky decomposition.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: iQuAB_Hidden, iQuAB, nnBstRSh_Hidden, nnBstRSh,
     &                                       iiBstRSh_Hidden, iiBstRSh,
     &                                         InfRed_Hidden,   InfRed
#include "implicit.fh"
      LOGICAL SKIP_PRESCREEN
      LOGICAL ALLOCATE_BOOKMARKS
#include "choorb.fh"
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "choptr2.fh"
#include "chosp.fh"
#include "chosubscr.fh"
#include "chobkm.fh"
#include "stdalloc.fh"

      DIMENSION XXB(8)

      CHARACTER*1  LINE
      CHARACTER*8  SECNAM
      CHARACTER*17 STRING
      PARAMETER (LINE = '=', SECNAM = 'CHO_INIT',
     &           STRING = 'Information from ')

      PARAMETER (GBLIM = 2.147483648D9)

      MULD2H(I,J)=IEOR(I-1,J-1)+1


C     Check settings for parallel runs.
C     Return code: 3 will cause verification to accept this as a passed
C     test (certain options are not available in parallel runs).
C     -----------------------------------------------------------------

      IRC = -1
      CALL CHO_P_CHECK(IRC)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_P_CHECK returned ',irc
C        CALL CHO_QUIT('Error in '//SECNAM,102)
         CALL CHO_QUIT('Parallel option conflicts in '//SECNAM,3)
      END IF

C     Allocate array for tracing idle procs.
C     --------------------------------------

      IF (TRACE_IDLE) THEN
         CALL CHO_TRCIDL_INIT()
      END IF

C     Set diagonal prescreening threshold.
C     ------------------------------------

      IF (SKIP_PRESCREEN) CHO_PRESCREEN = .FALSE.
      IF (CHO_PRESCREEN) THEN
         IF (THR_PRESCREEN .LT. 0.0D0) THEN
            THR_PRESCREEN = MIN(1.0d-14,THRCOM)
         END IF
      END IF

C     Get info from Seward.
C     ---------------------

      CALL CHO_MCA_INIT(SKIP_PRESCREEN)

C     Initialize chosp.fh (enabling use of function CHO_F2SP).
C     ---------------------------------------------------------

      NNSHL_SP = NNSHL

C     Set damping.
C     ------------

      CALL CHO_SETDAMP()

C     Allocate memory for reduced set index arrays.
C     ---------------------------------------------

      l_MYSP     = NNSHL
      Call mma_allocate(iiBstRSh_Hidden,nSym,nnShl,3,
     &                  Label='iiBstRSh_Hidden')
      iiBstRSh => iiBstRSh_Hidden
      Call mma_allocate(nnBstRSh_Hidden,nSym,nnShl,3,
     &                  Label='nnBstRSh_Hidden')
      nnBstRSh => nnBstRSh_Hidden
      Call mma_allocate(IntMap,nnShl,Label='IntMap')
      CALL CHO_MEM('mySP','ALLO','INTE',ip_MYSP,l_MYSP)

C     Initialize timings etc.
C     -----------------------

      CALL CHO_DZERO(TDECDRV,2)
      CALL CHO_DZERO(TINTEG,2*NINTEG)
      CALL CHO_DZERO(TDECOM,2*NDECOM)
      CALL CHO_DZERO(TMISC,2*NMISC)
      CALL CHO_IZERO(ICHKQ,4*(NCHKQ+1))
      CALL CHO_IZERO(NVECRS1,NSYM)

      DID_DECDRV = .FALSE.

      DIAMNZ = 0.0D0
      IABMNZ = 0
      NNZTOT = 0

      NSYS_CALL = 0
      NDGM_CALL = 0

C     Open files for vector and reduced set storage.
C     Open restart files.
C     ----------------------------------------------

      CALL CHO_UNINI()
      CALL CHO_P_OPENVR(1)

C     Initialize integral SP counter.
C     -------------------------------

      CALL CHO_INIMAP()

C     Allocate memory for INFRED and INFVEC arrays.
C     In so doing, determine the max. #vectors and #reduced sets.
C     -----------------------------------------------------------

      IF (MAXRED.LT.1 .OR. MAXVEC.LT.1) THEN
         XXBMX = -1.0D8
         XXBT  = 0.0D0
         DO ISYM = 1,NSYM
            XXB(ISYM) = 0.0D0
            DO ISYMB = 1,NSYM
               ISYMA = MULD2H(ISYMB,ISYM)
               IF (ISYMA .EQ. ISYMB) THEN
                  XA = DBLE(NBAS(ISYMA))
                  XXB(ISYM) = XXB(ISYM) + XA*(XA + 1.0D0)/2.0D0
               ELSE IF (ISYMA .GT. ISYMB) THEN
                  XA = DBLE(NBAS(ISYMA))
                  XB = DBLE(NBAS(ISYMB))
                  XXB(ISYM) = XXB(ISYM) + XA*XB
               END IF
            END DO
            XXBT  = XXBT + XXB(ISYM)     ! total diag. dim.
            XXBMX = MAX(XXBMX,XXB(ISYM)) ! max. diag. block
         END DO
         IF (MAXVEC .LT. 1) THEN
            NBSMAX = NBAS(1)
            DO ISYM = 2,NSYM
               NBSMAX = MAX(NBSMAX,NBAS(ISYM))
            END DO
            MAXVEC = 20*NBSMAX ! default max. #vectors
            IF (XXBMX .LT. GBLIM) THEN
               NNBMX  = INT(XXBMX)
               MAXVEC = MIN(MAXVEC,NNBMX) ! reset if less than default
            END IF
         END IF
         IF (MAXRED .LT. 1) THEN
            MAXRED = NSYM*MAXVEC ! default max. #red. sets
            IF (XXBT .LT. GBLIM) THEN
               NNBT   = INT(XXBT)
               MAXRED = MIN(MAXRED,NNBT) ! reset if less than default
            END IF
         END IF
      END IF

      IF (MAXRED.LT.1 .OR. MAXVEC.LT.1) THEN
         WRITE(LUPRI,*) SECNAM,': MAXRED = ',MAXRED
         WRITE(LUPRI,*) SECNAM,': MAXVEC = ',MAXVEC
         CALL CHO_QUIT('MAXRED/MAXVEC error in '//SECNAM,103)
      ELSE
         l_INFVEC = MAXVEC*INFVEC_N2*NSYM
         Call mma_allocate(InfRed_Hidden,MaxRed,Label='InfRed_Hidden')
         InfRed => InfRed_Hidden
         CALL CHO_MEM('INFVEC','ALLO','INTE',ip_INFVEC,l_INFVEC)
         Call mma_allocate(nDimRS,NSYM,MAXRED,Label='nDimRS')
      END IF

C     Allocate bookmarks (accuracy and number of Cholesky vectors).
C     Not available with restart.
C     -------------------------------------------------------------

      If (Allocate_Bookmarks) Then
         If (RSTCHO) Then
            ip_BkmVec=0
            l_BkmVec=0
            nRow_BkmVec=0
            nCol_BkmVec=0
            ip_BkmThr=0
            l_BkmThr=0
            nRow_BkmThr=0
            nCol_BkmThr=0
         Else
            l_BkmVec=nSym*MaxRed
            Call GetMem('BkmVec','Allo','Inte',ip_BkmVec,l_BkmVec)
            nRow_BkmVec=nSym
            nCol_BkmVec=0
            l_BkmThr=nSym*MaxRed
            Call GetMem('BkmThr','Allo','Real',ip_BkmThr,l_BkmThr)
            nRow_BkmThr=nSym
            nCol_BkmThr=0
         End If
      Else
         ip_BkmVec=0
         l_BkmVec=0
         nRow_BkmVec=0
         nCol_BkmVec=0
         ip_BkmThr=0
         l_BkmThr=0
         nRow_BkmThr=0
         nCol_BkmThr=0
      End If

C     Initialize INFRED, INFVEC, vector counter, etc.
C     Special handling depending on Cholesky restart.
C     -----------------------------------------------

      CALL CHO_INIT1()

C     Set threshold for screening in vector subtraction.
C     --------------------------------------------------

      IF (CHO_SSCREEN) THEN
         IF (SSTAU .LT. 0.0D0) THEN
            SSTAU = THRCOM*1.0D-6
         END IF
      END IF

C     Print header and configuration.
C     -------------------------------

      IF (IPRINT .GE. 1) THEN
         CALL CHO_PRTHEAD(.FALSE.)
         CALL CHO_FLUSH(LUPRI)
      END IF

C     Check configuration.
C     --------------------

      NCONFL = 0
      CALL CHO_CHKCONF(NCONFL,.TRUE.)
      IF (CHKONLY) THEN
         WRITE(LUPRI,'(A,A,I4,A)')
     &   SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
         CALL CHO_QUIT('End of configuration check in '//SECNAM,100)
      ELSE IF (NCONFL .NE. 0) THEN
         WRITE(LUPRI,'(A,A,I4,A)')
     &   SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
         CALL CHO_QUIT('Configuration conflicts in '//SECNAM,105)
      END IF

C     Allocate and set shell-to-center mapping for 1-center
C     decomposition.
C     -----------------------------------------------------

      IF (CHO_1CENTER) THEN
         Call mma_allocate(iAtomShl,nShell,Label='iAtomShl')

         IRC = -1
         CALL CHO_SETATOMSHL(IRC,IATOMSHL,SIZE(IATOMSHL))
         IF (IRC .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': CHO_SETATOMSHL returned ',IRC
            CALL CHO_QUIT(SECNAM//': shell-to-atom init failed!',102)
         END IF
      END IF

C     Allocate IQUAB array for qualification.
C     Allocate IQUAB_L array for parallel runs.
C     -----------------------------------------

      Call mma_allocate(iQuAB_Hidden,MaxQual,nSym,Label='iQuAB_Hidden')
      iQuAB => iQuAB_Hidden
      CALL CHO_P_INILQ(MAXQUAL,NSYM)

C     Set screening mode.
C     -------------------

      IF (CHO_DECALG.EQ.2 .OR. CHO_DECALG.EQ.3 .OR.
     &    CHO_DECALG.EQ.5 .OR. CHO_DECALG.EQ.6) THEN
         IF (CHO_1CENTER) THEN
            IF (CHO_NO2CENTER) THEN ! 2-c removed at diag. calc.
               MODE_SCREEN = 2 ! remove diagonals < THRCOM
            ELSE
               MODE_SCREEN = 3 ! remove 2-c diags and diags < THRCOM
            END IF
         ELSE
            MODE_SCREEN = 2 ! remove diagonals < THRCOM
         END IF
      ELSE
         MODE_SCREEN = 1 ! damped screening
      END IF

C     Print section.
C     --------------

      IF (IPRINT .GE. INF_INIT) THEN

         CALL CHO_HEAD(STRING//SECNAM,LINE,80,LUPRI)

         WRITE(LUPRI,'(/,2X,A,I10)')
     &   'Number of irreps        : ',NSYM
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of SOs           : ',NBAST
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of shells        : ',NSHELL
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Number of shell pairs   : ',NNSHL_TOT
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Contributing shell pairs: ',NNSHL
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Max. shell dimension    : ',MXORSH
         WRITE(LUPRI,'(2X,A,I10)')
     &   'Max. shell pair dim.    : ',MX2SH

         IF (IPRINT .GE. 4) THEN ! debug print

C           Basis size info.
C           ----------------

            WRITE(LUPRI,'(/,2X,A,/,2X,A)')
     &      '  Symmetry        NBAS        IBAS',
     &      '----------------------------------'
            DO ISYM = 1,NSYM
               WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10)')
     &         ISYM,NBAS(ISYM),IBAS(ISYM)
            END DO
            WRITE(LUPRI,'(2X,A)')
     &      '----------------------------------'

C           Shell info.
C           -----------

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')
     &     '     Shell   Dimension    Symmetry   Dimension      Offset',
     &     '             (NBSTSH)                (NBASSH)     (IBASSH)',
     &     '----------------------------------------------------------'
            DO ISHL = 1,NSHELL
               DO ISYM = 1,NSYM
                  IF (ISYM .EQ. 1) THEN
                     WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10,2X,I10,2X,I10)')
     &               ISHL,NBSTSH(ISHL),
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  ELSE
                     WRITE(LUPRI,'(26X,I10,2X,I10,2X,I10)')
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  END IF
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')
     &     '----------------------------------------------------------'

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')
     &      '    SO        SO    sym    Shell     Index ',
     &      ' (global) (reduced)      (ISOSHL)  (ISHLSO)',
     &      '-------------------------------------------'
            DO ISYM = 1,NSYM
               DO I = 1,NBAS(ISYM)
                  IA = IBAS(ISYM) + I
                  WRITE(LUPRI,'(2X,I9,1X,I9,1X,I3,1X,I9,1X,I9)')
     &                 IA,I,ISYM,ISOSHL(IA),ISHLSO(IA)
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')
     &      '-------------------------------------------'

         END IF

      END IF


      END
