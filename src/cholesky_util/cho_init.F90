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
      SUBROUTINE CHO_INIT(SKIP_PRESCREEN,ALLOCATE_BOOKMARKS)
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iAtomShl,       &
     &                  iShlSO, IntMap
!
!     Purpose: initializations.
!
!              IF (SKIP_PRESCREEN): skip prescreening of diagonal.
!              In this case, NNSHL and array iSP2F must be set
!              externally (the allocation is checked here).
!
!              IF (ALLOCATE_BOOKMARKS): allocate arrays needed to
!              record bookmarks during Cholesky decomposition.
!
      use ChoArr, only: nDimRS, MySP
      use ChoSwp, only: iQuAB_Hidden, iQuAB, nnBstRSh_Hidden, nnBstRSh, &
     &                                       iiBstRSh_Hidden, iiBstRSh, &
     &                                         InfRed_Hidden,   InfRed, &
     &                                         InfVec_Hidden,   InfVec
      use ChoBkm, only: BkmVec, BkmThr, nRow_BkmVec, nCol_BkmVec,       &
     &                   nRow_BkmThr, nCol_BkmThr
      use ChoSubScr, only: Cho_SScreen, SSTau
      use ChoSP, only: nnShl_SP
      use stdalloc, only: mma_allocate
      Implicit None
      LOGICAL SKIP_PRESCREEN
      LOGICAL ALLOCATE_BOOKMARKS
#include "choorb.fh"
#include "cholesky.fh"
#include "choprint.fh"

      Real*8 XXB(8)

      CHARACTER(LEN=1), PARAMETER:: LINE='='
      CHARACTER(LEN=8), PARAMETER:: SECNAM='CHO_INIT'
      CHARACTER(LEN=17), PARAMETER:: STRING='Information from '

      REAL*8, PARAMETER::GBLIM = 2.147483648D9
      Integer :: I, J, MulD2h
      Integer :: IA, IRC, ISHL, ISYM, ISYMA, ISYMB, nBsMax, nConfl,     &
     &           nnBMx, nnBT
      Real*8 :: XA, XB, XXBMx, XXBT

      MULD2H(I,J)=IEOR(I-1,J-1)+1


!     Check settings for parallel runs.
!     Return code: 3 will cause verification to accept this as a passed
!     test (certain options are not available in parallel runs).
!     -----------------------------------------------------------------

      IRC = -1
      CALL CHO_P_CHECK(IRC)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_P_CHECK returned ',irc
!        CALL CHO_QUIT('Error in '//SECNAM,102)
         CALL CHO_QUIT('Parallel option conflicts in '//SECNAM,3)
      END IF

!     Allocate array for tracing idle procs.
!     --------------------------------------

      IF (TRACE_IDLE) THEN
         CALL CHO_TRCIDL_INIT()
      END IF

!     Set diagonal prescreening threshold.
!     ------------------------------------

      IF (SKIP_PRESCREEN) CHO_PRESCREEN = .FALSE.
      IF (CHO_PRESCREEN) THEN
         IF (THR_PRESCREEN .LT. 0.0D0) THEN
            THR_PRESCREEN = MIN(1.0d-14,THRCOM)
         END IF
      END IF

!     Get info from Seward.
!     ---------------------

      CALL CHO_MCA_INIT(SKIP_PRESCREEN)

!     Initialize module ChoSP (enabling use of function CHO_F2SP).
!     ---------------------------------------------------------

      NNSHL_SP = NNSHL

!     Set damping.
!     ------------

      CALL CHO_SETDAMP()

!     Allocate memory for reduced set index arrays.
!     ---------------------------------------------

      Call mma_allocate(iiBstRSh_Hidden,nSym,nnShl,3,                   &
     &                  Label='iiBstRSh_Hidden')
      iiBstRSh => iiBstRSh_Hidden
      Call mma_allocate(nnBstRSh_Hidden,nSym,nnShl,3,                   &
     &                  Label='nnBstRSh_Hidden')
      nnBstRSh => nnBstRSh_Hidden
      Call mma_allocate(IntMap,nnShl,Label='IntMap')
      Call mma_allocate(MySP,nnShl,Label='MySP')

!     Initialize timings etc.
!     -----------------------

      CALL FZERO(TDECDRV,2)
      CALL FZERO(TINTEG,2*NINTEG)
      CALL FZERO(TDECOM,2*NDECOM)
      CALL FZERO(TMISC,2*NMISC)
      CALL IZERO(ICHKQ,4*(NCHKQ+1))
      CALL IZERO(NVECRS1,NSYM)

      DID_DECDRV = .FALSE.

      DIAMNZ = 0.0D0
      IABMNZ = 0
      NNZTOT = 0

      NSYS_CALL = 0
      NDGM_CALL = 0

!     Open files for vector and reduced set storage.
!     Open restart files.
!     ----------------------------------------------

      CALL CHO_UNINI()
      CALL CHO_P_OPENVR(1)

!     Initialize integral SP counter.
!     -------------------------------

      CALL CHO_INIMAP()

!     Allocate memory for INFRED and INFVEC arrays.
!     In so doing, determine the max. #vectors and #reduced sets.
!     -----------------------------------------------------------

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
         Call mma_allocate(InfRed_Hidden,MaxRed,Label='InfRed_Hidden')
         InfRed => InfRed_Hidden
         Call mma_allocate(InfVec_Hidden,MaxVec,INFVEC_N2,nSym,         &
     &                     Label='InfVec_Hidden')
         InfVec => InfVec_Hidden
         Call mma_allocate(nDimRS,NSYM,MAXRED,Label='nDimRS')
      END IF

!     Allocate bookmarks (accuracy and number of Cholesky vectors).
!     Not available with restart.
!     -------------------------------------------------------------

      If (Allocate_Bookmarks) Then
         If (RSTCHO) Then
            nRow_BkmVec=0
            nCol_BkmVec=0
            nRow_BkmThr=0
            nCol_BkmThr=0
         Else
            Call mma_allocate(BkmVec,nSym,MaxRed,Label='BkmVec')
            nRow_BkmVec=nSym
            nCol_BkmVec=0
            Call mma_allocate(BkmThr,nSym,MaxRed,Label='BkmThr')
            nRow_BkmThr=nSym
            nCol_BkmThr=0
         End If
      Else
         nRow_BkmVec=0
         nCol_BkmVec=0
         nRow_BkmThr=0
         nCol_BkmThr=0
      End If

!     Initialize INFRED, INFVEC, vector counter, etc.
!     Special handling depending on Cholesky restart.
!     -----------------------------------------------

      CALL CHO_INIT1()

!     Set threshold for screening in vector subtraction.
!     --------------------------------------------------

      IF (CHO_SSCREEN) THEN
         IF (SSTAU .LT. 0.0D0) THEN
            SSTAU = THRCOM*1.0D-6
         END IF
      END IF

!     Print header and configuration.
!     -------------------------------

      IF (IPRINT .GE. 1) THEN
         CALL CHO_PRTHEAD(.FALSE.)
         CALL CHO_FLUSH(LUPRI)
      END IF

!     Check configuration.
!     --------------------

      NCONFL = 0
      CALL CHO_CHKCONF(NCONFL,.TRUE.)
      IF (CHKONLY) THEN
         WRITE(LUPRI,'(A,A,I4,A)')                                      &
     &   SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
         CALL CHO_QUIT('End of configuration check in '//SECNAM,100)
      ELSE IF (NCONFL .NE. 0) THEN
         WRITE(LUPRI,'(A,A,I4,A)')                                      &
     &   SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
         CALL CHO_QUIT('Configuration conflicts in '//SECNAM,105)
      END IF

!     Allocate and set shell-to-center mapping for 1-center
!     decomposition.
!     -----------------------------------------------------

      IF (CHO_1CENTER) THEN
         Call mma_allocate(iAtomShl,nShell,Label='iAtomShl')

         IRC = -1
         CALL CHO_SETATOMSHL(IRC,IATOMSHL,SIZE(IATOMSHL))
         IF (IRC .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': CHO_SETATOMSHL returned ',IRC
            CALL CHO_QUIT(SECNAM//': shell-to-atom init failed!',102)
         END IF
      END IF

!     Allocate IQUAB array for qualification.
!     Allocate IQUAB_L array for parallel runs.
!     -----------------------------------------

      Call mma_allocate(iQuAB_Hidden,MaxQual,nSym,Label='iQuAB_Hidden')
      iQuAB => iQuAB_Hidden
      CALL CHO_P_INILQ(MAXQUAL,NSYM)

!     Set screening mode.
!     -------------------

      IF (CHO_DECALG.EQ.2 .OR. CHO_DECALG.EQ.3 .OR.                     &
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

!     Print section.
!     --------------

      IF (IPRINT .GE. INF_INIT) THEN

         CALL CHO_HEAD(STRING//SECNAM,LINE,80,LUPRI)

         WRITE(LUPRI,'(/,2X,A,I10)')                                    &
     &   'Number of irreps        : ',NSYM
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Number of SOs           : ',NBAST
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Number of shells        : ',NSHELL
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Number of shell pairs   : ',NNSHL_TOT
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Contributing shell pairs: ',NNSHL
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Max. shell dimension    : ',MXORSH
         WRITE(LUPRI,'(2X,A,I10)')                                      &
     &   'Max. shell pair dim.    : ',MX2SH

         IF (IPRINT .GE. 4) THEN ! debug print

!           Basis size info.
!           ----------------

            WRITE(LUPRI,'(/,2X,A,/,2X,A)')                              &
     &      '  Symmetry        NBAS        IBAS',                       &
     &      '----------------------------------'
            DO ISYM = 1,NSYM
               WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10)')                    &
     &         ISYM,NBAS(ISYM),IBAS(ISYM)
            END DO
            WRITE(LUPRI,'(2X,A)')                                       &
     &      '----------------------------------'

!           Shell info.
!           -----------

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')                       &
     &     '     Shell   Dimension    Symmetry   Dimension      Offset',&
     &     '             (NBSTSH)                (NBASSH)     (IBASSH)',&
     &     '----------------------------------------------------------'
            DO ISHL = 1,NSHELL
               DO ISYM = 1,NSYM
                  IF (ISYM .EQ. 1) THEN
                     WRITE(LUPRI,'(2X,I10,2X,I10,2X,I10,2X,I10,2X,I10)')&
     &               ISHL,NBSTSH(ISHL),                                 &
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  ELSE
                     WRITE(LUPRI,'(26X,I10,2X,I10,2X,I10)')             &
     &               ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
                  END IF
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')                                       &
     &     '----------------------------------------------------------'

            WRITE(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')                       &
     &      '    SO        SO    sym    Shell     Index ',              &
     &      ' (global) (reduced)      (ISOSHL)  (ISHLSO)',              &
     &      '-------------------------------------------'
            DO ISYM = 1,NSYM
               DO I = 1,NBAS(ISYM)
                  IA = IBAS(ISYM) + I
                  WRITE(LUPRI,'(2X,I9,1X,I9,1X,I3,1X,I9,1X,I9)')        &
     &                 IA,I,ISYM,ISOSHL(IA),ISHLSO(IA)
               END DO
            END DO
            WRITE(LUPRI,'(2X,A)')                                       &
     &      '-------------------------------------------'

         END IF

      END IF


      END
