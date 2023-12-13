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

subroutine CHO_INIT(SKIP_PRESCREEN,ALLOCATE_BOOKMARKS)
!
! Purpose: initializations.
!
! IF (SKIP_PRESCREEN): skip prescreening of diagonal.
! In this case, NNSHL and array iSP2F must be set
! externally (the allocation is checked here).
!
! IF (ALLOCATE_BOOKMARKS): allocate arrays needed to
! record bookmarks during Cholesky decomposition.

use Symmetry_Info, only: Mul
use Cholesky, only: BkmThr, BkmVec, CHKONLY, CHO_1CENTER, CHO_DECALG, CHO_NO2CENTER, CHO_PRESCREEN, Cho_SScreen, DIAMNZ, &
                    DID_DECDRV, IABMNZ, iAtomShl, iBas, iBasSh, ICHKQ, iiBstRSh, iiBstRSh_Hidden, INF_INIT, InfRed, InfRed_Hidden, &
                    InfVec, InfVec_Hidden, INFVEC_N2, IntMap, IPRINT, iQuAB, iQuAB_Hidden, iShlSO, iSOShl, LuCho, LuMap, LuPri, &
                    LuRed, LuRst, MaxQual, MaxRed, MaxVec, MODE_SCREEN, MX2SH, MXORSH, MySP, nBas, nBasSh, nBasT, nBstSh, &
                    nCol_BkmThr, nCol_BkmVec, nDGM_call, nDimRS, nnBstRSh, nnBstRSh_Hidden, nnShl, nnShl_SP, nnShl_tot, NNZTOT, &
                    nRow_BkmThr, nRow_BkmVec, nShell, nSym, nSys_call, nVecRS1, RstCho, SSTau, TDECDRV, TDECOM, Thr_PreScreen, &
                    ThrCom, TINTEG, TMISC, Trace_Idle
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: SKIP_PRESCREEN, ALLOCATE_BOOKMARKS
integer(kind=iwp) :: I, IA, IRC, ISHL, ISYM, ISYMA, ISYMB, nBsMax, nConfl, nnBMx, nnBT
real(kind=wp) :: XA, XB, XXB(8), XXBMx, XXBT
real(kind=wp), parameter :: GBLIM = 2.147483648e9_wp
character(len=*), parameter :: LINE = '=', SECNAM = 'CHO_INIT', STRING = 'Information from '

! Check settings for parallel runs.
! Return code: 3 will cause verification to accept this as a passed
! test (certain options are not available in parallel runs).
!     -----------------------------------------------------------------

IRC = -1
call CHO_P_CHECK(IRC)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_P_CHECK returned ',irc
  !call CHO_QUIT('Error in '//SECNAM,102)
  call CHO_QUIT('Parallel option conflicts in '//SECNAM,3)
end if

! Allocate array for tracing idle procs.
! --------------------------------------

if (TRACE_IDLE) call CHO_TRCIDL_INIT()

! Set diagonal prescreening threshold.
! ------------------------------------

if (SKIP_PRESCREEN) CHO_PRESCREEN = .false.
if (CHO_PRESCREEN) then
  if (THR_PRESCREEN < Zero) THR_PRESCREEN = min(1.0e-14_wp,THRCOM)
end if

! Get info from Seward.
! ---------------------

call CHO_MCA_INIT(SKIP_PRESCREEN)

! Initialize nnShl_SP (enabling use of function CHO_F2SP).
! --------------------------------------------------------

NNSHL_SP = NNSHL

! Set damping.
! ------------

call CHO_SETDAMP()

! Allocate memory for reduced set index arrays.
! ---------------------------------------------

call mma_allocate(iiBstRSh_Hidden,nSym,nnShl,3,Label='iiBstRSh_Hidden')
iiBstRSh => iiBstRSh_Hidden
call mma_allocate(nnBstRSh_Hidden,nSym,nnShl,3,Label='nnBstRSh_Hidden')
nnBstRSh => nnBstRSh_Hidden
call mma_allocate(IntMap,nnShl,Label='IntMap')
call mma_allocate(MySP,nnShl,Label='MySP')

! Initialize timings etc.
! -----------------------

TDECDRV(:) = Zero
TINTEG(:,:) = Zero
TDECOM(:,:) = Zero
TMISC(:,:) = Zero
ICHKQ(:,:) = 0
NVECRS1(1:NSYM) = 0

DID_DECDRV = .false.

DIAMNZ = Zero
IABMNZ = 0
NNZTOT = 0

NSYS_CALL = 0
NDGM_CALL = 0

! Open files for vector and reduced set storage.
! Open restart files.
! ----------------------------------------------

LURED = 0
LUCHO(1:NSYM) = 0
LURST = 0
LUMAP = 0
call CHO_P_OPENVR(1)

! Initialize integral SP counter.
! -------------------------------

call CHO_INIMAP()

! Allocate memory for INFRED and INFVEC arrays.
! In so doing, determine the max. #vectors and #reduced sets.
! -----------------------------------------------------------

if ((MAXRED < 1) .or. (MAXVEC < 1)) then
  XXBMX = -1.0e8_wp
  XXBT = Zero
  do ISYM=1,NSYM
    XXB(ISYM) = Zero
    do ISYMB=1,NSYM
      ISYMA = MUL(ISYMB,ISYM)
      if (ISYMA == ISYMB) then
        XA = real(NBAS(ISYMA),kind=wp)
        XXB(ISYM) = XXB(ISYM)+XA*(XA+One)*Half
      else if (ISYMA > ISYMB) then
        XA = real(NBAS(ISYMA),kind=wp)
        XB = real(NBAS(ISYMB),kind=wp)
        XXB(ISYM) = XXB(ISYM)+XA*XB
      end if
    end do
    XXBT = XXBT+XXB(ISYM)     ! total diag. dim.
    XXBMX = max(XXBMX,XXB(ISYM)) ! max. diag. block
  end do
  if (MAXVEC < 1) then
    NBSMAX = NBAS(1)
    do ISYM=2,NSYM
      NBSMAX = max(NBSMAX,NBAS(ISYM))
    end do
    MAXVEC = 20*NBSMAX ! default max. #vectors
    if (XXBMX < GBLIM) then
      NNBMX = int(XXBMX)
      MAXVEC = min(MAXVEC,NNBMX) ! reset if less than default
    end if
  end if
  if (MAXRED < 1) then
    MAXRED = NSYM*MAXVEC ! default max. #red. sets
    if (XXBT < GBLIM) then
      NNBT = int(XXBT)
      MAXRED = min(MAXRED,NNBT) ! reset if less than default
    end if
  end if
end if

if ((MAXRED < 1) .or. (MAXVEC < 1)) then
  write(LUPRI,*) SECNAM,': MAXRED = ',MAXRED
  write(LUPRI,*) SECNAM,': MAXVEC = ',MAXVEC
  call CHO_QUIT('MAXRED/MAXVEC error in '//SECNAM,103)
else
  call mma_allocate(InfRed_Hidden,MaxRed,Label='InfRed_Hidden')
  InfRed => InfRed_Hidden
  call mma_allocate(InfVec_Hidden,MaxVec,INFVEC_N2,nSym,Label='InfVec_Hidden')
  InfVec => InfVec_Hidden
  call mma_allocate(nDimRS,NSYM,MAXRED,Label='nDimRS')
end if

! Allocate bookmarks (accuracy and number of Cholesky vectors).
! Not available with restart.
! -------------------------------------------------------------

if (Allocate_Bookmarks) then
  if (RSTCHO) then
    nRow_BkmVec = 0
    nCol_BkmVec = 0
    nRow_BkmThr = 0
    nCol_BkmThr = 0
  else
    call mma_allocate(BkmVec,nSym,MaxRed,Label='BkmVec')
    nRow_BkmVec = nSym
    nCol_BkmVec = 0
    call mma_allocate(BkmThr,nSym,MaxRed,Label='BkmThr')
    nRow_BkmThr = nSym
    nCol_BkmThr = 0
  end if
else
  nRow_BkmVec = 0
  nCol_BkmVec = 0
  nRow_BkmThr = 0
  nCol_BkmThr = 0
end if

! Initialize INFRED, INFVEC, vector counter, etc.
! Special handling depending on Cholesky restart.
! -----------------------------------------------

call CHO_INIT1()

! Set threshold for screening in vector subtraction.
! --------------------------------------------------

if (CHO_SSCREEN) then
  if (SSTAU < Zero) SSTAU = THRCOM*1.0e-6_wp
end if

! Print header and configuration.
! -------------------------------

if (IPRINT >= 1) then
  call CHO_PRTHEAD(.false.)
  call XFLUSH(LUPRI)
end if

! Check configuration.
! --------------------

NCONFL = 0
call CHO_CHKCONF(NCONFL,.true.)
if (CHKONLY) then
  write(LUPRI,'(A,A,I4,A)') SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
  call CHO_QUIT('End of configuration check in '//SECNAM,100)
else if (NCONFL /= 0) then
  write(LUPRI,'(A,A,I4,A)') SECNAM,':',NCONFL,' conflicts detected in Cholesky config'
  call CHO_QUIT('Configuration conflicts in '//SECNAM,105)
end if

! Allocate and set shell-to-center mapping for 1-center decomposition.
! --------------------------------------------------------------------

if (CHO_1CENTER) then
  call mma_allocate(iAtomShl,nShell,Label='iAtomShl')

  IRC = -1
  call CHO_SETATOMSHL(IRC,IATOMSHL,size(IATOMSHL))
  if (IRC /= 0) then
    write(LUPRI,*) SECNAM,': CHO_SETATOMSHL returned ',IRC
    call CHO_QUIT(SECNAM//': shell-to-atom init failed!',102)
  end if
end if

! Allocate IQUAB array for qualification.
! Allocate IQUAB_L array for parallel runs.
! -----------------------------------------

call mma_allocate(iQuAB_Hidden,MaxQual,nSym,Label='iQuAB_Hidden')
iQuAB => iQuAB_Hidden
call CHO_P_INILQ(MAXQUAL,NSYM)

! Set screening mode.
! -------------------

if ((CHO_DECALG == 2) .or. (CHO_DECALG == 3) .or. (CHO_DECALG == 5) .or. (CHO_DECALG == 6)) then
  if (CHO_1CENTER) then
    if (CHO_NO2CENTER) then ! 2-c removed at diag. calc.
      MODE_SCREEN = 2 ! remove diagonals < THRCOM
    else
      MODE_SCREEN = 3 ! remove 2-c diags and diags < THRCOM
    end if
  else
    MODE_SCREEN = 2 ! remove diagonals < THRCOM
  end if
else
  MODE_SCREEN = 1 ! damped screening
end if

! Print section.
! --------------

if (IPRINT >= INF_INIT) then

  call CHO_HEAD(STRING//SECNAM,LINE,80,LUPRI)

  write(LUPRI,'(/,2X,A,I10)') 'Number of irreps        : ',NSYM
  write(LUPRI,'(2X,A,I10)') 'Number of SOs           : ',NBAST
  write(LUPRI,'(2X,A,I10)') 'Number of shells        : ',NSHELL
  write(LUPRI,'(2X,A,I10)') 'Number of shell pairs   : ',NNSHL_TOT
  write(LUPRI,'(2X,A,I10)') 'Contributing shell pairs: ',NNSHL
  write(LUPRI,'(2X,A,I10)') 'Max. shell dimension    : ',MXORSH
  write(LUPRI,'(2X,A,I10)') 'Max. shell pair dim.    : ',MX2SH

  if (IPRINT >= 4) then ! debug print

    ! Basis size info.
    ! ----------------

    write(LUPRI,'(/,2X,A,/,2X,A)') '  Symmetry        NBAS        IBAS','----------------------------------'
    do ISYM=1,NSYM
      write(LUPRI,'(2X,I10,2X,I10,2X,I10)') ISYM,NBAS(ISYM),IBAS(ISYM)
    end do
    write(LUPRI,'(2X,A)') '----------------------------------'

    ! Shell info.
    ! -----------

    write(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)') '     Shell   Dimension    Symmetry   Dimension      Offset', &
                                          '             (NBSTSH)                (NBASSH)     (IBASSH)', &
                                          '----------------------------------------------------------'
    do ISHL=1,NSHELL
      do ISYM=1,NSYM
        if (ISYM == 1) then
          write(LUPRI,'(2X,I10,2X,I10,2X,I10,2X,I10,2X,I10)') ISHL,NBSTSH(ISHL),ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
        else
          write(LUPRI,'(26X,I10,2X,I10,2X,I10)') ISYM,NBASSH(ISYM,ISHL),IBASSH(ISYM,ISHL)
        end if
      end do
    end do
    write(LUPRI,'(2X,A)') '----------------------------------------------------------'

    write(LUPRI,'(/,2X,A,/,2X,A,/,2X,A)') '    SO        SO    sym    Shell     Index ', &
                                          ' (global) (reduced)      (ISOSHL)  (ISHLSO)', &
                                          '-------------------------------------------'
    do ISYM=1,NSYM
      do I=1,NBAS(ISYM)
        IA = IBAS(ISYM)+I
        write(LUPRI,'(2X,I9,1X,I9,1X,I3,1X,I9,1X,I9)') IA,I,ISYM,ISOSHL(IA),ISHLSO(IA)
      end do
    end do
    write(LUPRI,'(2X,A)') '-------------------------------------------'

  end if

end if

end subroutine CHO_INIT
