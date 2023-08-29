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

subroutine CHO_INP(DFONLY,LUNIT,LUOUT)
!
! Purpose: If DFONLY, set defaults only.
!          Else, read and process input for Cholesky decomposition
!          from unit LUNIT. LUOUT is the unit of the output file
!          which is stored internally in the Cholesky program as
!          LUPRI

use Para_Info, only: Is_Real_Par
use Cholesky, only: BLOCKSIZE, CHKONLY, CHO_1CENTER, CHO_ADRVEC, CHO_DECALG, CHO_DECALG_DEF, CHO_DIACHK, CHO_FAKE_PAR, CHO_INTCHK, &
                    CHO_IOVEC, CHO_MINCHK, CHO_NO2CENTER, CHO_PRESCREEN, Cho_Real_Par, CHO_REORD, CHO_SIMP, CHO_SIMRI, &
                    Cho_SScreen, CHO_TRCNEG, CHO_TSTSCREEN, CHO_USEABS, Damp, FRAC_CHVBUF, HALTIT, IALQUA, IFCSEW, IPRINT, LBUF, &
                    LuPri, MaxQual, MaxRed, MaxVec, MinQual, ModRst, MXSHPR, N1_Qual, N1_VecRd, N2_Qual, N2_VecRd, N_Subtr, &
                    NCOL_CHK, n_MySP, RstCho, RstDia, RUN_INTERNAL, RUN_MODE, SCDIAG, Span, SSNorm, SSTau, SubScrStat, &
                    Thr_PreScreen, THR_SIMRI, ThrCom, ThrDef, ThrDiag, Tol_DiaChk, Trace_Idle, ThrNeg, TOONEG, WARNEG, XlDiag
use RICD_Info, only: Do_RI, Thrshld_CD
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: DFONLY
integer(kind=iwp), intent(in) :: LUNIT, LUOUT
integer(kind=iwp), parameter :: LINTSCR = 2, NOPTION = 58
integer(kind=iwp) :: IDKEY, INTSCR(LINTSCR), IRC, MXSHPR_DEF
logical(kind=iwp) :: DECALG_USRDEF, FAKE_USRDEF, FORCEPARALLEL, MXSHPR_USRDEF
character(len=180) :: LINE
character(len=*), parameter :: OPTION(NOPTION) = ['set decomposition threshold                       ', &
                                                  'set print level                                   ', &
                                                  'set buffer length for diagonal calculation        ', &
                                                  'set threshold for discarding initial diag. elem.  ', &
                                                  'set damping factor for first reduced set          ', &
                                                  'set damping factor for later reduced sets         ', &
                                                  'set span                                          ', &
                                                  'set minimum number of qualifieds for decomposition', &
                                                  'set maximum number of shell pair distributions    ', &
                                                  'turn on diag. screening (using damping factors)   ', &
                                                  'turn off diag. screening                          ', &
                                                  'set qualification algorithm                       ', &
                                                  'set threshold for zeroing negative diagonals      ', &
                                                  'set threshold for warning about negative diagonals', &
                                                  'set threshold for shutdown due to negative diag.  ', &
                                                  'check configuration only                          ', &
                                                  'check all integrals (debug)                       ', &
                                                  'diagonal restart                                  ', &
                                                  'decomposition restart with default restart model  ', &
                                                  'set restart model                                 ', &
                                                  'set maximum #qualifieds per symmetry              ', &
                                                  'set maximum #Cholesky vectors per symmetry        ', &
                                                  'set max. #reduced sets (i.e., integral passes)    ', &
                                                  'check specified #columns (debug)                  ', &
                                                  'minimal integral check (debug)                    ', &
                                                  'use abs. value to set up reduced sets             ', &
                                                  'do not use abs. value to set up red. sets         ', &
                                                  'turn on tracing of negative diagonals (debug)     ', &
                                                  'set algorithm for reading Cholesky vectors        ', &
                                                  'set vector reordering (for int. re-generation)    ', &
                                                  'halt execution after decomposition                ', &
                                                  'set memory fraction used as buffer for vec. read  ', &
                                                  'set fraction of memory to be used for qual. col.  ', &
                                                  'set maximum number of vectors in subtraction      ', &
                                                  'set file address mode (WA or DA for vector files) ', &
                                                  'set model used to interface to integral code      ', &
                                                  'one-step (conventional) decomposition             ', &
                                                  'two-step (generate map, then vectors) decomp.     ', &
                                                  'naive decomposition                               ', &
                                                  'set memory fraction used for global vector buffer ', &
                                                  'set diagonal checking and set tolerance (debug)   ', &
                                                  'test vector subtraction screening (statistics)    ', &
                                                  'use screening in vector subtraction               ', &
                                                  'do not use screening in vector subtraction        ', &
                                                  'threshold for screening in vector subtraction     ', &
                                                  'norm to use in vector subtraction screening       ', &
                                                  'decompose 1-center diagonals only                 ', &
                                                  'decompose 1-center diagonals only; excl. 2-center ', &
                                                  'turn off diagonal prescreening                    ', &
                                                  'turn on diagonal prescreening                     ', &
                                                  'set diagonal prescreening threshold               ', &
                                                  'use parallel decomposition algorithm              ', &
                                                  'simulate parallel algorithm (i.e. modified span)  ', &
                                                  'simulate RI (input 1-center threshold)            ', &
                                                  'activate fake parallel decomposition              ', &
                                                  'activate true parallel decomposition (deact. fake)', &
                                                  'set block size for Z vectors (parallel two-step)  ', &
                                                  'activate tracing of idle nodes (parallel run)     '], &
                               SECNAM = 'CHO_INP'
integer(kind=iwp), external :: IPRINTLEVEL
character(len=180), external :: GET_LN

! Define all entries in common blocks.
! Set output unit.
! Set run mode to "internal".
! ------------------------------------

IRC = 0
call CHO_X_SETINC(IRC)
LUPRI = LUOUT
RUN_MODE = RUN_INTERNAL
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_X_SETINC returned error code ',IRC
  write(LUPRI,*) '(most likely due to a programming error...)'
  call CHO_QUIT('Include file initialization error in '//SECNAM,102)
end if
n_MySP = 0
call CHO_SETGLOB()

! Set default parallel configuration.
! CHO_FAKE_PAR = .true. : fake parallel.
! CHO_FAKE_PAR = .false.: true parallel.
! --------------------------------------

CHO_FAKE_PAR = .false.
Cho_Real_Par = Is_Real_Par() .and. (.not. CHO_FAKE_PAR)

! Set default decomposition algorithm (depends on parallel or serial
! installation).
! ------------------------------------------------------------------

call CHO_SETDECALG_DEF()

! Set default MXSHPR (depends on parallel or serial).
! ---------------------------------------------------

call CHO_SETMXSHPR_DEF(MXSHPR_DEF)

! Set defaults.
! -------------

DECALG_USRDEF = .false.
MXSHPR_USRDEF = .false.
FAKE_USRDEF = .false.

CHO_DECALG = CHO_DECALG_DEF ! default decomposition algorithm
CHO_SIMP = .false.          ! do not simulate parallel algorithm
THRCOM = THRDEF             ! default decomposition threshold
LBUF = 1000000              ! buffer length for initial diagonal
THRDIAG = Zero              ! threshold for initial diag. screening
XLDIAG = Zero               ! just an initialization
DAMP(1) = -1.0e9_wp         ! damping for initial screening (neg.=>generic)
DAMP(2) = -1.0e9_wp         ! damping for later screenings (neg.=>generic)
SPAN = 1.0e-2_wp            ! span factor
MINQUAL = 50                ! min. #qual. needed for proceeding
MAXQUAL = 100               ! max. #qual. per symmetry
MXSHPR = MXSHPR_DEF         ! max. #sh. pairs for proceeding (0=>generic)
SCDIAG = .true.             ! screen diagonal during decom.
IALQUA = 2                  ! qualification algorithm
THRNEG = -1.0e-40_wp        ! diag<THRNEG => diag=0
WARNEG = -1.0e-10_wp        ! diag<WARNEG => diag=0, issue warning
TOONEG = -1.0e-8_wp         ! diag<TOONEG => shutdown
CHKONLY = .false.           ! flag for "check input only"
RSTDIA = .false.            ! flag for diagonal restart
RSTCHO = .false.            ! flag for decomposition restart
MODRST = -1                 ! default restart model (if restart at all)
MAXVEC = 0                  ! max. #Cholesky vectors (0=>generic)
MAXRED = 0                  ! max. #reduced sets (0=>generic)
NCOL_CHK = 0                ! #columns to check (0=>all if requested)
CHO_USEABS = .false.        ! do not use abs. value to set up red. sets
CHO_ADRVEC = 1              ! use WA addressing of vector files
CHO_IOVEC = 3               ! i/o model for reading Cholesky vectors
CHO_REORD = .false.         ! reorder Cholesky vectors to full storage
N1_VECRD = 2                ! numerator for fraction of mem. for iovec
N2_VECRD = 3                ! denominator for fraction of mem. for iovec
N1_QUAL = 1                 ! numerator for fraction of mem. for qual.
N2_QUAL = 3                 ! denominator for fraction of mem. for qual.
N_SUBTR = MAXQUAL           ! max. #vectors in subtraction of prev. vecs.
HALTIT = .false.            ! halt execution after decomposition
IFCSEW = 2                  ! get integrals directly in reduced set
FRAC_CHVBUF = 0.35_wp       ! memory fraction used for vector buffer
CHO_SSCREEN = .false.       ! screening in vector subtraction
SSTAU = -1.0e9_wp           ! threshold in vector subtraction (neg=>generic)
SSNORM = 'Max'              ! norm used for vector screening
SUBSCRSTAT(1) = Zero        ! vector screening statistics
SUBSCRSTAT(2) = Zero        ! vector screening statistics
CHO_1CENTER = .false.       ! do not decompose 1-center diagonals only
CHO_NO2CENTER = .false.     ! do not exclude 2-center diagonals
CHO_PRESCREEN = .true.      ! prescreen diagonal
THR_PRESCREEN = -1.0e9_wp   ! diag. prescreen threshold (neg=>generic)
CHO_SIMRI = .false.         ! simulate RI
THR_SIMRI = -1.0e9_wp       ! threshold for qualifying diags. in RI sim
#ifdef _DEBUGPRINT_
CHO_INTCHK = .true.         ! check integrals after decomposition
CHO_MINCHK = .true.         ! minimal integral check
CHO_TRCNEG = .true.         ! tracing of negative diagonals
CHO_DIACHK = .true.         ! check diagonals in qualified columns
TOL_DIACHK = 1.0e-10_wp     ! tolerance for dia. checking (if requested)
CHO_TSTSCREEN = .true.      ! test vector subtraction screening
#else
CHO_INTCHK = .false.        ! check integrals after decomposition
CHO_MINCHK = .false.        ! minimal integral check
CHO_TRCNEG = .false.        ! tracing of negative diagonals
CHO_DIACHK = .false.        ! check diagonals in qualified columns
TOL_DIACHK = 1.0e-14_wp     ! tolerance for dia. checking (if requested)
CHO_TSTSCREEN = .false.     ! test vector subtraction screening
#endif
BLOCKSIZE = 500             ! #vecs in each Z vector block
TRACE_IDLE = .false.        ! trace idle processes

! Set default print level using molcas environment (or whatever was
! set in seward initially).
! -----------------------------------------------------------------

IPRINT = IPRINTLEVEL(-1)
if (IPRINT <= 2) IPRINT = 1
IPRINT = IPRINT-1

! Return if defaults only.
! ------------------------

FORCEPARALLEL = .false.
if (.not. DFONLY) then

  ! Loop through keyword input.
  ! ---------------------------

  do

    ! Read next keyword and get corresponding ID.
    ! -------------------------------------------

    IDKEY = 0
    call CHO_MCA_GETKEY(LUNIT,OPTION,len(OPTION),NOPTION,IDKEY,LUPRI)

    if ((IDKEY >= 1) .and. (IDKEY <= NOPTION)) then  ! key found

      ! Branch for further processing.
      ! ------------------------------

      select case (IDKEY)

        case (1)
          ! Read decomposition threshold.
          ! -----------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,THRCOM)

        case (2)
          ! Read print level.
          ! -----------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,IPRINT)

        case (3)
          ! Read buffer length for diagonal calculation.
          ! --------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,LBUF)

        case (4)
          ! Read threshold for discarding initial diagonal elements.
          ! --------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,THRDIAG)

        case (5)
          ! Read damping factor for 1st reduced set.
          ! ----------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,DAMP(1))

        case (6)
          ! Read damping factor for later reduced set.
          ! ------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,DAMP(2))

        case (7)
          ! Read span factor.
          ! -----------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,SPAN)
          SPAN = min(abs(SPAN),One)

        case (8)
          ! Read minimum number of qualifieds for decomposition.
          ! ----------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,MINQUAL)

        case (9)
          ! Read maximum number of shell pair distributions that should
          ! be calculated before passing to decomposition.
          ! -----------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,MXSHPR)

        case (10)
          ! Turn on screening.
          ! ------------------

          SCDIAG = .true.

        case (11)
          ! Turn off screening.
          ! -------------------

          SCDIAG = .false.

        case (12)
          ! Read algorithm ID for qualification procedure,
          ! 1: qualify from 1 to max.
          ! 2: qualify the largest qualifiables
          ! ----------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,IALQUA)

        case (13)
          ! Read threshold for zeroing (small) negative diagonals.
          ! Diag < THRNEG => Diag = 0.
          ! ------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,THRNEG)

        case (14)
          ! Read threshold for warning about negative diagonals.
          ! Diag < WARNEG => Diag = 0, but issue warning.
          ! ----------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,WARNEG)

        case (15)
          ! Read threshold for shutdown due to negative diagonals.
          ! Diag < TOONEG => shutdown.
          ! ------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,TOONEG)

        case (16)
          ! Check configuration only.
          ! -------------------------

          CHKONLY = .true.

        case (17)
          ! Check all integrals after decomposition.
          ! ----------------------------------------

          CHO_INTCHK = .true.
          CHO_MINCHK = .false.
          NCOL_CHK = 0  ! => all

        case (18)
          ! Restart diagonal.
          ! -----------------

          RSTDIA = .true.

        case (19)
          ! Restart decomposition using default model.
          ! ------------------------------------------

          if (.not. RSTCHO) then
            RSTCHO = .true.
            MODRST = -1  ! use configuration from restart file
          end if

        case (20)
          ! Read restart model and set restart.
          ! -----------------------------------

          RSTCHO = .true.
          LINE = GET_LN(LUNIT)
          call GET_I1(1,MODRST)
          if (MODRST < 0) then
            MODRST = -1  ! use configuration from restart file
          else if (MODRST > 0) then
            MODRST = 1   ! use configuration from input
          end if

        case (21)
          ! Read max. #qualifieds per symmetry.
          ! -----------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,MAXQUAL)

        case (22)
          ! Read max. #Cholesky vectors per symmetry.
          ! -----------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,MAXVEC)

        case (23)
          ! Read max. #reduced sets.
          ! ------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,MAXRED)

        case (24)
          ! Read number of integral columns to check.
          ! -----------------------------------------

          CHO_INTCHK = .true.
          CHO_MINCHK = .false.
          LINE = GET_LN(LUNIT)
          call GET_I1(1,NCOL_CHK)

        case (25)
          ! Minimal integral check.
          ! -----------------------

          CHO_INTCHK = .true.
          CHO_MINCHK = .true.

        case (26)
          ! Use abs. value to set up reduced sets.
          ! --------------------------------------

          CHO_USEABS = .true.

        case (27)
          ! Do not use abs. value to set up reduced sets.
          ! ---------------------------------------------

          CHO_USEABS = .false.

        case (28)
          ! Turn on tracing of negative diagonals.
          ! --------------------------------------

          CHO_TRCNEG = .true.

        case (29)
          ! Read I/O model used for reading Cholesky vectors.
          ! -------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,CHO_IOVEC)

        case (30)
          ! Set vector reordering.
          ! ----------------------

          CHO_REORD = .true.

        case (31)
          ! Halt execution after decomposition.
          ! -----------------------------------

          HALTIT = .true.

        case (32)
          ! Set fraction of memory to be used as vector buffer during subtraction.
          ! ----------------------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I(1,INTSCR,2)
          N1_VECRD = INTSCR(1)
          N2_VECRD = INTSCR(2)

        case (33)
          ! Set fraction of memory to be used for qualified columns.
          ! --------------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I(1,INTSCR,2)
          N1_QUAL = INTSCR(1)
          N2_QUAL = INTSCR(2)

        case (34)
          ! Set max. #vectors in subtraction of prev. vectors.
          ! --------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,N_SUBTR)

        case (35)
          ! Set address mode for vector I/O: WA (=1) or DA (=2).
          ! ----------------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,CHO_ADRVEC)

        case (36)
          ! Set model used to interface to integral code.
          ! ---------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,IFCSEW)
          if (IFCSEW < 1) then
            IFCSEW = 1
          else if (IFCSEW > 2) then
            IFCSEW = 2
          end if

        case (37)
          ! One-step (conventional) decomposition.
          ! --------------------------------------

          CHO_DECALG = 1
          DECALG_USRDEF = .true.

        case (38)
          ! Two-step (generate map, then vectors) decomposition.
          ! ----------------------------------------------------

          CHO_DECALG = 2
          DECALG_USRDEF = .true.

        case (39)
          ! Naive decomposition (remove diagonals as soon as smaller than threshold).
          ! -------------------------------------------------------------------------

          CHO_DECALG = 3
          DECALG_USRDEF = .true.

        case (40)
          ! Set memory fraction used for vector buffer.
          ! -------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,FRAC_CHVBUF)

        case (41)
          ! Turn on diagonal checking of qualified columns both before
          ! and after subtraction of previous vectors. Set tolerance in checking.
          ! ---------------------------------------------------------------------

          CHO_DIACHK = .true.
          LINE = GET_LN(LUNIT)
          call GET_F1(1,TOL_DIACHK)

        case (42)
          ! Test vector subtraction screening.
          ! ----------------------------------

          CHO_TSTSCREEN = .true.

        case (43)
          ! Use screening in vector subtraction.
          ! ------------------------------------

          CHO_SSCREEN = .true.

        case (44)
          ! Do not use screening in vector subtraction.
          ! -------------------------------------------

          CHO_SSCREEN = .false.

        case (45)
          ! Threshold for screening in vector subtraction.
          ! ----------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_F1(1,SSTAU)

        case (46)
          ! Threshold for screening in vector subtraction.
          ! ----------------------------------------------

          LINE = GET_LN(LUNIT)
          call GET_S(1,SSNORM,1)

        case (47)
          ! Decompose 1-center diagonals only. 2-center diagonals are
          ! still included, but never explicitly decomposed.
          ! ---------------------------------------------------------

          CHO_1CENTER = .true.

        case (48)
          ! Decompose 1-center diagonals only. 2-center diagonals are
          ! never included.
          ! ---------------------------------------------------------

          CHO_1CENTER = .true.
          CHO_NO2CENTER = .true.

        case (49)
          ! Turn off diagonal prescreening.
          ! -------------------------------

          CHO_PRESCREEN = .false.

        case (50)
          ! Turn on diagonal prescreening.
          ! ------------------------------

          CHO_PRESCREEN = .true.

        case (51)
          ! Set diagonal prescreening threshold.
          ! ------------------------------------

          CHO_PRESCREEN = .true.
          LINE = GET_LN(LUNIT)
          call GET_F1(1,THR_PRESCREEN)

        case (52)
          ! Use "parallel" decomposition algorithm.
          ! ---------------------------------------

          FORCEPARALLEL = .true.

        case (53)
          ! Use "parallel" decomposition algorithm.
          ! ---------------------------------------

          CHO_SIMP = .true.

        case (54)
          ! Simulate RI.
          ! Activate 1-center decomposition.
          ! Read threshold for removing 1-center diagonals in map generation.
          ! -----------------------------------------------------------------

          CHO_SIMRI = .true.
          CHO_1CENTER = .true.
          LINE = GET_LN(LUNIT)
          call GET_F1(1,THR_SIMRI)

        case (55)
          ! Activate fake parallel:
          ! Run parallel decomposition, then distribute vectors if parallel.
          ! ----------------------------------------------------------------

          CHO_FAKE_PAR = .true.
          FAKE_USRDEF = .true.

        case (56)
          ! Activate true parallel decomposition (deactivate fake).
          ! -------------------------------------------------------

          CHO_FAKE_PAR = .false.
          FAKE_USRDEF = .true.

        case (57)
          ! Read BlockSize (Z vectors).
          ! ---------------------------

          LINE = GET_LN(LUNIT)
          call GET_I1(1,BLOCKSIZE)

        case (58)
          ! Turn on tracing of idle nodes.
          ! ------------------------------

          TRACE_IDLE = .true.

        case default

          ! If this section is executed, it's most likely a bug.
          ! ----------------------------------------------------

          write(LUPRI,*) 'IDKEY = ',IDKEY,' is formally legal, but:'
          write(LUPRI,*) 'Did you forget to change the select case ',SECNAM,'?'
          call CHO_QUIT('Illegal address in select case in '//SECNAM,105)

      end select

    else if (IDKEY == NOPTION+1) then  ! end of input
      exit
    else if (IDKEY == -5) then ! some internal error
      write(LUPRI,*) SECNAM,': internal error detected, IDKEY = ',IDKEY
      IRC = 103
      call CHO_QUIT('Error in '//SECNAM,IRC)
    else if (IDKEY == -1) then ! keyword error
      write(LUPRI,*) SECNAM,': keyword error detected, IDKEY = ',IDKEY
      IRC = 105
      call CHO_QUIT('Error in '//SECNAM,IRC)
    else  ! catch all out-of-bound errors
      write(LUPRI,*) SECNAM,': IDKEY out of bounds: ',IDKEY
      IRC = 104
      call CHO_QUIT('Error in '//SECNAM,IRC)
    end if
  end do
end if

! Post processing.
! ----------------

if (FAKE_USRDEF) then
  ! reset parallel configuration
  Cho_Real_Par = Is_Real_Par() .and. (.not. CHO_FAKE_PAR)
  if (.not. DECALG_USRDEF) then ! reset default decomposition alg
    call CHO_SETDECALG_DEF()
    CHO_DECALG = CHO_DECALG_DEF
  end if
  if (.not. MXSHPR_USRDEF) then ! reset default MXSHPR
    call CHO_SETMXSHPR_DEF(MXSHPR_DEF)
    MXSHPR = MXSHPR_DEF
  end if
end if
call CHO_INP_SETDECALG(FORCEPARALLEL)
if (.not. MXSHPR_USRDEF) then
  if ((CHO_DECALG == 4) .or. (CHO_DECALG == 5) .or. (CHO_DECALG == 6)) then
    MXSHPR = 1
  else
    MXSHPR = 0
  end if
end if

! Put decomposition threshold in RICD_Info.
! -----------------------------------------

if (.not. Do_RI) Thrshld_CD = THRCOM

! Normal exit.
! ------------

return

! "Line" is unused, but the function it comes from (Get_Ln) has side effects
#include "macros.fh"
unused_var(Line)

end subroutine CHO_INP
