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
!          LUPRI (in cholesky.fh).

use ChoSubScr, only: Cho_SScreen, SSTau, SubScrStat, SSNorm

implicit real*8(a-h,o-z)
logical DFONLY
#include "cholesky.fh"
#include "choprint.fh"
character*7 SECNAM
parameter(SECNAM='CHO_INP')
logical DECALG_USRDEF, MXSHPR_USRDEF, FAKE_USRDEF
logical FORCEPARALLEL
character*180 LINE
character*180 GET_LN
external GET_LN
integer IPRINTLEVEL
external IPRINTLEVEL
parameter(LINTSCR=2)
integer INTSCR(LINTSCR)
parameter(NOPTION=58,LOPTION=50)  ! #options
character*(LOPTION) OPTION(NOPTION)

data OPTION/'set decomposition threshold                       ','set print level                                   ', &
            'set buffer length for diagonal calculation        ','set threshold for discarding initial diag. elem.  ', &
            'set damping factor for first reduced set          ','set damping factor for later reduced sets         ', &
            'set span                                          ','set minimum number of qualifieds for decomposition', &
            'set maximum number of shell pair distributions    ','turn on diag. screening (using damping factors)   ', &
            'turn off diag. screening                          ','set qualification algorithm                       ', &
            'set threshold for zeroing negative diagonals      ','set threshold for warning about negative diagonals', &
            'set threshold for shutdown due to negative diag.  ','check configuration only                          ', &
            'check all integrals (debug)                       ','diagonal restart                                  ', &
            'decomposition restart with default restart model  ','set restart model                                 ', &
            'set maximum #qualifieds per symmetry              ','set maximum #Cholesky vectors per symmetry        ', &
            'set max. #reduced sets (i.e., integral passes)    ','check specified #columns (debug)                  ', &
            'minimal integral check (debug)                    ','use abs. value to set up reduced sets             ', &
            'do not use abs. value to set up red. sets         ','turn on tracing of negative diagonals (debug)     ', &
            'set algorithm for reading Cholesky vectors        ','set vector reordering (for int. re-generation)    ', &
            'halt execution after decomposition                ','set memory fraction used as buffer for vec. read  ', &
            'set fraction of memory to be used for qual. col.  ','set maximum number of vectors in subtraction      ', &
            'set file address mode (WA or DA for vector files) ','set model used to interface to integral code      ', &
            'one-step (conventional) decomposition             ','two-step (generate map, then vectors) decomp.     ', &
            'naive decomposition                               ','set memory fraction used for global vector buffer ', &
            'set diagonal checking and set tolerance (debug)   ','test vector subtraction screening (statistics)    ', &
            'use screening in vector subtraction               ','do not use screening in vector subtraction        ', &
            'threshold for screening in vector subtraction     ','norm to use in vector subtraction screening       ', &
            'decompose 1-center diagonals only                 ','decompose 1-center diagonals only; excl. 2-center ', &
            'turn off diagonal prescreening                    ','turn on diagonal prescreening                     ', &
            'set diagonal prescreening threshold               ','use parallel decomposition algorithm              ', &
            'simulate parallel algorithm (i.e. modified span)  ','simulate RI (input 1-center threshold)            ', &
            'activate fake parallel decomposition              ','activate true parallel decomposition (deact. fake)', &
            'set block size for Z vectors (parallel two-step)  ','activate tracing of idle nodes (parallel run)     '/

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
call CHO_SETPTR2()
call CHO_SETGLOB()

! Set default parallel configuration.
! CHO_FAKE_PAR = .TRUE. : fake parallel.
! CHO_FAKE_PAR = .FALSE.: true parallel.
! --------------------------------------

CHO_FAKE_PAR = .false.
call CHO_PARCONF(CHO_FAKE_PAR)

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
CHO_SIMP = .false. ! do not simulate parallel algorithm
THRCOM = THRDEF   ! default decomposition threshold
LBUF = 1000000  ! buffer length for initial diagonal
THRDIAG = 0.0d0    ! threshold for initial diag. screening
XLDIAG = 0.0d0    ! just an initialization
DAMP(1) = -1.0d9   ! damping for initial screening (neg.=>generic)
DAMP(2) = -1.0d9   ! damping for later screenings (neg.=>generic)
SPAN = 1.0D-2   ! span factor
MINQUAL = 50       ! min. #qual. needed for proceeding
MAXQUAL = 100      ! max. #qual. per symmetry
MXSHPR = MXSHPR_DEF ! max. #sh. pairs for proceeding (0=>generic)
SCDIAG = .true.   ! screen diagonal during decom.
IALQUA = 2        ! qualification algorithm
THRNEG = -1.0D-40 ! diag<THRNEG => diag=0
WARNEG = -1.0D-10 ! diag<WARNEG => diag=0, issue warning
TOONEG = -1.0D-8  ! diag<TOONEG => shutdown
CHKONLY = .false.  ! flag for "check input only"
RSTDIA = .false.  ! flag for diagonal restart
RSTCHO = .false.  ! flag for decomposition restart
MODRST = -1       ! default restart model (if restart at all)
MAXVEC = 0        ! max. #Cholesky vectors (0=>generic)
MAXRED = 0        ! max. #reduced sets (0=>generic)
NCOL_CHK = 0        ! #columns to check (0=>all if requested)
CHO_USEABS = .false. ! do not use abs. value to set up red. sets
CHO_ADRVEC = 1       ! use WA addressing of vector files
CHO_IOVEC = 3       ! i/o model for reading Cholesky vectors
CHO_REORD = .false. ! reorder Cholesky vectors to full storage
N1_VECRD = 2       ! nominator for fraction of mem. for iovec
N2_VECRD = 3       ! denominator for fraction of mem. for iovec
N1_QUAL = 1       ! nominator for fraction of mem. for qual.
N2_QUAL = 3       ! denominator for fraction of mem. for qual.
N_SUBTR = MAXQUAL ! max. #vectors in subtraction of prev. vecs.
HALTIT = .false. ! halt execution after decomposition
IFCSEW = 2       ! get integrals directly in reduced set
FRAC_CHVBUF = 0.35d0 ! memory fraction used for vector buffer
CHO_SSCREEN = .false. ! screening in vector subtraction
SSTAU = -1.0d9    ! threshold in vector subtraction (neg=>generic)
SSNORM = 'Max' ! norm used for vector screening
SUBSCRSTAT(1) = 0.0d0 ! vector screening statistics
SUBSCRSTAT(2) = 0.0d0 ! vector screening statistics
CHO_1CENTER = .false. ! do not decompose 1-center diagonals only
CHO_NO2CENTER = .false. ! do not exclude 2-center diagonals
CHO_PRESCREEN = .true.  ! prescreen diagonal
THR_PRESCREEN = -1.0d9 ! diag. prescreen threshold (neg=>generic)
CHO_SIMRI = .false.    ! simulate RI
THR_SIMRI = -1.0d9     ! threshold for qualifying diags. in RI sim
#ifdef _DEBUGPRINT_
CHO_INTCHK = .true.  ! check integrals after decomposition
CHO_MINCHK = .true.  ! minimal integral check
CHO_TRCNEG = .true.  ! tracing of negative diagonals
CHO_DIACHK = .true.  ! check diagonals in qualified columns
TOL_DIACHK = 1.0D-10 ! tolerance for dia. checking (if requested)
CHO_TSTSCREEN = .true. ! test vector subtraction screening
#else
CHO_INTCHK = .false. ! check integrals after decomposition
CHO_MINCHK = .false. ! minimal integral check
CHO_TRCNEG = .false. ! tracing of negative diagonals
CHO_DIACHK = .false. ! check diagonals in qualified columns
TOL_DIACHK = 1.0D-14 ! tolerance for dia. checking (if requested)
CHO_TSTSCREEN = .false. ! test vector subtraction screening
#endif
BLOCKSIZE = 500  ! #vecs in each Z vector block
TRACE_IDLE = .false. ! trace idle processes

! Set default print level using molcas environment (or whatever was
! set in seward initially).
! -----------------------------------------------------------------

IPRINT = IPRINTLEVEL(-1)
if (IPRINT <= 2) IPRINT = 1
IPRINT = IPRINT-1

! Return if defaults only.
! ------------------------

FORCEPARALLEL = .false.
if (DFONLY) GO TO 2000

! Loop through keyword input.
! ---------------------------

1000 continue

! Read next keyword and get corresponding ID.
! -------------------------------------------

IDKEY = 0
call CHO_MCA_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,LUPRI)

if ((IDKEY >= 1) .and. (IDKEY <= NOPTION)) then  ! key found

  ! Branch for further processing.
  ! ------------------------------

  if (IDKEY == 1) GO TO 1
  if (IDKEY == 2) GO TO 2
  if (IDKEY == 3) GO TO 3
  if (IDKEY == 4) GO TO 4
  if (IDKEY == 5) GO TO 5
  if (IDKEY == 6) GO TO 6
  if (IDKEY == 7) GO TO 7
  if (IDKEY == 8) GO TO 8
  if (IDKEY == 9) GO TO 9
  if (IDKEY == 10) GO TO 10
  if (IDKEY == 11) GO TO 11
  if (IDKEY == 12) GO TO 12
  if (IDKEY == 13) GO TO 13
  if (IDKEY == 14) GO TO 14
  if (IDKEY == 15) GO TO 15
  if (IDKEY == 16) GO TO 16
  if (IDKEY == 17) GO TO 17
  if (IDKEY == 18) GO TO 18
  if (IDKEY == 19) GO TO 19
  if (IDKEY == 20) GO TO 20
  if (IDKEY == 21) GO TO 21
  if (IDKEY == 22) GO TO 22
  if (IDKEY == 23) GO TO 23
  if (IDKEY == 24) GO TO 24
  if (IDKEY == 25) GO TO 25
  if (IDKEY == 26) GO TO 26
  if (IDKEY == 27) GO TO 27
  if (IDKEY == 28) GO TO 28
  if (IDKEY == 29) GO TO 29
  if (IDKEY == 30) GO TO 30
  if (IDKEY == 31) GO TO 31
  if (IDKEY == 32) GO TO 32
  if (IDKEY == 33) GO TO 33
  if (IDKEY == 34) GO TO 34
  if (IDKEY == 35) GO TO 35
  if (IDKEY == 36) GO TO 36
  if (IDKEY == 37) GO TO 37
  if (IDKEY == 38) GO TO 38
  if (IDKEY == 39) GO TO 39
  if (IDKEY == 40) GO TO 40
  if (IDKEY == 41) GO TO 41
  if (IDKEY == 42) GO TO 42
  if (IDKEY == 43) GO TO 43
  if (IDKEY == 44) GO TO 44
  if (IDKEY == 45) GO TO 45
  if (IDKEY == 46) GO TO 46
  if (IDKEY == 47) GO TO 47
  if (IDKEY == 48) GO TO 48
  if (IDKEY == 49) GO TO 49
  if (IDKEY == 50) GO TO 50
  if (IDKEY == 51) GO TO 51
  if (IDKEY == 52) GO TO 52
  if (IDKEY == 53) GO TO 53
  if (IDKEY == 54) GO TO 54
  if (IDKEY == 55) GO TO 55
  if (IDKEY == 56) GO TO 56
  if (IDKEY == 57) GO TO 57
  if (IDKEY == 58) GO TO 58

  ! If this section is executed, it's most likely a bug.
  ! ----------------------------------------------------

  write(LUPRI,*) 'IDKEY = ',IDKEY,' is formally legal, but:'
  write(LUPRI,*) 'Did you forget to change the computed goto in ',SECNAM,'?'
  call CHO_QUIT('Illegal address in computed GOTO in '//SECNAM,105)

  ! Read decomposition threshold.
  ! -----------------------------

1 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,THRCOM)
  GO TO 1000

  ! Read print level.
  ! -----------------

2 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,IPRINT)
  GO TO 1000

  ! Read buffer length for diagonal calculation.
  ! --------------------------------------------

3 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,LBUF)
  GO TO 1000

  ! Read threshold for discarding initial diagonal elements.
  ! --------------------------------------------------------

4 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,THRDIAG)
  GO TO 1000

  ! Read damping factor for 1st reduced set.
  ! ----------------------------------------

5 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,DAMP(1))
  GO TO 1000

  ! Read damping factor for later reduced set.
  ! ------------------------------------------

6 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,DAMP(2))
  GO TO 1000

  ! Read span factor.
  ! -----------------

7 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,SPAN)
  SPAN = min(abs(SPAN),1.0d0)
  GO TO 1000

  ! Read minimum number of qualifieds for decomposition.
  ! ----------------------------------------------------

8 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MINQUAL)
  GO TO 1000

  ! Read maximum number of shell pair distributions that should
  ! be calculated before passing to decomposition.
  ! -----------------------------------------------------------

9 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MXSHPR)
  GO TO 1000

  ! Turn on screening.
  ! ------------------

10 continue
  SCDIAG = .true.
  GO TO 1000

  ! Turn off screening.
  ! -------------------

11 continue
  SCDIAG = .false.
  GO TO 1000

  ! Read algorithm ID for qualification procedure,
  ! 1: qualify from 1 to max.
  ! 2: qualify the largest qualifiables
  ! ----------------------------------------------

12 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,IALQUA)
  GO TO 1000

  ! Read threshold for zeroing (small) negative diagonals.
  ! Diag < THRNEG => Diag = 0.
  ! ------------------------------------------------------

13 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,THRNEG)
  GO TO 1000

  ! Read threshold for warning about negative diagonals.
  ! Diag < WARNEG => Diag = 0, but issue warning.
  ! ----------------------------------------------------

14 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,WARNEG)
  GO TO 1000

  ! Read threshold for shutdown due to negative diagonals.
  ! Diag < TOONEG => shutdown.
  ! ------------------------------------------------------

15 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,TOONEG)
  GO TO 1000

  ! Check configuration only.
  ! -------------------------

16 continue
  CHKONLY = .true.
  GO TO 1000

  ! Check all integrals after decomposition.
  ! ----------------------------------------

17 continue
  CHO_INTCHK = .true.
  CHO_MINCHK = .false.
  NCOL_CHK = 0  ! => all
  GO TO 1000

  ! Restart diagonal.
  ! -----------------

18 continue
  RSTDIA = .true.
  GO TO 1000

  ! Restart decomposition using default model.
  ! ------------------------------------------

19 continue
  if (.not. RSTCHO) then
    RSTCHO = .true.
    MODRST = -1  ! use configuration from restart file
  end if
  GO TO 1000

  ! Read restart model and set restart.
  ! -----------------------------------

20 continue
  RSTCHO = .true.
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MODRST)
  if (MODRST < 0) then
    MODRST = -1  ! use configuration from restart file
  else if (MODRST > 0) then
    MODRST = 1   ! use configuration from input
  end if
  GO TO 1000

  ! Read max. #qualifieds per symmetry.
  ! -----------------------------------

21 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MAXQUAL)
  GO TO 1000

  ! Read max. #Cholesky vectors per symmetry.
  ! -----------------------------------------

22 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MAXVEC)
  GO TO 1000

  ! Read max. #reduced sets.
  ! ------------------------

23 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,MAXRED)
  GO TO 1000

  ! Read number of integral columns to check.
  ! -----------------------------------------

24 continue
  CHO_INTCHK = .true.
  CHO_MINCHK = .false.
  LINE = GET_LN(LUNIT)
  call GET_I1(1,NCOL_CHK)
  GO TO 1000

  ! Minimal integral check.
  ! -----------------------

25 continue
  CHO_INTCHK = .true.
  CHO_MINCHK = .true.
  GO TO 1000

  ! Use abs. value to set up reduced sets.
  ! --------------------------------------

26 continue
  CHO_USEABS = .true.
  GO TO 1000

  ! Do not use abs. value to set up reduced sets.
  ! ---------------------------------------------

27 continue
  CHO_USEABS = .false.
  GO TO 1000

  ! Turn on tracing of negative diagonals.
  ! --------------------------------------

28 continue
  CHO_TRCNEG = .true.
  GO TO 1000

  ! Read I/O model used for reading Cholesky vectors.
  ! -------------------------------------------------

29 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,CHO_IOVEC)
  GO TO 1000

  ! Set vector reordering.
  ! ----------------------

30 continue
  CHO_REORD = .true.
  GO TO 1000

  ! Halt execution after decomposition.
  ! -----------------------------------

31 continue
  HALTIT = .true.
  GO TO 1000

  ! Set fraction of memory to be used as vector buffer during subtraction.
  ! ----------------------------------------------------------------------

32 continue
  LINE = GET_LN(LUNIT)
  call GET_I(1,INTSCR,2)
  N1_VECRD = INTSCR(1)
  N2_VECRD = INTSCR(2)
  GO TO 1000

  ! Set fraction of memory to be used for qualified columns.
  ! --------------------------------------------------------

33 continue
  LINE = GET_LN(LUNIT)
  call GET_I(1,INTSCR,2)
  N1_QUAL = INTSCR(1)
  N2_QUAL = INTSCR(2)
  GO TO 1000

  ! Set max. #vectors in subtraction of prev. vectors.
  ! --------------------------------------------------

34 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,N_SUBTR)
  GO TO 1000

  ! Set address mode for vector I/O: WA (=1) or DA (=2).
  ! ----------------------------------------------------

35 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,CHO_ADRVEC)
  GO TO 1000

  ! Set model used to interface to integral code.
  ! ---------------------------------------------

36 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,IFCSEW)
  if (IFCSEW < 1) then
    IFCSEW = 1
  else if (IFCSEW > 2) then
    IFCSEW = 2
  end if
  GO TO 1000

  ! One-step (conventional) decomposition.
  ! --------------------------------------

37 continue
  CHO_DECALG = 1
  DECALG_USRDEF = .true.
  GO TO 1000

  ! Two-step (generate map, then vectors) decomposition.
  ! ----------------------------------------------------

38 continue
  CHO_DECALG = 2
  DECALG_USRDEF = .true.
  GO TO 1000

  ! Naive decomposition (remove diagonals as soon as smaller than threshold).
  ! -------------------------------------------------------------------------

39 continue
  CHO_DECALG = 3
  DECALG_USRDEF = .true.
  GO TO 1000

  ! Set memory fraction used for vector buffer.
  ! -------------------------------------------

40 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,FRAC_CHVBUF)
  GO TO 1000

  ! Turn on diagonal checking of qualified columns both before
  ! and after subtraction of previous vectors. Set tolerance in checking.
  ! ---------------------------------------------------------------------

41 continue
  CHO_DIACHK = .true.
  LINE = GET_LN(LUNIT)
  call GET_F1(1,TOL_DIACHK)
  GO TO 1000

  ! Test vector subtraction screening.
  ! ----------------------------------

42 continue
  CHO_TSTSCREEN = .true.
  GO TO 1000

  ! Use screening in vector subtraction.
  ! ------------------------------------

43 continue
  CHO_SSCREEN = .true.
  GO TO 1000

  ! Do not use screening in vector subtraction.
  ! -------------------------------------------

44 continue
  CHO_SSCREEN = .false.
  GO TO 1000

  ! Threshold for screening in vector subtraction.
  ! ----------------------------------------------

45 continue
  LINE = GET_LN(LUNIT)
  call GET_F1(1,SSTAU)
  GO TO 1000

  ! Threshold for screening in vector subtraction.
  ! ----------------------------------------------

46 continue
  LINE = GET_LN(LUNIT)
  call GET_S(1,SSNORM,1)
  GO TO 1000

  ! Decompose 1-center diagonals only. 2-center diagonals are
  ! still included, but never explicitly decomposed.
  ! ---------------------------------------------------------

47 continue
  CHO_1CENTER = .true.
  GO TO 1000

  ! Decompose 1-center diagonals only. 2-center diagonals are
  ! never included.
  ! ---------------------------------------------------------

48 continue
  CHO_1CENTER = .true.
  CHO_NO2CENTER = .true.
  GO TO 1000

  ! Turn off diagonal prescreening.
  ! -------------------------------

49 continue
  CHO_PRESCREEN = .false.
  GO TO 1000

  ! Turn on diagonal prescreening.
  ! ------------------------------

50 continue
  CHO_PRESCREEN = .true.
  GO TO 1000

  ! Set diagonal prescreening threshold.
  ! ------------------------------------

51 continue
  CHO_PRESCREEN = .true.
  LINE = GET_LN(LUNIT)
  call GET_F1(1,THR_PRESCREEN)
  GO TO 1000

  ! Use "parallel" decomposition algorithm.
  ! ---------------------------------------

52 continue
  FORCEPARALLEL = .true.
  GO TO 1000

  ! Use "parallel" decomposition algorithm.
  ! ---------------------------------------

53 continue
  CHO_SIMP = .true.
  GO TO 1000

  ! Simulate RI.
  ! Activate 1-center decomposition.
  ! Read threshold for removing 1-center diagonals in map generation.
  ! -----------------------------------------------------------------

54 continue
  CHO_SIMRI = .true.
  CHO_1CENTER = .true.
  LINE = GET_LN(LUNIT)
  call GET_F1(1,THR_SIMRI)
  GO TO 1000

  ! Activate fake parallel:
  ! Run parallel decomposition, then distribute vectors if parallel.
  ! ----------------------------------------------------------------

55 continue
  CHO_FAKE_PAR = .true.
  FAKE_USRDEF = .true.
  GO TO 1000

  ! Activate true parallel decomposition (deactivate fake).
  ! -------------------------------------------------------

56 continue
  CHO_FAKE_PAR = .false.
  FAKE_USRDEF = .true.
  GO TO 1000

  ! Read BlockSize (Z vectors).
  ! ---------------------------

57 continue
  LINE = GET_LN(LUNIT)
  call GET_I1(1,BLOCKSIZE)
  GO TO 1000

  ! Turn on tracing of idle nodes.
  ! ------------------------------

58 continue
  TRACE_IDLE = .true.
  GO TO 1000

else if (IDKEY == NOPTION+1) then  ! end of input
  GO TO 2000
else if (IDKEY == -5) then ! some internal error
  GO TO 1300
else if (IDKEY == -1) then ! keyword error
  GO TO 1200
else  ! catch all out-of-bound errors
  GO TO 1100
end if

! Post processing.
! ----------------

2000 continue
if (FAKE_USRDEF) then
  call CHO_PARCONF(CHO_FAKE_PAR) ! reset parallel configuration
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

call PUT_THR_CHO(THRCOM)

! Normal exit.
! ------------

return

! Error exits.
! ------------

1100 continue ! IDKEY out of bounds
write(LUPRI,*) SECNAM,': IDKEY out of bounds: ',IDKEY
IRC = 104
GO TO 9999

1200 continue ! keyword error
write(LUPRI,*) SECNAM,': keyword error detected, IDKEY = ',IDKEY
IRC = 105
GO TO 9999

1300 continue ! some internal error
write(LUPRI,*) SECNAM,': internal error detected, IDKEY = ',IDKEY
IRC = 103
GO TO 9999

9999 call CHO_QUIT('Error in '//SECNAM,IRC)

return

#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_character(Line)
#endif

end subroutine CHO_INP
