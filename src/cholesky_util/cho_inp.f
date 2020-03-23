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
      SUBROUTINE CHO_INP(DFONLY,LUNIT,LUOUT)
C
C     Purpose: If DFONLY, set defaults only.
C              Else, read and process input for Cholesky decomposition
C              from unit LUNIT. LUOUT is the unit of the output file
C              which is stored internally in the Cholesky program as
C              LUPRI (in cholesky.fh).
C
#include "implicit.fh"
      LOGICAL DFONLY
#include "cholesky.fh"
#include "chosubscr.fh"
#include "choprint.fh"
#include "chosimri.fh"

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CHO_INP')

      LOGICAL DECALG_USRDEF, MXSHPR_USRDEF, FAKE_USRDEF
      LOGICAL FORCEPARALLEL

      CHARACTER*180 LINE
      CHARACTER*180 GET_LN
      EXTERNAL GET_LN

      INTEGER  IPRINTLEVEL
      EXTERNAL IPRINTLEVEL

      PARAMETER (LINTSCR = 2)
      INTEGER INTSCR(LINTSCR)

      PARAMETER (NOPTION = 58, LOPTION = 50)  ! #options
      CHARACTER*(LOPTION) OPTION(NOPTION)

      DATA OPTION /'set decomposition threshold                       ',
     &             'set print level                                   ',
     &             'set buffer length for diagonal calculation        ',
     &             'set threshold for discarding initial diag. elem.  ',
     &             'set damping factor for first reduced set          ',
     &             'set damping factor for later reduced sets         ',
     &             'set span                                          ',
     &             'set minimum number of qualifieds for decomposition',
     &             'set maximum number of shell pair distributions    ',
     &             'turn on diag. screening (using damping factors)   ',
     &             'turn off diag. screening                          ',
     &             'set qualification algorithm                       ',
     &             'set threshold for zeroing negative diagonals      ',
     &             'set threshold for warning about negative diagonals',
     &             'set threshold for shutdown due to negative diag.  ',
     &             'check configuration only                          ',
     &             'check all integrals (debug)                       ',
     &             'diagonal restart                                  ',
     &             'decomposition restart with default restart model  ',
     &             'set restart model                                 ',
     &             'set maximum #qualifieds per symmetry              ',
     &             'set maximum #Cholesky vectors per symmetry        ',
     &             'set max. #reduced sets (i.e., integral passes)    ',
     &             'check specified #columns (debug)                  ',
     &             'minimal integral check (debug)                    ',
     &             'use abs. value to set up reduced sets             ',
     &             'do not use abs. value to set up red. sets         ',
     &             'turn on tracing of negative diagonals (debug)     ',
     &             'set algorithm for reading Cholesky vectors        ',
     &             'set vector reordering (for int. re-generation)    ',
     &             'halt execution after decomposition                ',
     &             'set memory fraction used as buffer for vec. read  ',
     &             'set fraction of memory to be used for qual. col.  ',
     &             'set maximum number of vectors in subtraction      ',
     &             'set file address mode (WA or DA for vector files) ',
     &             'set model used to interface to integral code      ',
     &             'one-step (conventional) decomposition             ',
     &             'two-step (generate map, then vectors) decomp.     ',
     &             'naive decomposition                               ',
     &             'set memory fraction used for global vector buffer ',
     &             'set diagonal checking and set tolerance (debug)   ',
     &             'test vector subtraction screening (statistics)    ',
     &             'use screening in vector subtraction               ',
     &             'do not use screening in vector subtraction        ',
     &             'threshold for screening in vector subtraction     ',
     &             'norm to use in vector subtraction screening       ',
     &             'decompose 1-center diagonals only                 ',
     &             'decompose 1-center diagonals only; excl. 2-center ',
     &             'turn off diagonal prescreening                    ',
     &             'turn on diagonal prescreening                     ',
     &             'set diagonal prescreening threshold               ',
     &             'use parallel decomposition algorithm              ',
     &             'simulate parallel algorithm (i.e. modified span)  ',
     &             'simulate RI (input 1-center threshold)            ',
     &             'activate fake parallel decomposition              ',
     &             'activate true parallel decomposition (deact. fake)',
     &             'set block size for Z vectors (parallel two-step)  ',
     &             'activate tracing of idle nodes (parallel run)     '/

C     Define all entries in common blocks.
C     Set output unit.
C     Set run mode to "internal".
C     ------------------------------------

      IRC = 0
      CALL CHO_X_SETINC(IRC)
      LUPRI = LUOUT
      RUN_MODE = RUN_INTERNAL
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_X_SETINC returned error code ',IRC
         WRITE(LUPRI,*) '(most likely due to a programming error...)'
         CALL CHO_QUIT('Include file initialization error in '//SECNAM,
     &                 102)
      END IF
      CALL CHO_SETPTR2()
      CALL CHO_SETGLOB()

C     Set default parallel configuration.
C     CHO_FAKE_PAR = .TRUE. : fake parallel.
C     CHO_FAKE_PAR = .FALSE.: true parallel.
C     --------------------------------------

      CHO_FAKE_PAR = .FALSE.
      CALL CHO_PARCONF(CHO_FAKE_PAR)

C     Set default decomposition algorithm (depends on parallel or serial
C     installation).
C     ------------------------------------------------------------------

      CALL CHO_SETDECALG_DEF()

C     Set default MXSHPR (depends on parallel or serial).
C     ---------------------------------------------------

      CALL CHO_SETMXSHPR_DEF(MXSHPR_DEF)

C     Set defaults.
C     -------------

      DECALG_USRDEF = .FALSE.
      MXSHPR_USRDEF = .FALSE.
      FAKE_USRDEF = .FALSE.

      CHO_DECALG = CHO_DECALG_DEF ! default decomposition algorithm
      CHO_SIMP = .FALSE. ! do not simulate parallel algorithm
      THRCOM  = THRDEF   ! default decomposition threshold
      LBUF    = 1000000  ! buffer length for initial diagonal
      THRDIAG = 0.0D0    ! threshold for initial diag. screening
      XLDIAG  = 0.0D0    ! just an initialization
      DAMP(1) = -1.0D9   ! damping for initial screening (neg.=>generic)
      DAMP(2) = -1.0D9   ! damping for later screenings (neg.=>generic)
      SPAN    = 1.0D-2   ! span factor
      MINQUAL = 50       ! min. #qual. needed for proceeding
      MAXQUAL = 100      ! max. #qual. per symmetry
      MXSHPR  = MXSHPR_DEF ! max. #sh. pairs for proceeding (0=>generic)
      SCDIAG  = .TRUE.   ! screen diagonal during decom.
      IALQUA  = 2        ! qualification algorithm
      THRNEG  = -1.0D-40 ! diag<THRNEG => diag=0
      WARNEG  = -1.0D-10 ! diag<WARNEG => diag=0, issue warning
      TOONEG  = -1.0D-8  ! diag<TOONEG => shutdown
      CHKONLY = .FALSE.  ! flag for "check input only"
      RSTDIA  = .FALSE.  ! flag for diagonal restart
      RSTCHO  = .FALSE.  ! flag for decomposition restart
      MODRST  = -1       ! default restart model (if restart at all)
      MAXVEC  = 0        ! max. #Cholesky vectors (0=>generic)
      MAXRED  = 0        ! max. #reduced sets (0=>generic)
      NCOL_CHK= 0        ! #columns to check (0=>all if requested)
      CHO_USEABS = .FALSE. ! do not use abs. value to set up red. sets
      CHO_ADRVEC = 1       ! use WA addressing of vector files
      CHO_IOVEC  = 3       ! i/o model for reading Cholesky vectors
      CHO_REORD  = .FALSE. ! reorder Cholesky vectors to full storage
      N1_VECRD   = 2       ! nominator for fraction of mem. for iovec
      N2_VECRD   = 3       ! denominator for fraction of mem. for iovec
      N1_QUAL    = 1       ! nominator for fraction of mem. for qual.
      N2_QUAL    = 3       ! denominator for fraction of mem. for qual.
      N_SUBTR    = MAXQUAL ! max. #vectors in subtraction of prev. vecs.
      HALTIT     = .FALSE. ! halt execution after decomposition
      IFCSEW     = 2       ! get integrals directly in reduced set
      FRAC_CHVBUF = 0.35D0 ! memory fraction used for vector buffer
      CHO_SSCREEN = .FALSE. ! screening in vector subtraction
      SSTAU = -1.0d9    ! threshold in vector subtraction (neg=>generic)
      SSNORM = 'Max' ! norm used for vector screening
      SUBSCRSTAT(1) = 0.0D0 ! vector screening statistics
      SUBSCRSTAT(2) = 0.0D0 ! vector screening statistics
      CHO_1CENTER = .FALSE. ! do not decompose 1-center diagonals only
      CHO_NO2CENTER = .FALSE. ! do not exclude 2-center diagonals
      CHO_PRESCREEN = .TRUE.  ! prescreen diagonal
      THR_PRESCREEN = -1.0d9 ! diag. prescreen threshold (neg=>generic)
      CHO_SIMRI = .FALSE.    ! simulate RI
      THR_SIMRI = -1.0D9     ! threshold for qualifying diags. in RI sim
#if defined (_DEBUG_)
      CHO_INTCHK = .TRUE.  ! check integrals after decomposition
      CHO_MINCHK = .TRUE.  ! minimal integral check
      CHO_TRCNEG = .TRUE.  ! tracing of negative diagonals
      CHO_DIACHK = .TRUE.  ! check diagonals in qualified columns
      TOL_DIACHK = 1.0D-10 ! tolerance for dia. checking (if requested)
      CHO_TSTSCREEN = .TRUE. ! test vector subtraction screening
#else
      CHO_INTCHK = .FALSE. ! check integrals after decomposition
      CHO_MINCHK = .FALSE. ! minimal integral check
      CHO_TRCNEG = .FALSE. ! tracing of negative diagonals
      CHO_DIACHK = .FALSE. ! check diagonals in qualified columns
      TOL_DIACHK = 1.0D-14 ! tolerance for dia. checking (if requested)
      CHO_TSTSCREEN = .FALSE. ! test vector subtraction screening
#endif
      BLOCKSIZE = 500  ! #vecs in each Z vector block
      TRACE_IDLE = .FALSE. ! trace idle processes

C     Set default print level using molcas environment (or whatever was
C     set in seward initially).
C     -----------------------------------------------------------------

      IPRINT = IPRINTLEVEL(-1)
      IF (IPRINT .LE. 2) IPRINT = 1
      IPRINT = IPRINT - 1

C     Return if defaults only.
C     ------------------------

      FORCEPARALLEL = .FALSE.
      IF (DFONLY) GO TO 2000

C     Loop through keyword input.
C     ---------------------------

 1000 CONTINUE

C        Read next keyword and get corresponding ID.
C        -------------------------------------------

         IDKEY = 0
         CALL CHO_MCA_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,LUPRI)

         IF (IDKEY.GE.1 .AND. IDKEY.LE.NOPTION) THEN  ! key found

C           Branch for further processing.
C           ------------------------------

            IF (IDKEY .EQ.  1) GO TO  1
            IF (IDKEY .EQ.  2) GO TO  2
            IF (IDKEY .EQ.  3) GO TO  3
            IF (IDKEY .EQ.  4) GO TO  4
            IF (IDKEY .EQ.  5) GO TO  5
            IF (IDKEY .EQ.  6) GO TO  6
            IF (IDKEY .EQ.  7) GO TO  7
            IF (IDKEY .EQ.  8) GO TO  8
            IF (IDKEY .EQ.  9) GO TO  9
            IF (IDKEY .EQ. 10) GO TO 10
            IF (IDKEY .EQ. 11) GO TO 11
            IF (IDKEY .EQ. 12) GO TO 12
            IF (IDKEY .EQ. 13) GO TO 13
            IF (IDKEY .EQ. 14) GO TO 14
            IF (IDKEY .EQ. 15) GO TO 15
            IF (IDKEY .EQ. 16) GO TO 16
            IF (IDKEY .EQ. 17) GO TO 17
            IF (IDKEY .EQ. 18) GO TO 18
            IF (IDKEY .EQ. 19) GO TO 19
            IF (IDKEY .EQ. 20) GO TO 20
            IF (IDKEY .EQ. 21) GO TO 21
            IF (IDKEY .EQ. 22) GO TO 22
            IF (IDKEY .EQ. 23) GO TO 23
            IF (IDKEY .EQ. 24) GO TO 24
            IF (IDKEY .EQ. 25) GO TO 25
            IF (IDKEY .EQ. 26) GO TO 26
            IF (IDKEY .EQ. 27) GO TO 27
            IF (IDKEY .EQ. 28) GO TO 28
            IF (IDKEY .EQ. 29) GO TO 29
            IF (IDKEY .EQ. 30) GO TO 30
            IF (IDKEY .EQ. 31) GO TO 31
            IF (IDKEY .EQ. 32) GO TO 32
            IF (IDKEY .EQ. 33) GO TO 33
            IF (IDKEY .EQ. 34) GO TO 34
            IF (IDKEY .EQ. 35) GO TO 35
            IF (IDKEY .EQ. 36) GO TO 36
            IF (IDKEY .EQ. 37) GO TO 37
            IF (IDKEY .EQ. 38) GO TO 38
            IF (IDKEY .EQ. 39) GO TO 39
            IF (IDKEY .EQ. 40) GO TO 40
            IF (IDKEY .EQ. 41) GO TO 41
            IF (IDKEY .EQ. 42) GO TO 42
            IF (IDKEY .EQ. 43) GO TO 43
            IF (IDKEY .EQ. 44) GO TO 44
            IF (IDKEY .EQ. 45) GO TO 45
            IF (IDKEY .EQ. 46) GO TO 46
            IF (IDKEY .EQ. 47) GO TO 47
            IF (IDKEY .EQ. 48) GO TO 48
            IF (IDKEY .EQ. 49) GO TO 49
            IF (IDKEY .EQ. 50) GO TO 50
            IF (IDKEY .EQ. 51) GO TO 51
            IF (IDKEY .EQ. 52) GO TO 52
            IF (IDKEY .EQ. 53) GO TO 53
            IF (IDKEY .EQ. 54) GO TO 54
            IF (IDKEY .EQ. 55) GO TO 55
            IF (IDKEY .EQ. 56) GO TO 56
            IF (IDKEY .EQ. 57) GO TO 57
            IF (IDKEY .EQ. 58) GO TO 58

C           If this section is executed, it's most likely a bug.
C           ----------------------------------------------------

            WRITE(LUPRI,*) 'IDKEY = ',IDKEY,' is formally legal, but:'
            WRITE(LUPRI,*) 'Did you forget to change the ',
     &                     'computed goto in ',SECNAM,'?'
            CALL CHO_QUIT('Illegal address in computed GOTO in '
     &                    //SECNAM,105)

C           Read decomposition threshold.
C           -----------------------------

    1       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,THRCOM)
            GO TO 1000

C           Read print level.
C           -----------------

    2       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,IPRINT)
            GO TO 1000

C           Read buffer length for diagonal calculation.
C           --------------------------------------------

    3       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,LBUF)
            GO TO 1000

C           Read threshold for discarding initial diagonal elements.
C           --------------------------------------------------------

    4       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,THRDIAG)
            GO TO 1000

C           Read damping factor for 1st reduced set.
C           ----------------------------------------

    5       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,DAMP(1))
            GO TO 1000

C           Read damping factor for later reduced set.
C           ------------------------------------------

    6       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,DAMP(2))
            GO TO 1000

C           Read span factor.
C           -----------------

    7       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,SPAN)
               SPAN=MIN(ABS(SPAN),1.0d0)
            GO TO 1000

C           Read minimum number of qualifieds for decomposition.
C           ----------------------------------------------------

    8       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MINQUAL)
            GO TO 1000

C           Read maximum number of shell pair distributions that should
C           be calculated before passing to decomposition.
C           -----------------------------------------------------------

    9       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MXSHPR)
            GO TO 1000

C           Turn on screening.
C           ------------------

   10       CONTINUE
               SCDIAG = .TRUE.
            GO TO 1000

C           Turn off screening.
C           -------------------

   11       CONTINUE
               SCDIAG = .FALSE.
            GO TO 1000

C           Read algorithm ID for qualification procedure,
C           1: qualify from 1 to max.
C           2: qualify the largest qualifiables
C           ----------------------------------------------

   12       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,IALQUA)
            GO TO 1000

C           Read threshold for zeroing (small) negative diagonals.
C           Diag < THRNEG => Diag = 0.
C           ------------------------------------------------------

   13       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,THRNEG)
            GO TO 1000

C           Read threshold for warning about negative diagonals.
C           Diag < WARNEG => Diag = 0, but issue warning.
C           ----------------------------------------------------

   14       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,WARNEG)
            GO TO 1000

C           Read threshold for shutdown due to negative diagonals.
C           Diag < TOONEG => shutdown.
C           ------------------------------------------------------

   15       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,TOONEG)
            GO TO 1000

C           Check configuration only.
C           -------------------------

   16       CONTINUE
               CHKONLY = .TRUE.
            GO TO 1000

C           Check all integrals after decomposition.
C           ----------------------------------------

   17       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .FALSE.
               NCOL_CHK   = 0  ! => all
            GO TO 1000

C           Restart diagonal.
C           -----------------

   18       CONTINUE
               RSTDIA = .TRUE.
            GO TO 1000

C           Restart decomposition using default model.
C           ------------------------------------------

   19       CONTINUE
               IF (.NOT. RSTCHO) THEN
                  RSTCHO = .TRUE.
                  MODRST = -1  ! use configuration from restart file
               END IF
            GO TO 1000

C           Read restart model and set restart.
C           -----------------------------------

   20       CONTINUE
               RSTCHO = .TRUE.
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MODRST)
               IF (MODRST .LT. 0) THEN
                  MODRST = -1  ! use configuration from restart file
               ELSE IF (MODRST .GT. 0) THEN
                  MODRST = 1   ! use configuration from input
               END IF
            GO TO 1000

C           Read max. #qualifieds per symmetry.
C           -----------------------------------

   21       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MAXQUAL)
            GO TO 1000

C           Read max. #Cholesky vectors per symmetry.
C           -----------------------------------------

   22       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MAXVEC)
            GO TO 1000

C           Read max. #reduced sets.
C           ------------------------

   23       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,MAXRED)
            GO TO 1000

C           Read number of integral columns to check.
C           -----------------------------------------

   24       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .FALSE.
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,NCOL_CHK)
            GO TO 1000

C           Minimal integral check.
C           -----------------------

   25       CONTINUE
               CHO_INTCHK = .TRUE.
               CHO_MINCHK = .TRUE.
            GO TO 1000

C           Use abs. value to set up reduced sets.
C           --------------------------------------

   26       CONTINUE
               CHO_USEABS = .TRUE.
            GO TO 1000

C           Do not use abs. value to set up reduced sets.
C           ---------------------------------------------

   27       CONTINUE
               CHO_USEABS = .FALSE.
            GO TO 1000

C           Turn on tracing of negative diagonals.
C           --------------------------------------

   28       CONTINUE
               CHO_TRCNEG = .TRUE.
            GO TO 1000

C           Read I/O model used for reading Cholesky vectors.
C           -------------------------------------------------

   29       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,CHO_IOVEC)
            GO TO 1000

C           Set vector reordering.
C           ----------------------

   30       CONTINUE
               CHO_REORD = .TRUE.
            GO TO 1000

C           Halt execution after decomposition.
C           -----------------------------------

   31       CONTINUE
               HALTIT = .TRUE.
            GO TO 1000

C           Set fraction of memory to be used as vector buffer during
C           subtraction.
C           ---------------------------------------------------------

   32       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I(1,INTSCR,2)
               N1_VECRD = INTSCR(1)
               N2_VECRD = INTSCR(2)
            GO TO 1000

C           Set fraction of memory to be used for qualified columns.
C           --------------------------------------------------------

   33       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I(1,INTSCR,2)
               N1_QUAL = INTSCR(1)
               N2_QUAL = INTSCR(2)
            GO TO 1000

C           Set max. #vectors in subtraction of prev. vectors.
C           --------------------------------------------------

   34       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,N_SUBTR)
            GO TO 1000

C           Set address mode for vector I/O: WA (=1) or DA (=2).
C           ----------------------------------------------------

   35       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,CHO_ADRVEC)
            GO TO 1000

C           Set model used to interface to integral code.
C           ---------------------------------------------

   36       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_I1(1,IFCSEW)
               IF (IFCSEW .LT. 1) THEN
                  IFCSEW = 1
               ELSE IF (IFCSEW .GT. 2) THEN
                  IFCSEW = 2
               END IF
            GO TO 1000

C           One-step (conventional) decomposition.
C           --------------------------------------

   37       CONTINUE
               CHO_DECALG = 1
               DECALG_USRDEF = .TRUE.
            GO TO 1000

C           Two-step (generate map, then vectors) decomposition.
C           ----------------------------------------------------

   38       CONTINUE
               CHO_DECALG = 2
               DECALG_USRDEF = .TRUE.
            GO TO 1000

C           Naive decomposition (remove diagonals as soon as smaller
C           than threshold).
C           --------------------------------------------------------

   39       CONTINUE
               CHO_DECALG = 3
               DECALG_USRDEF = .TRUE.
            GO TO 1000

C           Set memory fraction used for vector buffer.
C           -------------------------------------------

   40       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,FRAC_CHVBUF)
            GO TO 1000

C           Turn on diagonal checking of qualified columns both before
C           and after subtraction of previous vectors. Set tolerance in
C           checking.
C           -----------------------------------------------------------

   41       CONTINUE
               CHO_DIACHK = .TRUE.
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,TOL_DIACHK)
            GO TO 1000

C           Test vector subtraction screening.
C           ----------------------------------

   42       CONTINUE
               CHO_TSTSCREEN = .TRUE.
            GO TO 1000

C           Use screening in vector subtraction.
C           ------------------------------------

   43       CONTINUE
               CHO_SSCREEN = .TRUE.
            GO TO 1000

C           Do not use screening in vector subtraction.
C           -------------------------------------------

   44       CONTINUE
               CHO_SSCREEN = .FALSE.
            GO TO 1000

C           Threshold for screening in vector subtraction.
C           ----------------------------------------------

   45       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,SSTAU)
            GO TO 1000

C           Threshold for screening in vector subtraction.
C           ----------------------------------------------

   46       CONTINUE
               LINE = GET_LN(LUNIT)
               CALL GET_S(1,SSNORM,1)
            GO TO 1000

C           Decompose 1-center diagonals only. 2-center diagonals are
C           still included, but never explicitly decomposed.
C           ---------------------------------------------------------

   47       CONTINUE
               CHO_1CENTER = .TRUE.
            GO TO 1000

C           Decompose 1-center diagonals only. 2-center diagonals are
C           never included.
C           ---------------------------------------------------------

   48       CONTINUE
               CHO_1CENTER = .TRUE.
               CHO_NO2CENTER = .TRUE.
            GO TO 1000

C           Turn off diagonal prescreening.
C           -------------------------------

   49       CONTINUE
               CHO_PRESCREEN = .FALSE.
            GO TO 1000

C           Turn on diagonal prescreening.
C           ------------------------------

   50       CONTINUE
               CHO_PRESCREEN = .TRUE.
            GO TO 1000

C           Set diagonal prescreening threshold.
C           ------------------------------------

   51       CONTINUE
               CHO_PRESCREEN = .TRUE.
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,THR_PRESCREEN)
            GO TO 1000

C           Use "parallel" decomposition algorithm.
C           ---------------------------------------

   52       CONTINUE
               FORCEPARALLEL = .TRUE.
            GO TO 1000

C           Use "parallel" decomposition algorithm.
C           ---------------------------------------

   53       CONTINUE
               CHO_SIMP = .TRUE.
            GO TO 1000

C           Simulate RI.
C           Activate 1-center decomposition.
C           Read threshold for removing 1-center diagonals in map
C           generation.
C           -----------------------------------------------------

   54       CONTINUE
               CHO_SIMRI = .TRUE.
               CHO_1CENTER = .TRUE.
               LINE = GET_LN(LUNIT)
               CALL GET_F1(1,THR_SIMRI)
            GO TO 1000

C           Activate fake parallel:
C           Run parallel decomposition, then distribute vectors if
C           parallel.
C           ------------------------------------------------------

   55       CONTINUE
               CHO_FAKE_PAR = .TRUE.
               FAKE_USRDEF = .TRUE.
            GO TO 1000

C           Activate true parallel decomposition (deactivate fake).
C           -------------------------------------------------------

   56       CONTINUE
               CHO_FAKE_PAR = .FALSE.
               FAKE_USRDEF = .TRUE.
            GO TO 1000

C           Read BlockSize (Z vectors).
C           ---------------------------

   57       CONTINUE
                LINE = GET_LN(LUNIT)
                CALL GET_I1(1,BLOCKSIZE)
            GO TO 1000

C           Turn on tracing of idle nodes.
C           ------------------------------

   58       CONTINUE
               TRACE_IDLE = .TRUE.
            GO TO 1000

         ELSE IF (IDKEY .EQ. NOPTION+1) THEN  ! end of input
            GO TO 2000
         ELSE IF (IDKEY .EQ. -5) THEN ! some internal error
            GO TO 1300
         ELSE IF (IDKEY .EQ. -1) THEN ! keyword error
            GO TO 1200
         ELSE  ! catch all out-of-bound errors
            GO TO 1100
         END IF

C     Post processing.
C     ----------------

 2000 CONTINUE
      IF (FAKE_USRDEF) THEN
         CALL CHO_PARCONF(CHO_FAKE_PAR) ! reset parallel configuration
         IF (.NOT.DECALG_USRDEF) THEN ! reset default decomposition alg
            CALL CHO_SETDECALG_DEF()
            CHO_DECALG = CHO_DECALG_DEF
         END IF
         IF (.NOT.MXSHPR_USRDEF) THEN ! reset default MXSHPR
            CALL CHO_SETMXSHPR_DEF(MXSHPR_DEF)
            MXSHPR = MXSHPR_DEF
         END IF
      END IF
      CALL CHO_INP_SETDECALG(FORCEPARALLEL)
      IF (.NOT.MXSHPR_USRDEF) THEN
         IF (CHO_DECALG.EQ.4 .OR. CHO_DECALG.EQ.5 .OR.
     &       CHO_DECALG.EQ.6) THEN
            MXSHPR=1
         ELSE
            MXSHPR=0
         END IF
      END IF

C     Put decomposition threshold in info.fh.
C     ---------------------------------------

      CALL PUT_THR_CHO(THRCOM)

C     Normal exit.
C     ------------

      RETURN

C     Error exits.
C     ------------

 1100 CONTINUE ! IDKEY out of bounds
      WRITE(LUPRI,*) SECNAM,': IDKEY out of bounds: ',IDKEY
      IRC = 104
      GO TO 9999

 1200 CONTINUE ! keyword error
      WRITE(LUPRI,*) SECNAM,': keyword error detected, IDKEY = ',IDKEY
      IRC = 105
      GO TO 9999

 1300 CONTINUE ! some internal error
      WRITE(LUPRI,*) SECNAM,': internal error detected, IDKEY = ',IDKEY
      IRC = 103
      GO TO 9999

 9999 CALL CHO_QUIT('Error in '//SECNAM,IRC)
      END
************************************************************************
*                                                                      *
************************************************************************
      Subroutine Put_thr_Cho(ThrCom)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
*
      If (.not. Do_RI) Thrshld_CD = ThrCom
*
      Return
      End
************************************************************************
*                                                                      *
************************************************************************
      SubRoutine Cho_Inp_SetDecAlg(ForceParallel)
      Implicit None
      Logical ForceParallel
#include "cholesky.fh"
#include "cho_para_info.fh"

      If (Cho_Real_Par .or. ForceParallel) Then
         If (Cho_DecAlg .eq. 1) Then
            Cho_DecAlg = 4  ! parallel one-step
         Else If (Cho_DecAlg .eq. 2) Then
            Cho_DecAlg = 5 ! parallel two-step
         Else If (Cho_DecAlg .eq. 3) Then
            Cho_DecAlg = 6  ! parallel naive
         End If
      End If

      End
