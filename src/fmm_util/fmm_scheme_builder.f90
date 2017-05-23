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
MODULE fmm_scheme_builder

   USE fmm_global_paras
   USE fmm_stats

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_scheme,  &
             fmm_get_scheme

   ! "scheme" contains all the input (and default) parameter data to
   ! completely define a unique MM run-type
   TYPE(scheme_paras), TARGET, SAVE :: scheme

   ! Safety check for initialisation
   LOGICAL, SAVE :: scheme_initialised = .FALSE.

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_scheme(scheme_ptr)

      IMPLICIT NONE
      TYPE(scheme_paras), POINTER :: scheme_ptr
      INTEGER(INTK) :: iteration = 0

      IF (.NOT. scheme_initialised) CALL fmm_quit('fmm scheme uninitialised!')

      NULLIFY(scheme_ptr)
      scheme_ptr => scheme

      iteration = iteration + 1
      stat_iteration = iteration

   END SUBROUTINE fmm_get_scheme

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_scheme(job_type)

!      USE fmm_md4_globals, ONLY: fmm_grain_inv, fmm_extent_min

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: job_type

      INTEGER(INTK) :: flag
      INTEGER(INTK) :: LMAX, TLMAX, ALGORITHM, FEdim, LIPN
      REAL(REALK)   :: GRAIN, DENS_SCREEN, EXTENT_MIN

      INTEGER, PARAMETER :: IO = 5
      NAMELIST /FMM/ LMAX, TLMAX, ALGORITHM, GRAIN, DENS_SCREEN,   &
                     EXTENT_MIN, FEdim, LIPN

      !
      ! Set default scheme
      !
      scheme%job_type = job_type
      scheme%branch_free = BRFREE_DF
      scheme%pack_LHS = PACK_LHS_DF
      scheme%pack_RHS = PACK_RHS_DF
      ! We assume a symmetric T-matrix in all cases so
      ! use auxiliary (prefactor-adapted) moment types for RHS.
      scheme%T_con%LHS_mm_type = USE_RAW_QLM
      scheme%T_con%RHS_mm_type = USE_T_SYM_QLM

      !
      ! Change defaults with user input
      !
!fixme Put defaults in parameter module
      ALGORITHM = DO_FMM
      LMAX  = 4
      TLMAX = 12
      GRAIN = 1.0_REALK
      DENS_SCREEN = 1e-15_REALK
      EXTENT_MIN = EXTENT_MIN_DF
      FEdim = 10
      LIPN = 2
      REWIND(IO)
      READ(IO, FMM, END=9010, ERR=9000)
      GO TO 9010
 9000 WRITE(6,*) 'o Check NAMELIST FMM'
      CALL Abend()
 9010 CONTINUE
      scheme%algorithm = ALGORITHM
      scheme%raw_LMAX = LMAX
      scheme%trans_LMAX = TLMAX
      scheme%grain = GRAIN
!      fmm_grain_inv = one/scheme%grain
      scheme%dens_screen_thr = DENS_SCREEN
      scheme%extent_min = EXTENT_MIN
      scheme%FEdim = FEdim
      scheme%lipn = LIPN

      !
      ! Set up (W) translator and (T) contractor options
      !
      SELECT CASE( scheme%job_type )
         CASE( GFC_FMM )
            scheme%include_near_field = .TRUE.
            scheme%W_con%W_buffer = TREE_W_BUFFER
            scheme%W_con%ID = W_CONTRACTOR_FAST
            scheme%W_con%BR_ID = W_CONTRACTOR_BOUNDARY
            scheme%W_con%sort_para = SORT_BY_SCALE
            scheme%T_con%NF_ID = T_CONTRACTOR_BOUNDARY
            scheme%T_con%NF_T_buffer = NULL_T_BUFFER
            IF (scheme%algorithm == DO_FQ) THEN
               scheme%T_con%FF_ID = T_CONTRACTOR_BOUNDARY
               scheme%T_con%FF_T_buffer = NULL_T_BUFFER
            ELSE
               scheme%T_con%FF_T_buffer = SCALE_T_BUFFER
               scheme%T_con%FF_ID = T_CONTRACTOR_SCALE
            END IF
         CASE( MD4_FMM, FE_FMM )
            scheme%include_near_field = .FALSE.
            scheme%W_con%W_buffer = TREE_W_BUFFER
            scheme%W_con%ID = W_CONTRACTOR_FAST
            scheme%W_con%BR_ID = W_CONTRACTOR_FAST
            scheme%W_con%sort_para = SORT_BY_SCALE
            scheme%T_con%NF_T_buffer = NULL_T_BUFFER
            scheme%T_con%NF_ID = T_CONTRACTOR_FULL
            flag = 2
            IF (scheme%algorithm == DO_FQ) flag = 1
            SELECT CASE( flag )
               CASE( 0 )
                  ! we do no contractions (for diagnostics)
                  scheme%T_con%FF_T_buffer = SKIP_T_BUFFER
                  scheme%T_con%FF_ID = T_CONTRACTOR_DIRECT
               CASE( 1 )
                  scheme%T_con%FF_T_buffer = NULL_T_BUFFER
                  scheme%T_con%FF_ID = T_CONTRACTOR_FULL
               CASE( 2 )
                  scheme%T_con%FF_T_buffer = SCALE_T_BUFFER
                  scheme%T_con%FF_ID = T_CONTRACTOR_SCALE
      !            scheme%T_con%sort_para = SORT_BY_SCALE
      !         CASE( 3 )
      !            scheme%T_con%T_buffer = TREE_T_BUFFER
      !            scheme%T_con%sort_para = SORT_BY_SCALE
      !            scheme%T_con%ID = T_CONTRACTOR_TREE
      !         CASE( 4 )
      !            scheme%T_con%T_buffer = TREE_T_BUFFER
      !            scheme%T_con%ID = T_CONTRACTOR_SCALE_TREE
               CASE DEFAULT
                  CALL fmm_quit('invalid T-contractor specified!')
            END SELECT
         CASE DEFAULT
            CALL fmm_quit('invalid FMM run-type requested!')
      END SELECT

      CALL fmm_verify_scheme
      CALL fmm_print_scheme
      scheme_initialised = .TRUE.

   END SUBROUTINE fmm_init_scheme

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_verify_scheme

      IMPLICIT NONE

      SELECT CASE( scheme%algorithm )
      CASE (MD4_FMM)
         IF (WS_MIN < 2 * CEILING(scheme%extent_min/scheme%grain*half)) THEN
            WRITE(LUPRI,*) 'WS_MIN = ', WS_MIN
            WRITE(LUPRI,*) 'Extent_min = ', scheme%extent_min
            WRITE(LUPRI,*) 'Grain  = ', scheme%grain
            CALL fmm_quit('RPQ cut off too large or boxes too small!')
         END IF
      CASE DEFAULT
         CONTINUE
      END SELECT

      IF (scheme%trans_LMAX < scheme%raw_LMAX) CALL fmm_quit('increase TLMAX!')

   END SUBROUTINE fmm_verify_scheme

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_print_scheme

!      USE fmm_md4_globals

      IMPLICIT NONE

      WRITE(LUPRI,'(/,A)') " -----------------------------------------"
      WRITE(LUPRI,'(A)')   " |  Multipole module runtime parameters  |"
      WRITE(LUPRI,'(A,/)') " -----------------------------------------"

      SELECT CASE( scheme%job_type )
         CASE( GFC_FMM )
            WRITE(LUPRI,*) "Computing classical boundary potential."
         CASE( MD4_FMM )
            WRITE(LUPRI,*) "Computing multipole contribution to J-matrix."
         CASE( FE_FMM )
            WRITE(LUPRI,*) "Computing full J-matrix via FE-FMM."
         CASE DEFAULT
            CALL fmm_quit ('MM job type not recognised!')
      END SELECT

      SELECT CASE( scheme%algorithm )
         CASE( DO_NULL )
            WRITE(LUPRI,*) "Skipping all FF interactions."
         CASE( DO_FQ )
            WRITE(LUPRI,*) "Running simple O(N^2) algorithm."
         CASE( DO_BQ )
            WRITE(LUPRI,*) "Running fast O(N^2) algorithm with boxes."
         CASE( DO_NLOGN )
            WRITE(LUPRI,*) "Running hierarchical O(NlogN) algorithm."
         CASE( DO_FMM )
            WRITE(LUPRI,*) "Running hierarchical O(N) FMM algorithm."
         CASE DEFAULT
            CALL fmm_quit ('far-field algorithm type not recognised!')
      END SELECT

      IF ( scheme%include_near_field ) &
            WRITE(LUPRI,*) "Including all classical NF interactions."

      WRITE(LUPRI,'(A,I4)') " LMAX   =", scheme%raw_LMAX
      IF (scheme%algorithm /= DO_FQ) &
         WRITE(LUPRI,'(A,I4)') " TLMAX  =", scheme%trans_LMAX
      IF (scheme%branch_free) THEN
         WRITE(LUPRI,*) "Running in branch-free mode."
         WRITE(LUPRI,'(A,F8.4)')" Minimum extent =", scheme%extent_min
      END IF
      WRITE(LUPRI,'(A,F8.4)')" Smallest box dimension =", scheme%grain

!      IF (scheme%job_type == MD4_FMM ) THEN
!         WRITE(LUPRI,'(A,E12.4)')" Short-range threshold =", EXP(-fmm_ThrInt)
!         WRITE(LUPRI,'(A,F9.4)')" Extent X0 parameter   =", fmm_X0
!      END IF
!      IF (scheme%job_type == FE_FMM ) THEN
!         WRITE(LUPRI,'(A,F8.4)')" Width of finite element =",    &
!                                  scheme%grain/(scheme%FEdim-1)
!      END IF

      IF (scheme%include_near_field) THEN
         SELECT CASE( scheme%T_con%NF_T_buffer )
            CASE( TREE_T_BUFFER )
               WRITE(LUPRI,*) "Using Tree Buffer for NF T matrices."
            CASE( NULL_T_BUFFER )
               WRITE(LUPRI,*) "Building all NF T matrices on the fly."
            CASE( SKIP_T_BUFFER )
               WRITE(LUPRI,*) "Skipping all NF T matrix contractions."
            CASE( MULTI_T_BUFFER )
               WRITE(LUPRI,*) "Using buffer for multiple NF T matrix build."
            CASE( SCALE_T_BUFFER )
               WRITE(LUPRI,*) "Using buffer for scaled NF T matrix build."
            CASE DEFAULT
               CALL fmm_quit('invalid T-vector buffer in fmm_print_scheme!')
         END SELECT
      END IF

      SELECT CASE( scheme%T_con%FF_T_buffer )
         CASE( TREE_T_BUFFER )
            WRITE(LUPRI,*) "Using Tree Buffer for FF T matrices."
         CASE( NULL_T_BUFFER )
            WRITE(LUPRI,*) "Building all FF T matrices on the fly."
         CASE( SKIP_T_BUFFER )
            WRITE(LUPRI,*) "Skipping all FF T matrix contractions."
         CASE( MULTI_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for multiple FF T matrix build."
         CASE( SCALE_T_BUFFER )
            WRITE(LUPRI,*) "Using buffer for scaled FF T matrix build."
         CASE DEFAULT
            CALL fmm_quit ('invalid T-vector buffer in fmm_print_scheme!')
      END SELECT

      SELECT CASE( scheme%W_con%W_buffer )
         CASE( TREE_W_BUFFER )
            WRITE(LUPRI,*) "Using Tree Buffer for W matrices."
         CASE( NULL_W_BUFFER )
            WRITE(LUPRI,*) "Building all W matrices on the fly."
         CASE( SKIP_W_BUFFER )
            WRITE(LUPRI,*) "Skipping all W matrix contractions."
         CASE DEFAULT
            CALL fmm_quit ('invalid W-vector buffer in fmm_print_scheme!')
      END SELECT

      WRITE(LUPRI,'(/,A,/)') " -----------------------------------------"

   END SUBROUTINE fmm_print_scheme

!-------------------------------------------------------------------------------

END MODULE fmm_scheme_builder
