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

module fmm_scheme_builder

use fmm_global_paras, only: INTK, REALK, LUPRI, LURD, scheme_paras, GFC_FMM, MD4_FMM, FE_FMM, DO_NULL, DO_FQ, DO_NlogN, DO_BQ, &
                            DO_FMM, TREE_T_BUFFER, NULL_T_BUFFER, SKIP_T_BUFFER, MULTI_T_BUFFER, SCALE_T_BUFFER, TREE_W_BUFFER, &
                            NULL_W_BUFFER, SKIP_W_BUFFER, WS_MIN, BRFREE_DF, EXTENT_MIN_DF, PACK_LHS_DF, PACK_RHS_DF, &
                            SORT_BY_SCALE, T_CONTRACTOR_BOUNDARY, T_CONTRACTOR_DIRECT, T_CONTRACTOR_FULL, T_CONTRACTOR_SCALE, &
                            USE_RAW_QLM, USE_T_SYM_QLM, W_CONTRACTOR_BOUNDARY, W_CONTRACTOR_FAST, One, Half
use fmm_stats, only: stat_iteration
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_init_scheme, &
          fmm_get_scheme

! "scheme" contains all the input (and default) parameter data to
! completely define a unique MM run-type
type(scheme_paras), target, save :: scheme

! Safety check for initialisation
logical, save :: scheme_initialised = .false.

contains

!-------------------------------------------------------------------------------

subroutine fmm_get_scheme(scheme_ptr)

  implicit none
  type(scheme_paras), pointer :: scheme_ptr
  integer(INTK) :: iteration = 0

  if (.not. scheme_initialised) call fmm_quit('fmm scheme uninitialised!')

  nullify(scheme_ptr)
  scheme_ptr => scheme

  iteration = iteration+1
  stat_iteration = iteration

end subroutine fmm_get_scheme

!-------------------------------------------------------------------------------

subroutine fmm_init_scheme(job_type)

  !use fmm_md4_globals, only: fmm_grain_inv, fmm_extent_min

  implicit none
  integer(INTK), intent(in) :: job_type

  integer(INTK) :: flag
  integer(INTK) :: LMAX, TLMAX, ALGORITHM, FEdim, LIPN
  real(REALK) :: GRAIN, DENS_SCREEN, EXTENT_MIN

  integer(INTK) :: ERROR
  namelist /FMM/ LMAX,TLMAX,ALGORITHM,GRAIN,DENS_SCREEN,EXTENT_MIN,FEdim,LIPN

  ! Set default scheme

  scheme%job_type = job_type
  scheme%branch_free = BRFREE_DF
  scheme%pack_LHS = PACK_LHS_DF
  scheme%pack_RHS = PACK_RHS_DF
  ! We assume a symmetric T-matrix in all cases so
  ! use auxiliary (prefactor-adapted) moment types for RHS.
  scheme%T_con%LHS_mm_type = USE_RAW_QLM
  scheme%T_con%RHS_mm_type = USE_T_SYM_QLM

  ! Change defaults with user input

  !fixme Put defaults in parameter module
  ALGORITHM = DO_FMM
  LMAX = 4
  TLMAX = 12
  GRAIN = One
  DENS_SCREEN = 1.0e-15_REALK
  EXTENT_MIN = EXTENT_MIN_DF
  FEdim = 10
  LIPN = 2
  rewind(LURD)
  read(LURD,FMM,iostat=ERROR)
  if (ERROR > 0) then
    write(LUPRI,*) 'o Check NAMELIST FMM'
    call Abend()
  end if
  scheme%algorithm = ALGORITHM
  scheme%raw_LMAX = LMAX
  scheme%trans_LMAX = TLMAX
  scheme%grain = GRAIN
  !fmm_grain_inv = one/scheme%grain
  scheme%dens_screen_thr = DENS_SCREEN
  scheme%extent_min = EXTENT_MIN
  scheme%FEdim = FEdim
  scheme%lipn = LIPN

  ! Set up (W) translator and (T) contractor options

  select case(scheme%job_type)
    case(GFC_FMM)
      scheme%include_near_field = .true.
      scheme%W_con%W_buffer = TREE_W_BUFFER
      scheme%W_con%ID = W_CONTRACTOR_FAST
      scheme%W_con%BR_ID = W_CONTRACTOR_BOUNDARY
      scheme%W_con%sort_para = SORT_BY_SCALE
      scheme%T_con%NF_ID = T_CONTRACTOR_BOUNDARY
      scheme%T_con%NF_T_buffer = NULL_T_BUFFER
      if (scheme%algorithm == DO_FQ) then
        scheme%T_con%FF_ID = T_CONTRACTOR_BOUNDARY
        scheme%T_con%FF_T_buffer = NULL_T_BUFFER
      else
        scheme%T_con%FF_T_buffer = SCALE_T_BUFFER
        scheme%T_con%FF_ID = T_CONTRACTOR_SCALE
      end if
    case(MD4_FMM,FE_FMM)
      scheme%include_near_field = .false.
      scheme%W_con%W_buffer = TREE_W_BUFFER
      scheme%W_con%ID = W_CONTRACTOR_FAST
      scheme%W_con%BR_ID = W_CONTRACTOR_FAST
      scheme%W_con%sort_para = SORT_BY_SCALE
      scheme%T_con%NF_T_buffer = NULL_T_BUFFER
      scheme%T_con%NF_ID = T_CONTRACTOR_FULL
      flag = 2
      if (scheme%algorithm == DO_FQ) flag = 1
      select case(flag)
        case(0)
          ! we do no contractions (for diagnostics)
          scheme%T_con%FF_T_buffer = SKIP_T_BUFFER
          scheme%T_con%FF_ID = T_CONTRACTOR_DIRECT
        case(1)
          scheme%T_con%FF_T_buffer = NULL_T_BUFFER
          scheme%T_con%FF_ID = T_CONTRACTOR_FULL
        case(2)
          scheme%T_con%FF_T_buffer = SCALE_T_BUFFER
          scheme%T_con%FF_ID = T_CONTRACTOR_SCALE
          !scheme%T_con%sort_para = SORT_BY_SCALE
        !case(3)
        !  scheme%T_con%T_buffer = TREE_T_BUFFER
        !  scheme%T_con%sort_para = SORT_BY_SCALE
        !  scheme%T_con%ID = T_CONTRACTOR_TREE
        !case(4)
        !  scheme%T_con%T_buffer = TREE_T_BUFFER
        !  scheme%T_con%ID = T_CONTRACTOR_SCALE_TREE
        case default
          call fmm_quit('invalid T-contractor specified!')
      end select
    case default
      call fmm_quit('invalid FMM run-type requested!')
  end select

  call fmm_verify_scheme()
  call fmm_print_scheme()
  scheme_initialised = .true.

end subroutine fmm_init_scheme

!-------------------------------------------------------------------------------

subroutine fmm_verify_scheme()

  implicit none

  select case(scheme%algorithm)
    case(MD4_FMM)
      if (WS_MIN < 2*ceiling(scheme%extent_min/scheme%grain*half)) then
        write(LUPRI,*) 'WS_MIN = ',WS_MIN
        write(LUPRI,*) 'Extent_min = ',scheme%extent_min
        write(LUPRI,*) 'Grain  = ',scheme%grain
        call fmm_quit('RPQ cut off too large or boxes too small!')
      end if
    case default
      continue
  end select

  if (scheme%trans_LMAX < scheme%raw_LMAX) call fmm_quit('increase TLMAX!')

end subroutine fmm_verify_scheme

!-------------------------------------------------------------------------------

subroutine fmm_print_scheme()

  !use fmm_md4_globals

  implicit none

  write(LUPRI,'(/,A)') ' -----------------------------------------'
  write(LUPRI,'(A)') ' |  Multipole module runtime parameters  |'
  write(LUPRI,'(A,/)') ' -----------------------------------------'

  select case(scheme%job_type)
    case(GFC_FMM)
      write(LUPRI,*) 'Computing classical boundary potential.'
    case(MD4_FMM)
      write(LUPRI,*) 'Computing multipole contribution to J-matrix.'
    case(FE_FMM)
      write(LUPRI,*) 'Computing full J-matrix via FE-FMM.'
    case default
      call fmm_quit('MM job type not recognised!')
  end select

  select case(scheme%algorithm)
    case(DO_NULL)
      write(LUPRI,*) 'Skipping all FF interactions.'
    case(DO_FQ)
      write(LUPRI,*) 'Running simple O(N^2) algorithm.'
    case(DO_BQ)
      write(LUPRI,*) 'Running fast O(N^2) algorithm with boxes.'
    case(DO_NLOGN)
      write(LUPRI,*) 'Running hierarchical O(NlogN) algorithm.'
    case(DO_FMM)
      write(LUPRI,*) 'Running hierarchical O(N) FMM algorithm.'
    case default
      call fmm_quit('far-field algorithm type not recognised!')
  end select

  if (scheme%include_near_field) write(LUPRI,*) 'Including all classical NF interactions.'

  write(LUPRI,'(A,I4)') ' LMAX   =',scheme%raw_LMAX
  if (scheme%algorithm /= DO_FQ) write(LUPRI,'(A,I4)') ' TLMAX  =',scheme%trans_LMAX
  if (scheme%branch_free) then
    write(LUPRI,*) 'Running in branch-free mode.'
    write(LUPRI,'(A,F8.4)') ' Minimum extent =',scheme%extent_min
  end if
  write(LUPRI,'(A,F8.4)') ' Smallest box dimension =',scheme%grain

  !if (scheme%job_type == MD4_FMM) then
  !   write(LUPRI,'(A,E12.4)') 'Short-range threshold =',exp(-fmm_ThrInt)
  !   write(LUPRI,'(A,F9.4)') 'Extent X0 parameter   =',fmm_X0
  !end if
  !if (scheme%job_type == FE_FMM) then
  !   write(LUPRI,'(A,F8.4)') 'Width of finite element =',scheme%grain/(scheme%FEdim-1)
  !end if

  if (scheme%include_near_field) then
    select case(scheme%T_con%NF_T_buffer)
    case(TREE_T_BUFFER)
      write(LUPRI,*) 'Using Tree Buffer for NF T matrices.'
    case(NULL_T_BUFFER)
      write(LUPRI,*) 'Building all NF T matrices on the fly.'
    case(SKIP_T_BUFFER)
      write(LUPRI,*) 'Skipping all NF T matrix contractions.'
    case(MULTI_T_BUFFER)
      write(LUPRI,*) 'Using buffer for multiple NF T matrix build.'
    case(SCALE_T_BUFFER)
      write(LUPRI,*) 'Using buffer for scaled NF T matrix build.'
    case default
      call fmm_quit('invalid T-vector buffer in fmm_print_scheme!')
    end select
  end if

  select case(scheme%T_con%FF_T_buffer)
    case(TREE_T_BUFFER)
      write(LUPRI,*) 'Using Tree Buffer for FF T matrices.'
    case(NULL_T_BUFFER)
      write(LUPRI,*) 'Building all FF T matrices on the fly.'
    case(SKIP_T_BUFFER)
      write(LUPRI,*) 'Skipping all FF T matrix contractions.'
    case(MULTI_T_BUFFER)
      write(LUPRI,*) 'Using buffer for multiple FF T matrix build.'
    case(SCALE_T_BUFFER)
      write(LUPRI,*) 'Using buffer for scaled FF T matrix build.'
    case default
      call fmm_quit('invalid T-vector buffer in fmm_print_scheme!')
  end select

  select case(scheme%W_con%W_buffer)
    case(TREE_W_BUFFER)
      write(LUPRI,*) 'Using Tree Buffer for W matrices.'
    case(NULL_W_BUFFER)
      write(LUPRI,*) 'Building all W matrices on the fly.'
    case(SKIP_W_BUFFER)
      write(LUPRI,*) 'Skipping all W matrix contractions.'
    case default
      call fmm_quit('invalid W-vector buffer in fmm_print_scheme!')
  end select

  write(LUPRI,'(/,A,/)') ' -----------------------------------------'

end subroutine fmm_print_scheme

!-------------------------------------------------------------------------------

end module fmm_scheme_builder
