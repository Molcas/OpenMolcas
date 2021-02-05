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

module fmm_stats

use fmm_global_paras, only: INTK, REALK, LUPRI, fmm_stats_printed, MAX_LEVEL, TOP_LEVEL
use fmm_utils, only: fmm_quit
implicit none
private

! Public procedures
public :: fmm_init_buffer_stats, &
          fmm_init_matrix_stats, &
          fmm_print_stats

! Public variables
!-----------------

real(REALK), pointer, public, save :: stat_tpack_unique, &
                                      stat_tpack_total, &
                                      stat_tpack_chunks, &
                                      stat_T_mat_builds, &
                                      stat_W_mat_builds

integer(INTK), public, save :: stat_iteration = 0, &
                               stat_n_basis = 0, &
                               stat_points = 0, &
                               stat_raw_moms_RHS = 0, &
                               stat_pkd_moms_RHS = 0, &
                               stat_screened_moms_RHS = 0, &
                               stat_raw_moms_LHS = 0, &
                               stat_pkd_moms_LHS = 0, &
                               stat_deepest_level = 0, &
                               stat_max_branch = 0, &
                               stat_min_branch = 10000, &
                               stat_level_saturation = 0, &
                               stat_RHS_boxes(MAX_LEVEL) = 0, &
                               stat_LHS_boxes(MAX_LEVEL) = 0

logical, public, save :: stat_NF_not_FF = .false.

! Private variables
!------------------

real(REALK), target, save :: stat_T_direction_NF = 0, &
                             stat_T_matrix_NF = 0, &
                             stat_T_chunks_NF = 0, &
                             stat_T_total_NF = 0, &
                             stat_T_direction_FF = 0, &
                             stat_T_matrix_FF = 0, &
                             stat_T_chunks_FF = 0, &
                             stat_T_total_FF = 0

real(REALK), target, save :: stat_W_direction_RB = 0, &
                             stat_W_matrix_RB = 0, &
                             stat_W_total_RB = 0, &
                             stat_W_chunks_RB = 0, &
                             stat_W_direction_BB = 0, &
                             stat_W_matrix_BB = 0, &
                             stat_W_total_BB = 0, &
                             stat_W_chunks_BB = 0, &
                             stat_W_direction_BR = 0, &
                             stat_W_matrix_BR = 0, &
                             stat_W_total_BR = 0, &
                             stat_W_chunks_BR = 0

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_buffer_stats(mode,W_mode)

  implicit none
  character, intent(in)                  :: mode
  character(len=7), optional, intent(in) :: W_mode

  select case (mode)
    case ('T')
      if (stat_NF_not_FF) then
        stat_tpack_chunks => stat_T_chunks_NF
        stat_tpack_unique => stat_T_direction_NF
        stat_tpack_total => stat_T_total_NF
      else
        stat_tpack_chunks => stat_T_chunks_FF
        stat_tpack_unique => stat_T_direction_FF
        stat_tpack_total => stat_T_total_FF
      end if
    case ('W')
      select case (W_mode)
        case ('RAW_BOX')
          stat_tpack_chunks => stat_W_chunks_RB
          stat_tpack_unique => stat_W_direction_RB
          stat_tpack_total => stat_W_total_RB
        case ('BOX_BOX')
          stat_tpack_chunks => stat_W_chunks_BB
          stat_tpack_unique => stat_W_direction_BB
          stat_tpack_total => stat_W_total_BB
        case ('BOX_RAW')
          stat_tpack_chunks => stat_W_chunks_BR
          stat_tpack_unique => stat_W_direction_BR
          stat_tpack_total => stat_W_total_BR
        case default
          call fmm_quit('cannot reconcile W runtype!')
      end select
    case default
      call fmm_quit('cannot reconcile buffer statistics requested')
  end select

end subroutine fmm_init_buffer_stats

!-------------------------------------------------------------------------------

subroutine fmm_init_matrix_stats(mode,W_mode)

  implicit none
  character, intent(in)                  :: mode
  character(len=7), optional, intent(in) :: W_mode

  select case (mode)
    case ('T')
      if (stat_NF_not_FF) then
        stat_T_mat_builds => stat_T_matrix_NF
      else
        stat_T_mat_builds => stat_T_matrix_FF
      end if
    case ('W')
      select case (W_mode)
        case ('RAW_BOX')
          stat_W_mat_builds => stat_W_matrix_RB
        case ('BOX_BOX')
          stat_W_mat_builds => stat_W_matrix_BB
        case ('BOX_RAW')
          stat_W_mat_builds => stat_W_matrix_BR
        case default
          call fmm_quit('cannot reconcile W runtype!')
      end select
    case default
      call fmm_quit('cannot reconcile buffer statistics requested')
  end select

end subroutine fmm_init_matrix_stats

!-------------------------------------------------------------------------------

subroutine fmm_print_stats()

  implicit none
  integer(INTK) :: i, RHS_tot, LHS_tot

  if (.not. fmm_stats_printed) then
    return
  end if

  write(LUPRI,'(/,A)') ' Multipole method statistics'
  write(LUPRI,'(A,/)') ' ---------------------------'

  write(LUPRI,'(A,I8)') ' Number of AOs   =',stat_n_basis
  write(LUPRI,'(A,I6)') ' Deepest level   =',stat_deepest_level
  if (stat_level_saturation > 0) write(LUPRI,*) 'deepest level saturated!',stat_level_saturation

  write(LUPRI,'(A,I6)') ' Max branch      =',stat_max_branch
  write(LUPRI,'(A,I6)') ' Min branch      =',stat_min_branch
  write(LUPRI,'(/,A)') ' Unboxed moments:'
  write(LUPRI,'(A,7X,A,I9)') ' RHS raw','=',stat_raw_moms_RHS
  write(LUPRI,'(A,4X,A,I9)') ' RHS packed','=',stat_pkd_moms_RHS
  write(LUPRI,'(A,2X,A,I9)') ' RHS screened','=',stat_screened_moms_RHS
  write(LUPRI,'(A,7X,A,I9)') ' LHS raw','=',stat_raw_moms_LHS
  write(LUPRI,'(A,4X,A,I9)') ' LHS packed','=',stat_pkd_moms_LHS

  LHS_tot = 0
  RHS_tot = 0
  do i=stat_deepest_level,TOP_LEVEL,-1
    if (stat_RHS_boxes(i) /= -1) then
      write(LUPRI,'(A,I2)') ' Boxes at level: ',i
      write(LUPRI,'(A,5X,A,I9)') ' RHS boxes','=',stat_RHS_boxes(i)
      RHS_tot = RHS_tot+stat_RHS_boxes(i)
    end if
    if (stat_LHS_boxes(i) /= -1) then
      write(LUPRI,'(A,5X,A,I9)') ' LHS boxes','=',stat_LHS_boxes(i)
      LHS_tot = LHS_tot+stat_LHS_boxes(i)
    end if
  end do
  write(LUPRI,'(/,A,2X,A,I9)') ' Total RHS boxes','=',RHS_tot
  write(LUPRI,'(A,2X,A,I9)') ' Total LHS boxes','=',LHS_tot

  if (stat_T_total_NF > 0) then
    write(LUPRI,'(/,A)') ' NF contraction pairs'
    write(LUPRI,'(A,E17.9)') ' total interactions     =',stat_T_total_NF
    write(LUPRI,'(A,E17.9)') ' total directions       =',stat_T_direction_NF
    write(LUPRI,'(A,E17.9)') ' total buffer chunks    =',stat_T_chunks_NF
    write(LUPRI,'(A,E17.9)') ' T matrices built       =',stat_T_matrix_NF
  end if
  write(LUPRI,'(/,A)') ' FF contraction pairs'
  write(LUPRI,'(A,E17.9)') ' total interactions     =',stat_T_total_FF
  write(LUPRI,'(A,E17.9)') ' total directions       =',stat_T_direction_FF
  write(LUPRI,'(A,E17.9)') ' total buffer chunks    =',stat_T_chunks_FF
  write(LUPRI,'(A,E17.9)') ' T matrices built       =',stat_T_matrix_FF

  write(LUPRI,'(/,A)') ' raw-box translations'
  write(LUPRI,'(A,E17.9)') ' total number         =',stat_W_total_RB
  write(LUPRI,'(A,E17.9)') ' total directions     =',stat_W_direction_RB
  write(LUPRI,'(A,E17.9)') ' total buffer chunks  =',stat_W_chunks_RB
  write(LUPRI,'(A,E17.9)') ' W matrices built     =',stat_W_matrix_RB

  write(LUPRI,'(/,A)') ' box-box translations'
  write(LUPRI,'(A,E17.9)') ' total number         =',stat_W_total_BB
  write(LUPRI,'(A,E17.9)') ' total directions     =',stat_W_direction_BB
  write(LUPRI,'(A,E17.9)') ' total buffer chunks  =',stat_W_chunks_BB
  write(LUPRI,'(A,E17.9)') ' W matrices built     =',stat_W_matrix_BB

  write(LUPRI,'(/,A)') ' box-raw translations'
  write(LUPRI,'(A,E17.9)') ' total number         =',stat_W_total_BR
  write(LUPRI,'(A,E17.9)') ' total directions     =',stat_W_direction_BR
  write(LUPRI,'(A,E17.9)') ' total buffer chunks  =',stat_W_chunks_BR
  write(LUPRI,'(A,E17.9)') ' W matrices built     =',stat_W_matrix_BR

  write(LUPRI,'(/,A,/)') ' ------------------------------------------'

end subroutine fmm_print_stats

!-------------------------------------------------------------------------------

end module fmm_stats
