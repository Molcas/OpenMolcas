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
MODULE fmm_stats

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, fmm_stats_printed, MAX_LEVEL, TOP_LEVEL
   IMPLICIT NONE
   PRIVATE

   ! Public procedures
   PUBLIC :: fmm_init_buffer_stats,      &
             fmm_init_matrix_stats,      &
             fmm_print_stats

   ! Public variables
   !-----------------

   REAL(REALK), POINTER, PUBLIC, SAVE :: stat_tpack_unique,      &
                                         stat_tpack_total,       &
                                         stat_tpack_chunks,      &
                                         stat_T_mat_builds,      &
                                         stat_W_mat_builds

   INTEGER(INTK), PUBLIC, SAVE :: stat_iteration=0,             &
                                  stat_n_basis=0,               &
                                  stat_points=0,                &
                                  stat_raw_moms_RHS=0,          &
                                  stat_pkd_moms_RHS=0,          &
                                  stat_screened_moms_RHS=0,     &
                                  stat_raw_moms_LHS=0,          &
                                  stat_pkd_moms_LHS=0,          &
                                  stat_deepest_level=0,         &
                                  stat_max_branch=0,            &
                                  stat_min_branch=10000,        &
                                  stat_level_saturation=0,      &
                                  stat_RHS_boxes(MAX_LEVEL)=0,  &
                                  stat_LHS_boxes(MAX_LEVEL)=0

   LOGICAL, PUBLIC, SAVE :: stat_NF_not_FF=.FALSE.

   ! Private variables
   !------------------

   REAL(REALK), TARGET, SAVE :: stat_T_direction_NF=0,    &
                                stat_T_matrix_NF=0,       &
                                stat_T_chunks_NF=0,       &
                                stat_T_total_NF=0,        &
                                stat_T_direction_FF=0,    &
                                stat_T_matrix_FF=0,       &
                                stat_T_chunks_FF=0,       &
                                stat_T_total_FF=0

   REAL(REALK), TARGET, SAVE :: stat_W_direction_RB=0,    &
                                stat_W_matrix_RB=0,       &
                                stat_W_total_RB=0,        &
                                stat_W_chunks_RB=0,       &
                                stat_W_direction_BB=0,    &
                                stat_W_matrix_BB=0,       &
                                stat_W_total_BB=0,        &
                                stat_W_chunks_BB=0,       &
                                stat_W_direction_BR=0,    &
                                stat_W_matrix_BR=0,       &
                                stat_W_total_BR=0,        &
                                stat_W_chunks_BR=0

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_buffer_stats(mode,W_mode)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN)           :: mode
      CHARACTER(LEN=7), OPTIONAL, INTENT(IN) :: W_mode

      SELECT CASE (mode)
         CASE ('T')
            IF (stat_NF_not_FF) THEN
               stat_tpack_chunks => stat_T_chunks_NF
               stat_tpack_unique => stat_T_direction_NF
               stat_tpack_total => stat_T_total_NF
            ELSE
               stat_tpack_chunks => stat_T_chunks_FF
               stat_tpack_unique => stat_T_direction_FF
               stat_tpack_total => stat_T_total_FF
            END IF
         CASE ('W')
            SELECT CASE (W_mode)
               CASE ('RAW_BOX')
                  stat_tpack_chunks => stat_W_chunks_RB
                  stat_tpack_unique => stat_W_direction_RB
                  stat_tpack_total => stat_W_total_RB
               CASE ('BOX_BOX')
                  stat_tpack_chunks => stat_W_chunks_BB
                  stat_tpack_unique => stat_W_direction_BB
                  stat_tpack_total => stat_W_total_BB
               CASE ('BOX_RAW')
                  stat_tpack_chunks => stat_W_chunks_BR
                  stat_tpack_unique => stat_W_direction_BR
                  stat_tpack_total => stat_W_total_BR
               CASE DEFAULT
                  CALL fmm_quit('cannot reconcile W runtype!')
            END SELECT
         CASE DEFAULT
            CALL fmm_quit('cannot reconcile buffer statistics requested')
      END SELECT

   END SUBROUTINE fmm_init_buffer_stats

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_matrix_stats(mode,W_mode)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN)           :: mode
      CHARACTER(LEN=7), OPTIONAL, INTENT(IN) :: W_mode

      SELECT CASE (mode)
         CASE ('T')
            IF (stat_NF_not_FF) THEN
               stat_T_mat_builds => stat_T_matrix_NF
            ELSE
               stat_T_mat_builds => stat_T_matrix_FF
            END IF
         CASE ('W')
            SELECT CASE (W_mode)
               CASE ('RAW_BOX')
                  stat_W_mat_builds => stat_W_matrix_RB
               CASE ('BOX_BOX')
                  stat_W_mat_builds => stat_W_matrix_BB
               CASE ('BOX_RAW')
                  stat_W_mat_builds => stat_W_matrix_BR
               CASE DEFAULT
                  CALL fmm_quit('cannot reconcile W runtype!')
            END SELECT
         CASE DEFAULT
            CALL fmm_quit('cannot reconcile buffer statistics requested')
      END SELECT

   END SUBROUTINE fmm_init_matrix_stats

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_print_stats

      IMPLICIT NONE
      INTEGER(INTK) :: i, RHS_tot, LHS_tot

      IF (.NOT. fmm_stats_printed) THEN
         RETURN
      END IF

      WRITE(LUPRI,'(/,A)') " Multipole method statistics"
      WRITE(LUPRI,'(A,/)') " ---------------------------"

      WRITE(LUPRI,'(A,I8)') " Number of AOs   =", stat_n_basis
      WRITE(LUPRI,'(A,I6)') " Deepest level   =", stat_deepest_level
      IF (stat_level_saturation > 0) &
         WRITE(LUPRI,*) 'deepest level saturated!', stat_level_saturation

      WRITE(LUPRI,'(A,I6)') " Max branch      =", stat_max_branch
      WRITE(LUPRI,'(A,I6)') " Min branch      =", stat_min_branch
      WRITE(LUPRI,'(/,A)') " Unboxed moments:"
      WRITE(LUPRI,'(A,7X,A,I9)') " RHS raw","=", stat_raw_moms_RHS
      WRITE(LUPRI,'(A,4X,A,I9)') " RHS packed","=", stat_pkd_moms_RHS
      WRITE(LUPRI,'(A,2X,A,I9)') " RHS screened","=", stat_screened_moms_RHS
      WRITE(LUPRI,'(A,7X,A,I9)') " LHS raw","=", stat_raw_moms_LHS
      WRITE(LUPRI,'(A,4X,A,I9)') " LHS packed","=", stat_pkd_moms_LHS

      LHS_tot = 0
      RHS_tot = 0
      DO i = stat_deepest_level, TOP_LEVEL, -1
         IF (stat_RHS_boxes(i) /= -1) THEN
            WRITE(LUPRI,'(A,I2)') " Boxes at level: ", i
            WRITE(LUPRI,'(A,5X,A,I9)') " RHS boxes","=", stat_RHS_boxes(i)
            RHS_tot = RHS_tot + stat_RHS_boxes(i)
         END IF
         IF (stat_LHS_boxes(i) /= -1) THEN
            WRITE(LUPRI,'(A,5X,A,I9)') " LHS boxes","=", stat_LHS_boxes(i)
            LHS_tot = LHS_tot + stat_LHS_boxes(i)
         END IF
      END DO
      WRITE(LUPRI,'(/,A,2X,A,I9)') " Total RHS boxes","=", RHS_tot
      WRITE(LUPRI,'(A,2X,A,I9)') " Total LHS boxes","=", LHS_tot

      IF (stat_T_total_NF > 0) THEN
         WRITE(LUPRI,'(/,A)') " NF contraction pairs"
         WRITE(LUPRI,'(A,E17.9)') " total interactions     =", &
                                    stat_T_total_NF
         WRITE(LUPRI,'(A,E17.9)') " total directions       =", &
                                    stat_T_direction_NF
         WRITE(LUPRI,'(A,E17.9)') " total buffer chunks    =", &
                                    stat_T_chunks_NF
         WRITE(LUPRI,'(A,E17.9)') " T matrices built       =", &
                                    stat_T_matrix_NF
      END IF
      WRITE(LUPRI,'(/,A)') " FF contraction pairs"
      WRITE(LUPRI,'(A,E17.9)') " total interactions     =", stat_T_total_FF
      WRITE(LUPRI,'(A,E17.9)') " total directions       =", stat_T_direction_FF
      WRITE(LUPRI,'(A,E17.9)') " total buffer chunks    =", stat_T_chunks_FF
      WRITE(LUPRI,'(A,E17.9)') " T matrices built       =", stat_T_matrix_FF

      WRITE(LUPRI,'(/,A)') " raw-box translations"
      WRITE(LUPRI,'(A,E17.9)') " total number         =", stat_W_total_RB
      WRITE(LUPRI,'(A,E17.9)') " total directions     =", stat_W_direction_RB
      WRITE(LUPRI,'(A,E17.9)') " total buffer chunks  =", stat_W_chunks_RB
      WRITE(LUPRI,'(A,E17.9)') " W matrices built     =", stat_W_matrix_RB

      WRITE(LUPRI,'(/,A)') " box-box translations"
      WRITE(LUPRI,'(A,E17.9)') " total number         =", stat_W_total_BB
      WRITE(LUPRI,'(A,E17.9)') " total directions     =", stat_W_direction_BB
      WRITE(LUPRI,'(A,E17.9)') " total buffer chunks  =", stat_W_chunks_BB
      WRITE(LUPRI,'(A,E17.9)') " W matrices built     =", stat_W_matrix_BB

      WRITE(LUPRI,'(/,A)') " box-raw translations"
      WRITE(LUPRI,'(A,E17.9)') " total number         =", stat_W_total_BR
      WRITE(LUPRI,'(A,E17.9)') " total directions     =", stat_W_direction_BR
      WRITE(LUPRI,'(A,E17.9)') " total buffer chunks  =", stat_W_chunks_BR
      WRITE(LUPRI,'(A,E17.9)') " W matrices built     =", stat_W_matrix_BR

      WRITE(LUPRI,'(/,A,/)') " ------------------------------------------"

   END SUBROUTINE fmm_print_stats

!-------------------------------------------------------------------------------

END MODULE fmm_stats
