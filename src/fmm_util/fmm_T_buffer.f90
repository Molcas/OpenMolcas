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
MODULE fmm_T_buffer

   USE fmm_global_paras
   USE fmm_stats

   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: fmm_add_to_T_buffer,   &
             fmm_open_T_buffer,     &
             fmm_close_T_buffer

   ! only one T_buffer can be open at once
   INTEGER(INTK), SAVE :: buffer = -1

   ! diagnostic flag
   CHARACTER(4), SAVE :: T_buffer_stat

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_add_to_T_buffer(T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL fmm_selected_t_buffer
      EXTERNAL fmm_selected_t_contractor

      CALL fmm_selected_t_buffer(fmm_selected_t_contractor,T_pair)

   END SUBROUTINE fmm_add_to_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_close_T_buffer(scheme)

      USE fmm_T_contractors,  ONLY: fmm_lock_T_con
      USE fmm_tree_buffer,    ONLY: fmm_tree_buffer_finish
      USE fmm_multi_T_buffer, ONLY: fmm_free_multi_T_buffer
      USE fmm_scale_T_buffer, ONLY: fmm_free_scale_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      EXTERNAL fmm_selected_t_contractor

      IF (T_buffer_stat /= 'OPEN') CALL fmm_quit ('T_buffer already closed!')

      SELECT CASE (buffer)
         CASE (SKIP_T_BUFFER)
            ! do nothing
         CASE (NULL_T_BUFFER)
            ! do nothing
         CASE (TREE_T_BUFFER)
            CALL fmm_tree_buffer_finish(fmm_selected_t_contractor)
         CASE (MULTI_T_BUFFER)
            CALL fmm_free_multi_T_buffer(fmm_selected_t_contractor)
         CASE (SCALE_T_BUFFER)
            CALL fmm_free_scale_T_buffer(fmm_selected_t_contractor)
         CASE DEFAULT
            CALL fmm_quit ('cannot reconcile list type in fmm_close_T_buffer')
      END SELECT

      T_buffer_stat = 'FREE'
      fmm_lock_T_con = .FALSE.

! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real(scheme) ! not really real, but well...
   END SUBROUTINE fmm_close_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_null_T_buffer(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor

      stat_tpack_total = stat_tpack_total + one
      CALL T_contractor(T_pair)

   END SUBROUTINE fmm_null_T_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

   SUBROUTINE fmm_skip_T_buffer(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor

      stat_T_mat_builds = stat_T_mat_builds + 1
      RETURN

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_real(T_contractor)
         CALL Unused_real(T_pair) ! not really real, but well...
      END IF
   END SUBROUTINE fmm_skip_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_open_T_buffer(scheme)

      USE fmm_T_contractors,  ONLY: fmm_lock_T_con
      USE fmm_tree_buffer,    ONLY: fmm_tree_buffer_init,        &
                                    fmm_tree_buffer_add
      USE fmm_multi_T_buffer, ONLY: fmm_init_multi_T_buffer,     &
                                    fmm_multi_T_buffer_add
      USE fmm_scale_T_buffer, ONLY: fmm_init_scale_T_buffer,     &
                                    fmm_scale_T_buffer_add

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK) :: sort_para
      EXTERNAL fmm_store_t_buffer

      CALL fmm_init_buffer_stats('T')
      IF (T_buffer_stat == 'OPEN') CALL fmm_quit('cannot reopen T_buffer')

      IF (scheme%phase == NEAR_FIELD) THEN
         buffer = scheme%T_con%NF_T_buffer
         sort_para = scheme%T_con%NF_sort_para
      ELSE
         buffer = scheme%T_con%FF_T_buffer
         sort_para = scheme%T_con%FF_sort_para
      END IF

      SELECT CASE (buffer)
         CASE (SKIP_T_BUFFER)
            ! all T-contractions will be skipped by this choice of buffer
            CALL fmm_store_t_buffer(fmm_skip_T_buffer)
         CASE (NULL_T_BUFFER)
            CALL fmm_store_t_buffer(fmm_null_T_buffer)
         CASE (TREE_T_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL fmm_store_t_buffer(fmm_tree_buffer_add)
            CALL fmm_tree_buffer_init(TREE_LENGTH,sort_para)
         CASE (SCALE_T_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL fmm_store_t_buffer(fmm_scale_T_buffer_add)
            CALL fmm_init_scale_T_buffer
         CASE (MULTI_T_BUFFER)
             ! use buffer to drive multiple T matrix simultaneous build
            CALL fmm_store_t_buffer(fmm_multi_T_buffer_add)
            CALL fmm_init_multi_T_buffer(TMATM_DF)
         CASE DEFAULT
            CALL fmm_quit ('cannot reconcile list type in fmm_open_T_buffer')
      END SELECT

      T_buffer_stat = 'OPEN'
      fmm_lock_T_con = .TRUE.

   END SUBROUTINE fmm_open_T_buffer

!-------------------------------------------------------------------------------

END MODULE fmm_T_buffer
