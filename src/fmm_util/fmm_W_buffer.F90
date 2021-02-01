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
MODULE fmm_W_buffer

   USE fmm_global_paras

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_add_to_W_buffer,   &
             fmm_open_W_buffer,     &
             fmm_close_W_buffer

   ! diagnostic flag
   CHARACTER(LEN=4), SAVE :: W_buffer_stat

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_add_to_W_buffer(W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL fmm_selected_w_buffer
      EXTERNAL fmm_selected_w_contractor

      CALL fmm_selected_w_buffer(fmm_selected_w_contractor,W_pair)

   END SUBROUTINE fmm_add_to_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_close_W_buffer(scheme)

      USE fmm_W_contractors, ONLY: fmm_lock_W_con
      USE fmm_tree_buffer,   ONLY: fmm_tree_buffer_finish

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      EXTERNAL fmm_selected_w_contractor

      IF (W_buffer_stat /= 'OPEN') CALL fmm_quit ('W_buffer already closed!')
      SELECT CASE (scheme%W_con%W_buffer)
         CASE (SKIP_W_BUFFER)
            ! do nothing
         CASE (NULL_W_BUFFER)
            ! do nothing
         CASE (TREE_W_BUFFER)
            CALL fmm_tree_buffer_finish(fmm_selected_w_contractor)
         CASE DEFAULT
            CALL fmm_quit ('cannot reconcile list type in fmm_close_W_buffer')
      END SELECT
      W_buffer_stat = 'FREE'
      fmm_lock_W_con = .FALSE.

   END SUBROUTINE fmm_close_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_null_W_buffer(W_contractor,W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL W_contractor

      CALL W_contractor(W_pair)

   END SUBROUTINE fmm_null_W_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

   SUBROUTINE fmm_skip_W_buffer(W_contractor,W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL W_contractor

      RETURN

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL W_contractor(W_pair)
      END IF
   END SUBROUTINE fmm_skip_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_open_W_buffer(scheme)

      USE fmm_W_contractors, ONLY: fmm_lock_W_con
      USE fmm_tree_buffer,   ONLY: fmm_tree_buffer_init,      &
                                   fmm_tree_buffer_add

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      EXTERNAL fmm_store_w_buffer

      IF (W_buffer_stat == 'OPEN') CALL fmm_quit('cannot reopen W_buffer')

      SELECT CASE (scheme%W_con%W_buffer)
         CASE (SKIP_W_BUFFER)
            ! all W-contractions will be skipped by this choice of buffer
            CALL fmm_store_w_buffer(fmm_skip_W_buffer)
         CASE (NULL_W_BUFFER)
            CALL fmm_store_w_buffer(fmm_null_W_buffer)
         CASE (TREE_W_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL fmm_store_w_buffer(fmm_tree_buffer_add)
            CALL fmm_tree_buffer_init(TREE_LENGTH,scheme%W_con%sort_para)
         CASE DEFAULT
         CALL fmm_quit ('cannot reconcile list type in fmm_open_W_buffer')
      END SELECT

      W_buffer_stat = 'OPEN'
      fmm_lock_W_con = .TRUE.

   END SUBROUTINE fmm_open_W_buffer

!-------------------------------------------------------------------------------

END MODULE fmm_W_buffer
