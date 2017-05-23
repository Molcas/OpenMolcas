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
MODULE fmm_W_pair_builder

   USE fmm_global_paras
   USE fmm_stats
   USE fmm_W_contractors, ONLY: fmm_select_W_con

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_translate_raw_moments,          &
             fmm_translate_boxed_moments,        &
             fmm_translate_parents_Vff,          &
             fmm_get_raw_Vff_from_boxed_Vff

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_W_pair_builder(scheme)

      USE fmm_W_buffer, ONLY: fmm_open_W_buffer
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CALL fmm_open_W_buffer(scheme)

   END SUBROUTINE fmm_init_W_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_W_pair_builder(scheme)

      USE fmm_W_buffer, ONLY: fmm_close_W_buffer
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CALL fmm_close_W_buffer(scheme)

   END SUBROUTINE fmm_free_W_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_W_pair(addr,r_ab,new_LMAX,old_LMAX,object,W_pair)

      IMPLICIT NONE
      TYPE(old_new),       INTENT(IN)  :: addr
      REAL(REALK),         INTENT(IN)  :: r_ab(3)
      INTEGER(INTK),       INTENT(IN)  :: new_LMAX, old_LMAX
      CHARACTER(3),        INTENT(IN)  :: object
      TYPE(T_pair_single), INTENT(OUT) :: W_pair

      W_pair%paras%ratio = one ! i.e. W_pair%r_ab is the actual vector

      ! indices to map back to actual moments in separate array
      W_pair%paras%LHS_id = addr%new
      W_pair%paras%RHS_id = addr%old
      ! orders of contraction with T-matrix
      W_pair%paras%LHS_LMAX = new_LMAX
      W_pair%paras%RHS_LMAX = old_LMAX

      ! (see fmm_W_worker)
      SELECT CASE (object)
      CASE ('qlm')
         W_pair%r_ab(:) = r_ab(:)
         W_pair%N_or_T  = 'N'
      CASE ('Vff')
         ! For Vff translations, we require W'(-r_ab)
         W_pair%r_ab(:) = -r_ab(:)
         ! Use 'T' (transpose) for W_matrix contraction with DTRMV
         W_pair%N_or_T  = 'T'
      CASE DEFAULT
         CALL fmm_quit('cannot resolve translation object in fmm_get_W_pair!')
      END SELECT

      W_pair%LMAX   = MAX(W_pair%paras%LHS_LMAX,W_pair%paras%RHS_LMAX)
      W_pair%lm_max = (1+W_pair%LMAX)**2

   END SUBROUTINE fmm_get_W_pair

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_translate_raw_moments(scheme,mms_in,mms_out)

      USE fmm_W_contractors, ONLY: fmm_set_W_con_ptrs
      USE fmm_W_buffer,      ONLY: fmm_add_to_W_buffer
      USE fmm_qlm_utils,     ONLY: fmm_get_T_sym_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
!FIXME: do not need to pass all this in! (just the moments)
      TYPE(raw_mm_data),  INTENT(IN)    :: mms_in
      TYPE(box_mm_data),  INTENT(INOUT) :: mms_out

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr
      INTEGER(INTK)       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      CALL fmm_select_W_con(scheme%W_con%ID)
      old_LMAX = scheme%raw_LMAX
      new_LMAX = scheme%trans_LMAX
      ! set pointers in W_contractor
      CALL fmm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

      ! generate W_pairs and pass them to the W_buffer
      CALL fmm_init_buffer_stats('W','RAW_BOX')
      CALL fmm_init_matrix_stats('W','RAW_BOX')
      CALL fmm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(mms_in%paras)
         ! indices to actual moments to translate from and to
         addr%old = mms_in%paras(i)%id
         addr%new = mms_in%paras(i)%map_up
         IF (addr%new == 0) CALL fmm_quit ('parameter mappings incomplete! 1')
         old_centre = mms_in%paras(i)%cntr
         new_centre = mms_in%paras(i)%box_cntr
         ! make W_pair and throw into W-buffer
         trans_vec = new_centre - old_centre
         CALL fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
         CALL fmm_add_to_W_buffer(W_pair)
      END DO

      ! scale new translated moments for T_matrix symmetry
      ! but first ensure W_buffer is empty
      CALL fmm_free_W_pair_builder(scheme)
      CALL fmm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

   END SUBROUTINE fmm_translate_raw_moments

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_translate_boxed_moments(scheme,mms_in,mms_out)

      USE fmm_W_contractors, ONLY: fmm_set_W_con_ptrs
      USE fmm_W_buffer,      ONLY: fmm_add_to_W_buffer
      USE fmm_qlm_utils,     ONLY: fmm_get_T_sym_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(box_mm_data),  INTENT(IN)    :: mms_in
      TYPE(box_mm_data),  INTENT(INOUT) :: mms_out

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr
      INTEGER(INTK)       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      CALL fmm_select_W_con(scheme%W_con%ID)
      old_LMAX = scheme%trans_LMAX
      new_LMAX = scheme%trans_LMAX
      ! set pointers in W_contractor
      CALL fmm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

      ! generate W_pairs and pass them to the W_buffer
      CALL fmm_init_buffer_stats('W','BOX_BOX')
      CALL fmm_init_matrix_stats('W','BOX_BOX')
      CALL fmm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(mms_in%RHS_paras)
         ! indices to actual moments to translate from and to
         addr%old = mms_in%RHS_paras(i)%id
         addr%new = mms_in%RHS_paras(i)%map_up
         IF (addr%new == 0) CALL fmm_quit ('parameter mappings incomplete! 2')
         old_centre = mms_in%RHS_paras(i)%cntr
         new_centre = mms_in%RHS_paras(i)%cntr_up
         ! make W_pair and throw into W-buffer
         trans_vec = new_centre - old_centre
         CALL fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
         CALL fmm_add_to_W_buffer(W_pair)
      END DO

      ! scale new translated moments for T_matrix symmetry
      ! but first ensure W_buffer is empty
      CALL fmm_free_W_pair_builder(scheme)
      CALL fmm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

   END SUBROUTINE fmm_translate_boxed_moments

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_translate_parents_Vff(level,scheme,Vff_p,Vff_c,c_box_paras)

      USE fmm_W_contractors, ONLY: fmm_set_W_con_ptrs
      USE fmm_W_buffer,      ONLY: fmm_add_to_W_buffer

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)    :: level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(IN)    :: Vff_p(:,:)
      REAL(REALK),        INTENT(INOUT) :: Vff_c(:,:)
      TYPE(box_mm_paras), INTENT(IN)    :: c_box_paras(:)   ! child

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr
      INTEGER(INTK)       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      IF (level <= 2) RETURN

      new_LMAX = scheme%trans_LMAX
      old_LMAX = scheme%trans_LMAX

      CALL fmm_select_W_con(scheme%W_con%ID)
      CALL fmm_set_W_con_ptrs(Vff_p,Vff_c)
      CALL fmm_init_buffer_stats('W','BOX_BOX')
      CALL fmm_init_matrix_stats('W','BOX_BOX')
      CALL fmm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(Vff_c,2)
         addr%new = c_box_paras(i)%id
         addr%old = c_box_paras(i)%map_up
         IF (addr%old == 0) CALL fmm_quit ('parameter mappings incomplete! 3')
         new_centre = c_box_paras(i)%cntr
         old_centre = c_box_paras(i)%cntr_up
         trans_vec  = new_centre - old_centre
         CALL fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
         CALL fmm_add_to_W_buffer(W_pair)
      END DO
      CALL fmm_free_W_pair_builder(scheme)

   END SUBROUTINE fmm_translate_parents_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_raw_Vff_from_boxed_Vff(raw_paras,scheme,boxed_Vff,Vff)

      USE fmm_W_contractors, ONLY: fmm_set_W_con_ptrs
      USE fmm_W_buffer,      ONLY: fmm_add_to_W_buffer

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: raw_paras(:)
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(IN)    :: boxed_Vff(:,:)
      REAL(REALK),        INTENT(INOUT) :: Vff(:,:)

      TYPE(T_pair_single) :: W_pair
      TYPE(old_new)       :: addr
      INTEGER(INTK)       :: i, new_LMAX, old_LMAX
      REAL(REALK)         :: new_centre(3), old_centre(3), trans_vec(3)

      new_LMAX = scheme%raw_LMAX
      old_LMAX = scheme%trans_LMAX

      CALL fmm_select_W_con(scheme%W_con%BR_ID)
      CALL fmm_set_W_con_ptrs(boxed_Vff,Vff)
      CALL fmm_init_buffer_stats('W','BOX_RAW')
      CALL fmm_init_matrix_stats('W','BOX_RAW')
      CALL fmm_init_W_pair_builder(scheme)
      DO i = 1, SIZE(raw_paras)
         addr%new = raw_paras(i)%id
         addr%old = raw_paras(i)%map_up
         IF (addr%old == 0) CALL fmm_quit ('parameter mappings incomplete! 4')
         new_centre = raw_paras(i)%cntr
         old_centre = raw_paras(i)%box_cntr
         trans_vec = new_centre - old_centre
         CALL fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
         CALL fmm_add_to_W_buffer(W_pair)
      END DO
      CALL fmm_free_W_pair_builder(scheme)

   END SUBROUTINE fmm_get_raw_Vff_from_boxed_Vff

!-------------------------------------------------------------------------------

END MODULE fmm_W_pair_builder
