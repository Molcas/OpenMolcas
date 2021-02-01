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
! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs.


!FIXME:  ************************************************************
!FIXME:  *                                                          *
!FIXME:  *  Array bounds when W-vector is ZERO !!!!!!!!!!           *
!FIXME:  *                                                          *
!FIXME:  ************************************************************


MODULE fmm_W_contractors

   USE fmm_global_paras
   USE fmm_stats
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_W_contractors,           &
             fmm_free_W_contractors,           &
             fmm_select_W_con,                 &
             fmm_set_W_con_ptrs

   ! Public variable to stop the resetting of W_con pointers with open T-buffer
   PUBLIC :: fmm_lock_W_con

   REAL(REALK), ALLOCATABLE, SAVE :: W_matrix(:,:)
   INTEGER(INTK),            SAVE :: WLDA
   ! Pointers to actual moments and potentials described elsewhere
   REAL(REALK), POINTER, SAVE :: old_ptr(:,:)
   REAL(REALK), POINTER, SAVE :: new_ptr(:,:)
   ! Diagnostic variables
   CHARACTER(LEN=11), SAVE :: W_con_stat
   LOGICAL,       SAVE :: fmm_lock_W_con

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_W_contractors(LMAX)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX

      IF (ALLOCATED(W_matrix)) CALL fmm_quit('W_matrix not deallocated!')
      ALLOCATE(W_matrix((1+LMAX)**2,(1+LMAX)**2))
      WLDA = (1+LMAX)**2
      W_matrix = zero

   END SUBROUTINE fmm_init_W_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_W_contractors

      IMPLICIT NONE
      DEALLOCATE(W_matrix)

   END SUBROUTINE fmm_free_W_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_W_con_ptrs(old,new)

      IMPLICIT NONE
      REAL(REALK), TARGET, INTENT(IN) :: old(:,:), new(:,:)

      IF (W_con_stat /= 'initialised') STOP 'no W_contractor preselected!'
      IF (fmm_lock_W_con) STOP 'W_buffer not empty! Cannot reset W_con!'
      NULLIFY (old_ptr, new_ptr)
      old_ptr => old(:,:)
      new_ptr => new(:,:)

   END SUBROUTINE fmm_set_W_con_ptrs

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_check_W_status

      IMPLICIT NONE
      IF ((.NOT.ASSOCIATED(old_ptr)) .OR. (.NOT.ASSOCIATED(new_ptr))) THEN
         CALL fmm_quit('W_contractor pointers not associated as reqd.')
      END IF

   END SUBROUTINE fmm_check_W_status

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one W_pair at a time.

   SUBROUTINE fmm_W_con_direct(W_pair)

      USE fmm_W_worker, ONLY: fmm_get_ltsqr_W_matrix,       &
                             fmm_contract_Wq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair

      REAL(REALK)   :: arr_tmp(W_pair%lm_max)
      REAL(REALK)   :: r_mod2
      INTEGER(INTK) :: n,m, lm_dim, hi, p,q, LMAX,JMAX
      CHARACTER(LEN=1)  :: NT

      CALL fmm_check_W_status

      NT = W_pair%N_or_T
      p = W_pair%paras%LHS_id
      q = W_pair%paras%RHS_id
      lm_dim = W_pair%lm_max

      r_mod2 = DOT_PRODUCT(W_pair%r_ab,W_pair%r_ab)
      IF (r_mod2 > ZERO_VECT_TOL*ZERO_VECT_TOL) THEN
         LMAX = W_pair%paras%LHS_LMAX
         JMAX = W_pair%paras%RHS_LMAX
         stat_W_mat_builds = stat_W_mat_builds + one
         CALL fmm_get_ltsqr_W_matrix(LMAX,JMAX,W_pair%r_ab,W_matrix)
         IF (LMAX /= JMAX) THEN
            n = (1+JMAX)**2
            m = (1+LMAX)**2
            CALL fmm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
         ELSE
            arr_tmp(:lm_dim) = old_ptr(:lm_dim,q)
            CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
            new_ptr(:lm_dim,p) = new_ptr(:lm_dim,p) + arr_tmp(:lm_dim)
         END IF
      ELSE
         hi = (1+W_pair%paras%LHS_LMAX)**2
         new_ptr(:hi,p) = new_ptr(:hi,p) + old_ptr(:hi,q)
      END IF

   END SUBROUTINE fmm_W_con_direct

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_W_con_X(W_pairs)

      USE fmm_W_worker, ONLY: fmm_get_ltsqr_W_matrix

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: W_pairs

      REAL(REALK)   :: arr_tmp(W_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER(INTK) :: LMAX, i, p,q, lm_dim, hi
      CHARACTER(LEN=1)  :: NT

      CALL fmm_check_W_status

      LMAX = W_pairs%LMAX
      NT = W_pairs%N_or_T
      r_pq_mod(:) = W_pairs%r_ab(:)
      lm_dim = W_pairs%lm_max
      lastlen = zero

      DO i = 1, SIZE(W_pairs%paras)

         p  = W_pairs%paras(i)%LHS_id
         q  = W_pairs%paras(i)%RHS_id

         hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
         IF (hi < SIZE(arr_tmp)) arr_tmp = zero
         arr_tmp(:hi) = old_ptr(:hi,q)

         IF(ABS(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) THEN
            ! zero translation, so just copy over old to new
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)
            lastlen = W_pairs%paras(i)%ratio
            CYCLE
         END IF
         ! get matrix if different
         IF(ABS(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * W_pairs%paras(i)%ratio
            lastlen = W_pairs%paras(i)%ratio
            stat_W_mat_builds = stat_W_mat_builds + one
            CALL fmm_get_ltsqr_W_matrix(LMAX,LMAX,r_pq,W_matrix)
         END IF
         ! now translate
         hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
         IF (hi < lm_dim) arr_tmp = zero
         arr_tmp(:hi) = old_ptr(:hi,q)
         CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
         hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
         new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)

      END DO

   END SUBROUTINE fmm_W_con_X

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_W_con_FAST(W_pairs)

      USE fmm_W_worker, ONLY: fmm_get_ltsqr_W_matrix, fmm_contract_Wq

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: W_pairs

      REAL(REALK)   :: arr_tmp(W_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER(INTK) :: LMAX,JMAX, n,m,lm_dim, i, p,q, hi
      CHARACTER(LEN=1)  :: NT

      CALL fmm_check_W_status

      NT = W_pairs%N_or_T
      r_pq_mod(:) = W_pairs%r_ab(:)
      lastlen = zero
      LMAX = W_pairs%LHS_LMAX
      JMAX = W_pairs%RHS_LMAX
      n = (1+JMAX)**2
      m = (1+LMAX)**2
      lm_dim = W_pairs%lm_max

      DO i = 1, SIZE(W_pairs%paras)

         p  = W_pairs%paras(i)%LHS_id
         q  = W_pairs%paras(i)%RHS_id
         IF(ABS(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) THEN
            ! zero translation, so just copy over old to new
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + old_ptr(:hi,q)
            lastlen = W_pairs%paras(i)%ratio
            CYCLE
         END IF
         IF(ABS(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * W_pairs%paras(i)%ratio
            lastlen = W_pairs%paras(i)%ratio
            stat_W_mat_builds = stat_W_mat_builds + one
            CALL fmm_get_ltsqr_W_matrix(LMAX,JMAX,r_pq,W_matrix)
         END IF
         IF (LMAX /= JMAX) THEN
            CALL fmm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
         ELSE
            hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
            IF (hi < lm_dim) arr_tmp = zero
            arr_tmp(:hi) = old_ptr(:hi,q)
            CALL DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
            hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
            new_ptr(:hi,p) = new_ptr(:hi,p) + arr_tmp(:hi)
         END IF

      END DO

   END SUBROUTINE fmm_W_con_FAST

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_W_con_BOUNDARY(W_pairs)

      USE fmm_W_worker, ONLY: fmm_get_boundary_W_matrix

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: W_pairs

      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER(INTK) :: LMAX, n, i, p,q

      CALL fmm_check_W_status

      r_pq_mod(:) = W_pairs%r_ab(:)
      lastlen = zero
      LMAX = W_pairs%RHS_LMAX
      n = (1+LMAX)**2

      DO i = 1, SIZE(W_pairs%paras)

         p  = W_pairs%paras(i)%LHS_id
         q  = W_pairs%paras(i)%RHS_id
         IF(ABS(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) THEN
            ! zero translation, so just copy over old to new
            new_ptr(1,p) = new_ptr(1,p) + old_ptr(1,q)
            lastlen = W_pairs%paras(i)%ratio
            CYCLE
         END IF
         IF(ABS(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            r_pq = r_pq_mod * W_pairs%paras(i)%ratio
            lastlen = W_pairs%paras(i)%ratio
            stat_W_mat_builds = stat_W_mat_builds + one
            CALL fmm_get_boundary_W_matrix(LMAX,r_pq,W_matrix)
         END IF

         ! Always perform transpose contraction to get boundary potential
         new_ptr(1,p) = new_ptr(1,p)+DOT_PRODUCT(W_matrix(1:n,1),old_ptr(:,q))

      END DO

   END SUBROUTINE fmm_W_con_BOUNDARY

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_select_W_con(W_con_ID)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: W_con_ID
      EXTERNAL fmm_store_w_contractor

      IF (.NOT.ALLOCATED(W_matrix)) CALL fmm_quit('W_matrix not allocated!')

      SELECT CASE (W_con_ID)
      CASE (W_CONTRACTOR_DIRECT)
         CALL fmm_store_w_contractor(fmm_W_con_DIRECT)
      CASE (W_CONTRACTOR_X)
         CALL fmm_store_w_contractor(fmm_W_con_X)
      CASE (W_CONTRACTOR_FAST)
         CALL fmm_store_w_contractor(fmm_W_con_FAST)
      CASE (W_CONTRACTOR_BOUNDARY)
         CALL fmm_store_w_contractor(fmm_W_con_BOUNDARY)
      CASE DEFAULT
         CALL fmm_quit ('invalid W_contractor requested!')
      END SELECT
      ! initialise diagnostics
      W_con_stat = 'initialised'
      W_con_stat = 'initialised'
      fmm_lock_W_con = .FALSE.

   END SUBROUTINE fmm_select_W_con

!-------------------------------------------------------------------------------

END MODULE fmm_W_contractors

