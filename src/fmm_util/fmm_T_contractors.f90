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
MODULE fmm_T_contractors

! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs to generate
! an energy or J-matrix directly (e.g.to exploit symmetry).

   USE fmm_global_paras
   USE fmm_stats
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_T_contractors,          &
             fmm_free_T_contractors,          &
             fmm_select_T_con,                &
             fmm_set_T_con_ptrs

   ! Public variable to stop the resetting of T_con pointers with open T-buffer
   PUBLIC :: fmm_lock_T_con

   REAL(REALK), ALLOCATABLE, SAVE :: T_matrix(:,:)
   INTEGER(INTK),            SAVE :: TLDA
   ! when building multiple T matrices together
   REAL(REALK), ALLOCATABLE, SAVE :: T_mats(:,:,:)

   ! Pointers to actual moments and potentials described elsewhere
   REAL(REALK), POINTER, SAVE :: Vff_ptr(:,:)
   REAL(REALK), POINTER, SAVE :: qlm_ptr(:,:)

   ! Diagnostic variables
   CHARACTER(LEN=11), SAVE :: T_con_stat
   LOGICAL,       SAVE :: fmm_lock_T_con

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_T_contractors(scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK) :: LMAX, T_con
      INTEGER(INTK) :: lm_max

      LMAX = scheme%trans_lmax
      lm_max = (LMAX+1)**2
      IF (scheme%phase == NEAR_FIELD) THEN
         T_con = scheme%T_con%NF_id
      ELSE
         T_con = scheme%T_con%FF_id
      END IF

      SELECT CASE (T_con)
         CASE (T_CONTRACTOR_MULTI)
            IF (ALLOCATED(T_mats)) CALL fmm_quit('T_mats not deallocated!')
!FIXME: will need to change TMATM_DF if changed elsewhere
            ALLOCATE(T_mats(TMATM_DF,lm_max,lm_max))
            T_mats = zero
         CASE (T_CONTRACTOR_BOUNDARY)
            IF (ALLOCATED(T_matrix)) CALL fmm_quit('T_matrix not deallocated!')
            ALLOCATE(T_matrix(lm_max,1))
            T_matrix = zero
         CASE DEFAULT
            IF (ALLOCATED(T_matrix)) CALL fmm_quit('T_matrix not deallocated!')
            ALLOCATE(T_matrix(lm_max,lm_max))
            T_matrix = zero
      END SELECT
      TLDA = lm_max

      ! code for statistics only
      CALL fmm_init_matrix_stats('T')

   END SUBROUTINE fmm_init_T_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_T_contractors

      IMPLICIT NONE

      IF (ALLOCATED(T_matrix)) DEALLOCATE(T_matrix)
      IF (ALLOCATED(T_mats)) DEALLOCATE(T_mats)

   END SUBROUTINE fmm_free_T_contractors

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_T_con_ptrs(Vff,qlm)

      IMPLICIT NONE
      REAL(REALK), TARGET, INTENT(IN) :: Vff(:,:), qlm(:,:)

      IF (T_con_stat /= 'initialised') STOP 'no T_contractor preselected!'
      IF (fmm_lock_T_con) STOP 'T_buffer not empty! Cannot reset T_con!'
      NULLIFY (Vff_ptr, qlm_ptr)
      Vff_ptr => Vff(:,:)
      qlm_ptr => qlm(:,:)

   END SUBROUTINE fmm_set_T_con_ptrs

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one T_pair at a time.
! Only generates a square T-matrix with dim = MAX(LHS_LMAX,RHS_LMAX)
! Generates local potential to order MAX(LHS_LMAX,RHS_LMAX)
! Adds to Vff only up to LHS_LMAX.
! FIXME: because T_matrix is only built for (l+j)<=LMAX then with dynamic
! contraction we need LEXTRA>0 even for primitive basis sets and FQUAD.

   SUBROUTINE fmm_T_con_DIRECT(T_pair)

      USE fmm_T_worker, ONLY: fmm_get_SPLTSQ_T_matrix,  &
                              fmm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair

      REAL(REALK)   :: Vff_tmp(T_pair%lm_max)
      INTEGER(INTK) :: weight, iLHS,iRHS, hi

      stat_T_mat_builds = stat_T_mat_builds + one
      CALL fmm_get_SPLTSQ_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

      iRHS = T_pair%paras%RHS_id
      hi = T_pair%lm_max

      Vff_tmp = fmm_contract_Tq(T_pair%LMAX,qlm_ptr(:hi,iRHS),T_matrix)
      iLHS = T_pair%paras%LHS_id
      hi = (1+T_pair%paras%LHS_LMAX)**2
      weight = T_pair%paras%weight
      !WRITE(LUPRI,*) iRHS, qlm_ptr(1,iRHS), weight, Vff_tmp(:hi)
      Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + weight*Vff_tmp(:hi)

   END SUBROUTINE fmm_T_con_DIRECT

!-------------------------------------------------------------------------------
! Modification of T_con_DIRECT for evaluation of boundary potential
! when only the first element of Vff needs to be built

   SUBROUTINE fmm_T_con_BOUNDARY(T_pair)

      USE fmm_T_worker, ONLY: fmm_get_boundary_T_matrix

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair

      REAL(REALK)   :: Vff_tmp
      INTEGER(INTK) :: weight, iLHS,iRHS, hi

      stat_T_mat_builds = stat_T_mat_builds + one
      CALL fmm_get_boundary_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

      iRHS = T_pair%paras%RHS_id
      hi = T_pair%lm_max

      Vff_tmp = half*DOT_PRODUCT(qlm_ptr(:hi,iRHS),T_matrix(:hi,1))

      iLHS = T_pair%paras%LHS_id
      hi = (1+T_pair%paras%LHS_LMAX)**2
      weight = T_pair%paras%weight
      !WRITE(LUPRI,*) iRHS, qlm_ptr(1,iRHS), weight, Vff_tmp

      Vff_ptr(1,iLHS) = Vff_ptr(1,iLHS) + weight*Vff_tmp

   END SUBROUTINE fmm_T_con_BOUNDARY

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_T_con_TREE(T_pairs)

      USE fmm_T_worker, ONLY: fmm_get_SPLTSQ_T_matrix,    &
                              fmm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: T_pairs

      REAL(REALK)   :: Vff_tmp(T_pairs%lm_max)
      REAL(REALK)   :: lastlen, r_pq(3), r_pq_mod(3)
      INTEGER(INTK) :: LMAX, weight, i, p,q, hi

      r_pq_mod(:) = T_pairs%r_ab(:)
      LMAX = T_pairs%LMAX
      lastlen = zero

      DO i = 1, SIZE(T_pairs%paras)

         p = T_pairs%paras(i)%LHS_id
         q = T_pairs%paras(i)%RHS_id

         IF(ABS(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            stat_T_mat_builds = stat_T_mat_builds + one
            r_pq = r_pq_mod * T_pairs%paras(i)%ratio
            CALL fmm_get_SPLTSQ_T_matrix(LMAX,r_pq,T_matrix)
            lastlen = T_pairs%paras(i)%ratio
         END IF

         hi = T_pairs%lm_max
         Vff_tmp = fmm_contract_Tq(LMAX,qlm_ptr(:hi,q),T_matrix)
         hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
         weight = T_pairs%paras(i)%weight
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + weight*Vff_tmp(:hi)

      END DO

   END SUBROUTINE fmm_T_con_TREE

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_T_con_SCALE(T_pairs)

      USE fmm_T_worker, ONLY: fmm_get_SPLTSQ_T_matrix,    &
                              fmm_contract_Tq,            &
                              fmm_scale_vec

      IMPLICIT NONE
      TYPE(T_pair_batch), INTENT(IN) :: T_pairs

      INTEGER(INTK) :: LMAX, i, p,q, hi, lastq
      REAL(REALK)   :: ratio, lastlen, pref
      REAL(REALK)   :: Vff_tmp(T_pairs%items(1)%lm_max)
      REAL(REALK)   :: scaled_qlm(T_pairs%items(1)%lm_max)
      REAL(REALK)   :: scale_vec(T_pairs%items(1)%lm_max)
      LOGICAL       :: new_vec

      pref = one
      lastq = -1
      lastlen = one
      scale_vec(:) = one
      LMAX = T_pairs%items(1)%LMAX
      stat_T_mat_builds = stat_T_mat_builds + one
      CALL fmm_get_SPLTSQ_T_matrix(LMAX,T_pairs%items(1)%r_ab,T_matrix)

      iloop: DO i = 1, T_pairs%ndim

         p = T_pairs%items(i)%paras%LHS_id
         q = T_pairs%items(i)%paras%RHS_id

         new_vec = .FALSE.
         ratio = T_pairs%items(i)%paras%ratio
         IF (ABS(ratio-lastlen) > DISTINCT_T_TOL) THEN
            CALL fmm_scale_vec(LMAX,ratio,scale_vec,pref)
            lastlen = ratio
            new_vec = .TRUE.
         END IF

         hi = T_pairs%items(i)%lm_max
         IF ( new_vec .OR. (q /= lastq) ) THEN
            scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
            lastq = q
         END IF

         Vff_tmp = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

         hi = (1+T_pairs%items(i)%paras%LHS_LMAX)**2
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + pref*scale_vec(:hi)*Vff_tmp(:hi)

      END DO iloop

   END SUBROUTINE fmm_T_con_SCALE

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_T_con_SCALE_TREE(T_pairs)

      USE fmm_T_worker, ONLY: fmm_get_SPLTSQ_T_matrix,    &
                              fmm_contract_Tq,            &
                              fmm_scale_vec

      IMPLICIT NONE
      TYPE(T_pair_list), INTENT(IN) :: T_pairs

      INTEGER(INTK) :: LMAX, i, p,q, hi, lastq
      REAL(REALK)   :: weight, lastlen, pref
      REAL(REALK)   :: Vff_tmp(T_pairs%lm_max)
      REAL(REALK)   :: scaled_qlm(T_pairs%lm_max)
      REAL(REALK)   :: scale_vec(T_pairs%lm_max)
      LOGICAL       :: new_vec

      pref = one
      lastq = -1
      lastlen = one
      scale_vec(:) = one
      LMAX = T_pairs%LMAX
      stat_T_mat_builds = stat_T_mat_builds + one
      CALL fmm_get_SPLTSQ_T_matrix(LMAX,T_pairs%r_ab,T_matrix)

      iloop: DO i = 1, SIZE(T_pairs%paras)

         p = T_pairs%paras(i)%LHS_id
         q = T_pairs%paras(i)%RHS_id

         new_vec = .FALSE.
         IF(ABS(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) THEN
            CALL fmm_scale_vec(LMAX,T_pairs%paras(i)%ratio,scale_vec,pref)
            lastlen = T_pairs%paras(i)%ratio
            new_vec = .TRUE.
         END IF

         hi = T_pairs%lm_max
         IF ( new_vec .OR. (q /= lastq) ) THEN
            scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
            lastq = q
         END IF

         Vff_tmp = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

         hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
         weight = pref*T_pairs%paras(i)%weight
         Vff_ptr(:hi,p) = Vff_ptr(:hi,p) + weight*scale_vec(:hi)*Vff_tmp(:hi)

      END DO iloop

   END SUBROUTINE fmm_T_con_SCALE_TREE

!-------------------------------------------------------------------------------
! Special contractor designed to take batch of T-pairs ordered in pairs
! such that ( b T1 1; a -T1 2; c T2 1; a -T2 3 ....)
! to halve the number of T matrix builds (using fmm_scale...)
! But only the common RHS qlm can be contracted simultaneously which slows
! it down a lot.
!
!   SUBROUTINE fmm_T_con_MULTI(T_pairs)
!
!      USE fmm_T_worker_multi, ONLY: fmm_get_SPLTSQ_T_matrices,    &
!                                   fmm_contract_multi_Tq
!      USE fmm_T_worker,       ONLY: fmm_contract_Tq, fmm_scale_vec
!
!      IMPLICIT NONE
!      TYPE(T_pair_batch), INTENT(IN) :: T_pairs
!
!      REAL(REALK), ALLOCATABLE :: Vff_tmp(:,:)
!
!      REAL(REALK), ALLOCATABLE :: scaled_qlm(:)
!      REAL(REALK), ALLOCATABLE :: scale_vec(:)
!
!      REAL(REALK)   :: T_vectors((T_pairs%ndim/2),3)
!      REAL(REALK)   :: pref, weight
!      INTEGER(INTK) :: i,j, iLHS,iRHS,iRHS_last, hi, LMAX, nT
!
!      ! FIRST DO MULTIPLE T contracted with SAME RHS
!      !---------------------------------------------
!
!      IF (BTEST(T_pairs%ndim,0)) CALL fmm_quit('ndim not EVEN!')
!      nT = T_pairs%ndim/2
!
!      LMAX = 0
!      iRHS_last = T_pairs%items(1)%paras%RHS_id
!      DO i = 1, nT
!         j = 2*i-1
!         LMAX = MAX(LMAX,T_pairs%items(i)%LMAX)
!         iRHS = T_pairs%items(j)%paras%RHS_id
!         ! get *distinct* T-vectors as every second item in batch
!         T_vectors(i,:) = T_pairs%items(j)%r_ab(:)
!         IF (iRHS /= iRHS_last) CALL fmm_quit('must have same qlm on RHS')
!         iRHS_last = iRHS
!      END DO
!
!      stat_T_mat_builds = stat_T_mat_builds + nT
!      CALL fmm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats(:nT,:,:))
!
!      hi = (1+LMAX)**2
!      ALLOCATE( Vff_tmp(nT,hi) )
!      ALLOCATE( scaled_qlm(hi) )
!      ALLOCATE( scale_vec(hi) )
!      Vff_tmp(:,:) = zero
!      scale_vec(:) = one
!
!      iRHS = T_pairs%items(1)%paras%RHS_id
!
!      Vff_tmp(:,:hi) = fmm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),    &
!                                            T_mats(:nT,:,:),nT)
!
!      DO i = 1, nT
!         j = 2*i-1
!         iLHS = T_pairs%items(j)%paras%LHS_id
!         hi = (1+T_pairs%items(j)%paras%LHS_LMAX)**2
!         weight = T_pairs%items(j)%paras%weight
!         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + weight*Vff_tmp(i,:hi)
!      END DO
!
!      ! NOW DO remaining half of batched T-pairs corresponding to same LHS
!      !-------------------------------------------------------------------
!
!      iLHS = T_pairs%items(2)%paras%LHS_id
!      DO i = 1, nT
!         ! get T_matrix correpxonding to the minus T-vector
!         CALL fmm_scale_vec(LMAX,-one,scale_vec,pref)
!         iRHS = T_pairs%items(2*i)%paras%RHS_id
!         scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,iRHS)
!         Vff_tmp(1,:hi) = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_mats(i,:,:))
!         weight = pref
!         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)  &
!                             + weight*scale_vec(:hi)*Vff_tmp(1,:hi)
!      END DO
!
!      DEALLOCATE(Vff_tmp)
!      DEALLOCATE(scaled_qlm)
!      DEALLOCATE(scale_vec)
!
!   END SUBROUTINE fmm_T_con_MULTI

!-------------------------------------------------------------------------------
! Contractor designed to build multiple T matrices simultaneously.
! Exact performance will depend on architecture and choice of NDIM.
! Assumes RHS moments are the same, but T-vectors are different.

   SUBROUTINE fmm_T_con_MULTI(T_pairs)

      USE fmm_multiple_T_worker, ONLY: fmm_get_SPLTSQ_T_matrices,    &
                                       fmm_contract_multi_Tq

      IMPLICIT NONE
      TYPE(T_pair_batch), INTENT(IN) :: T_pairs

      REAL(REALK)   :: Vff_tmp(T_pairs%ndim,T_pairs%items(1)%lm_max)
      REAL(REALK)   :: T_vectors(T_pairs%ndim,3)
      INTEGER(INTK) :: i, iLHS,iRHS, hi, LMAX, nT

      nT = T_pairs%ndim
      LMAX = T_pairs%items(1)%LMAX

      DO i = 1, nT
         T_vectors(i,:) = T_pairs%items(i)%r_ab
         iRHS = T_pairs%items(i)%paras%RHS_id
         IF ( iRHS /=  T_pairs%items(MAX(i-1,1))%paras%RHS_id ) THEN
            CALL fmm_quit('RHS moments not sorted in fmm_T_con_MULTI')
         END IF
      END DO

      stat_T_mat_builds = stat_T_mat_builds + nT
      CALL fmm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats(:nt,:,:))

      iRHS = T_pairs%items(1)%paras%RHS_id
      hi = (1+LMAX)**2

      Vff_tmp(:,:hi) = fmm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),    &
                                            T_mats(:nT,:,:),nT)

      DO i = 1, nT
         iLHS = T_pairs%items(i)%paras%LHS_id
         Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + Vff_tmp(i,:hi)
      END DO

   END SUBROUTINE fmm_T_con_MULTI

!-------------------------------------------------------------------------------
! Builds full T-matrix for exact contraction at low orders

   SUBROUTINE fmm_T_con_FULL(T_pair)

      USE fmm_T_worker, ONLY: fmm_get_FLTSQ_T_matrix, fmm_postfac_Vff,  &
                              fmm_contract_Tq

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair

      REAL(REALK)   :: Vff_tmp(T_pair%lm_max)
      INTEGER(INTK) :: iLHS,iRHS, hi

      stat_T_mat_builds = stat_T_mat_builds + one
      CALL fmm_get_FLTSQ_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

      iRHS = T_pair%paras%RHS_id
      hi = T_pair%lm_max
      CALL DSYMV('L',hi,1d0,T_matrix,TLDA,qlm_ptr(:,iRHS),1,0D0,Vff_tmp,1)

      iLHS = T_pair%paras%LHS_id
      hi = (1+T_pair%paras%LHS_LMAX)**2
      CALL fmm_postfac_Vff(T_pair%paras%LHS_LMAX,Vff_tmp)
      Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS) + Vff_tmp(:hi)

   END SUBROUTINE fmm_T_con_FULL

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_select_T_con(scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK) :: T_con_ID
      EXTERNAL fmm_store_t_contractor

      IF (scheme%phase == NEAR_FIELD) THEN
         T_con_ID = scheme%T_con%NF_id
      ELSE
         T_con_ID = scheme%T_con%FF_id
      END IF

      SELECT CASE (T_con_ID)
         CASE (T_CONTRACTOR_DIRECT)
            CALL fmm_store_t_contractor(fmm_T_con_DIRECT)
         CASE (T_CONTRACTOR_BOUNDARY)
            CALL fmm_store_t_contractor(fmm_T_con_BOUNDARY)
         CASE (T_CONTRACTOR_TREE)
            CALL fmm_store_t_contractor(fmm_T_con_TREE)
         CASE (T_CONTRACTOR_SCALE_TREE)
            CALL fmm_store_t_contractor(fmm_T_con_SCALE_TREE)
         CASE (T_CONTRACTOR_SCALE)
            CALL fmm_store_t_contractor(fmm_T_con_SCALE)
         CASE (T_CONTRACTOR_MULTI)
            CALL fmm_store_t_contractor(fmm_T_con_MULTI)
         CASE (T_CONTRACTOR_FULL)
            CALL fmm_store_t_contractor(fmm_T_con_FULL)
         CASE DEFAULT
            CALL fmm_quit ('invalid T_contractor requested!')
      END SELECT
      ! initialise diagnostics
      T_con_stat = 'initialised'
      fmm_lock_T_con = .FALSE.

   END SUBROUTINE fmm_select_T_con

!-------------------------------------------------------------------------------

END MODULE fmm_T_contractors

