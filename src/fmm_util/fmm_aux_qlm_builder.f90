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
MODULE fmm_aux_qlm_builder

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_aux_qlm

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_aux_qlm(scheme,LHS_mms,RHS_mms)

      USE fmm_qlm_utils, ONLY: fmm_renormalise_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS_mms, RHS_mms

      ! We use the SCALED solid harmonic formulation throughout
      CALL fmm_renormalise_qlm(scheme%raw_LMAX,LHS_mms%qlm)
      CALL fmm_renormalise_qlm(scheme%raw_LMAX,RHS_mms%qlm)

      ! Identify batches of unique centre/extent and sort them in
      ! preparation for packing
      CALL sort_centres_and_assign_batches(scheme,LHS_mms,RHS_mms)

      ! Get parameter mappings and preconditioned W and T moments
      CALL get_RHS_data(scheme,RHS_mms)
      CALL get_LHS_data(scheme,LHS_mms)

      ! All relevant data should now be held in %qlm_T and %qlm_W
      DEALLOCATE (LHS_mms%qlm, RHS_mms%qlm)
      NULLIFY (LHS_mms%qlm, RHS_mms%qlm)

   END SUBROUTINE fmm_get_aux_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE sort_centres_and_assign_batches(scheme,LHS,RHS)

      USE fmm_qlm_utils, ONLY: fmm_sort_paras_wrt_centre, fmm_assign_batches

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS, RHS

      IF (scheme%pack_LHS) THEN
         CALL fmm_sort_paras_wrt_centre(1_INTK,LHS%paras)
         CALL fmm_assign_batches(LHS%paras)
      END IF
      IF (scheme%pack_RHS) THEN
         CALL fmm_sort_paras_wrt_centre(1_INTK,RHS%paras)
         CALL fmm_assign_batches(RHS%paras)
      END IF

   END SUBROUTINE sort_centres_and_assign_batches

!-------------------------------------------------------------------------------

   SUBROUTINE get_LHS_data(scheme,LHS)

      USE fmm_qlm_utils, ONLY: fmm_factor_in_dens, fmm_pack_raw_parameters

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: LHS
      INTEGER(INTK) :: i, ndim, qlm_dim, alloc_error

      ! Get LHS parameters from raw data
      IF (scheme%pack_LHS) THEN
         CALL fmm_pack_raw_parameters(LHS)
      END IF

      ! Sync (packed)parameters:(packed)moments mappings
      DO i = 1, SIZE(LHS%paras)
         LHS%paras(i)%id = i
      END DO

      ! LHS preconditioning for T-matrix
      SELECT CASE (scheme%T_con%LHS_mm_type)
         CASE (USE_RAW_QLM)
            ! Note that dimensions of paras and qlm,2 may be different
            qlm_dim = SIZE(LHS%qlm,2)
            ndim = SIZE(LHS%qlm,1)
            WRITE(LUPRI,*) 'LHS%qlm_T: Attempting to allocate',  &
                           MAX(1,qlm_dim*ndim*8/1000000), 'MB of memory...'
            ALLOCATE (LHS%qlm_T(ndim,qlm_dim), STAT=alloc_error)
            IF (alloc_error /= 0) WRITE(LUPRI,*) '... Failed!'
            LHS%qlm_T(:,:) = LHS%qlm(:,:)
         CASE DEFAULT
            CALL fmm_quit('cannot reconcile LHS_mm_type')
      END SELECT

      ! Factorise in density if required
      IF ( scheme%LHS_dens ) THEN
         CALL fmm_factor_in_dens(LHS%dens,LHS%qlm_T)
         DEALLOCATE (LHS%dens)
         NULLIFY (LHS%dens)
      END IF

   END SUBROUTINE get_LHS_data

!-------------------------------------------------------------------------------

   SUBROUTINE get_RHS_data(scheme,RHS)

      USE fmm_qlm_utils, ONLY: fmm_get_T_sym_qlm,     &
                               fmm_pack_raw_moments,  &
                               fmm_factor_in_dens

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: RHS
      INTEGER(INTK) :: LMAX, i, qlm_dim, ndim, alloc_error
      LOGICAL :: dens
      REAL(REALK) :: thr

      LMAX = scheme%raw_LMAX

      IF ( scheme%pack_RHS ) THEN
         ! We now pack all moments belonging to
         ! the same batch (same centre,extent);
         ! Also factor in density and perform
         ! density-based screening if flagged
         dens = scheme%RHS_dens
         thr = scheme%dens_screen_thr
         CALL fmm_pack_raw_moments(RHS,dens,thr)
      END IF

      ndim = (1+LMAX)**2
      qlm_dim = SIZE(RHS%qlm,2)
      WRITE(LUPRI,*) 'RHS%qlm_W: Attempting to allocate',  &
                     MAX(1,qlm_dim*ndim*8/1000000), 'MB of memory...'
      ALLOCATE (RHS%qlm_W(ndim,qlm_dim), STAT=alloc_error)
      IF (alloc_error /= 0) WRITE(LUPRI,*) '... Failed!'
      RHS%qlm_W(:,:) = RHS%qlm(:,:)

      IF ( .NOT. scheme%pack_RHS ) THEN
         IF ( scheme%RHS_dens ) CALL fmm_factor_in_dens(RHS%dens,RHS%qlm_W)
      END IF

      ! Density nolonger required
      IF ( scheme%RHS_dens ) THEN
         DEALLOCATE(RHS%dens)
         NULLIFY(RHS%dens)
      END IF

      IF (qlm_dim /= SIZE(RHS%paras)) CALL fmm_quit('error in RHS data')
      ! Resync parameters:moments mappings
      DO i = 1, qlm_dim
         RHS%paras(i)%id = i
      END DO

      ! RHS preconditioning for T-matrix
      SELECT CASE (scheme%T_con%RHS_mm_type)
         CASE (USE_RAW_QLM)
            RHS%qlm_T => RHS%qlm_W(:,:)
         CASE (USE_T_SYM_QLM)
            ALLOCATE (RHS%qlm_T(ndim,qlm_dim))
             ! build %qlm_T by rescaling significant %qlm_W
            CALL fmm_get_T_sym_qlm(LMAX,RHS%qlm_W,RHS%qlm_T)
         CASE DEFAULT
            CALL fmm_quit('cannot reconcile RHS_mm_type')
      END SELECT

   END SUBROUTINE get_RHS_data

!-------------------------------------------------------------------------------

END MODULE fmm_aux_qlm_builder
