!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2002, Mark A. Watson                                   *
!***********************************************************************
! fmm_driver:
! General module for driving multipole method calculations.

MODULE fmm_driver

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_build_J_matrix,                 &
             fmm_get_multipole_potential

   ! The multipole potential
   REAL(REALK), POINTER, SAVE :: Vff(:,:)
   ! The (modified) multipole expansions
   TYPE(raw_mm_data),    SAVE :: LHS_mms, RHS_mms

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_driver(scheme,dens)

      USE fmm_qlm_builder,     ONLY: fmm_get_raw_qlm
      USE fmm_aux_qlm_builder, ONLY: fmm_get_aux_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(IN)    :: dens(:,:)

      NULLIFY (Vff)

      ! Get basic multipole data
      CALL fmm_get_raw_qlm(scheme,dens,LHS_mms,RHS_mms)
      ! Now ensure appropriate prefactors, normalisation and packing
      CALL fmm_get_aux_qlm(scheme,LHS_mms,RHS_mms)
      ! Allocate the far field potential
      CALL fmm_allocate_Vff(scheme)

   END SUBROUTINE fmm_init_driver

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_driver

      USE fmm_qlm_builder, ONLY: fmm_deallocate_qlm

      IMPLICIT NONE
      DEALLOCATE(Vff)
      NULLIFY (Vff)
      CALL fmm_deallocate_qlm(LHS_mms,RHS_mms)

   END SUBROUTINE fmm_free_driver

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_allocate_Vff(scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK) :: lm_dim, mms_dim, alloc_error

      IF (.NOT.ASSOCIATED(LHS_mms%paras)) STOP 'mms ptrs not set in fmm_driver!'
      IF (ASSOCIATED(Vff)) CALL fmm_quit('Vff should NOT be allocated already!')

      mms_dim = SIZE(LHS_mms%paras)
      ! Note we shouldn't use this array for translated, BOXED potentials
      lm_dim = (1+ scheme%raw_LMAX)**2
      IF (scheme%job_type == GFC_FMM) lm_dim = 1

      WRITE(LUPRI,*) 'Vff: Attempting to allocate',  &
                     MAX(1,lm_dim*mms_dim*8/1000000), 'MB of memory...'
      ALLOCATE (Vff(lm_dim,mms_dim), STAT=alloc_error)
      IF (alloc_error /= 0) WRITE(LUPRI,*) '... Failed!'

      ! Must zero out since Vff is built additively in general
      Vff(:,:) = zero

   END SUBROUTINE fmm_allocate_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_J_via_raw_potentials(scheme,dens,J_matrix,energy,txt)

      USE fmm_Vff_driver, ONLY: fmm_get_Vff
      USE fmm_J_builder,  ONLY: fmm_get_J_from_Vff,       &
                                fmm_get_J_from_pkd_Vff,   &
                                fmm_get_E_from_pkd_Vff,   &
                                fmm_get_E_from_Vff
      USE fmm_qlm_utils,  ONLY: fmm_factor_in_dens

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(IN)    :: dens(:,:)
      REAL(REALK),        INTENT(OUT)   :: J_matrix(:,:)
      REAL(REALK),        INTENT(OUT)   :: energy
      CHARACTER(*),       INTENT(OUT)   :: txt

      ! We only have the density on RHS when getting J-matrix
      ! via far-field potential
      scheme%LHS_dens = .FALSE.
      scheme%RHS_dens = .TRUE.

      ! Prepare moments and allocate potential
      CALL fmm_init_driver(scheme,dens)
      ! Get potential
      CALL fmm_get_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

      ! Get J-matrix
      J_matrix = zero
      energy = zero
      IF (scheme%pack_LHS) THEN
         CALL fmm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! Get energy after factoring in density to LHS
         CALL fmm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL fmm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,txt)
      ELSE
         CALL fmm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! Get energy after factoring in density to LHS
         CALL fmm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL fmm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,txt)
      END IF

      CALL fmm_free_driver

   END SUBROUTINE fmm_get_J_via_raw_potentials

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_J_matrix(n_el,dens,J_matrix)

      USE fmm_stats
      USE fmm_scheme_builder, ONLY: fmm_get_scheme

      IMPLICIT NONE
      CHARACTER(6), INTENT(IN)  :: n_el
      REAL(REALK),  INTENT(IN)  :: dens(:,:)
      REAL(REALK),  INTENT(OUT) :: J_matrix(:,:)

      TYPE(scheme_paras), POINTER :: scheme
      CHARACTER(36) :: E_text
      REAL(REALK)   :: energy, T0, fmm_second, TTOT

      T0 = fmm_second()

      CALL fmm_get_scheme(scheme)

      SELECT CASE (n_el)
      CASE('ONE_EL')
         CALL fmm_quit('nuclear moments not available!')
         scheme%LHS_mm_range = ELECTRONIC_ONLY
         scheme%RHS_mm_range = NUCLEAR_ONLY
      CASE('TWO_EL')
         scheme%LHS_mm_range = ELECTRONIC_ONLY
         scheme%RHS_mm_range = ELECTRONIC_ONLY
      CASE('FULL_J')
         CALL fmm_quit('nuclear moments not available!')
         scheme%LHS_mm_range = ELECTRONIC_ONLY
         scheme%RHS_mm_range = ALL_MOMENTS
      CASE DEFAULT
         CALL fmm_quit ('require 1, 2, or full J_matrix build!')
      END SELECT

      CALL fmm_get_J_via_raw_potentials(scheme,dens,J_matrix,energy,E_text)
      WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy

      TTOT = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_J_matrix', TTOT, LUPRI)
      CALL fmm_print_stats

   END SUBROUTINE fmm_build_J_matrix

!-------------------------------------------------------------------------------
! Note that the potential returned by this routine can only be
! contracted with appropriately scaled multipole moments

   SUBROUTINE fmm_get_multipole_potential(mode,dens,potential)

      USE fmm_stats
      USE fmm_scheme_builder, ONLY: fmm_get_scheme
      USE fmm_boundary,       ONLY: fmm_opt_near_field
      USE fmm_Vff_driver,     ONLY: fmm_get_Vff

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: mode
      REAL(REALK),   INTENT(IN)  :: dens(:,:)
      REAL(REALK),   INTENT(OUT) :: potential(:,:)

      TYPE(scheme_paras), POINTER :: scheme
      REAL(REALK) :: T0, fmm_second, TTOT
      INTEGER(INTK) :: lmdim

      T0 = fmm_second()

      CALL fmm_get_scheme(scheme)

      scheme%LHS_mm_range = NUCLEAR_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
      scheme%LHS_dens = .FALSE.
      scheme%RHS_dens = .TRUE.
      scheme%pack_LHS = .FALSE.

      ! Prepare moments and allocate potential
      CALL fmm_init_driver(scheme,dens)

      ! Test if we can skip the near-field interactions
      IF (mode == GFC_FMM) THEN
         CALL fmm_opt_near_field(scheme,LHS_mms%paras,RHS_mms%paras)
      END IF

      ! Get potential
      CALL fmm_get_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

      ! Note we assume here that the LHS hasn't been rearranged
      lmdim = SIZE(potential,1)
      IF (SIZE(potential,2) /= SIZE(Vff,2)) CALL fmm_quit('bounds: potential')
      potential(:,:) = Vff(1:lmdim,:)

      CALL fmm_free_driver

      TTOT = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_multipole_potential', TTOT, LUPRI)
      CALL fmm_print_stats

   END SUBROUTINE fmm_get_multipole_potential

!-------------------------------------------------------------------------------

END MODULE fmm_driver

