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
MODULE fmm_T_pair_tests

   USE fmm_global_paras
   USE fmm_box_utils, ONLY: fmm_NF_boxes, fmm_RFF_boxes

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_T_pair_tests,         &
             fmm_close_T_pair_tests,        &
             fmm_test_and_buffer_T_pair

    ! Flag to test initialisation
   CHARACTER(11), SAVE :: init_tests

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_T_pair_tests(scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme

      IF (scheme%phase == NEAR_FIELD) THEN
         CALL fmm_store_test(fmm_test_NF_ext)
      ELSE
         SELECT CASE (scheme%algorithm)
            CASE (DO_FQ)
               CALL fmm_store_test(fmm_test_raw_FF)
            CASE (DO_BQ)
               CALL fmm_store_test(fmm_test_FF)
            CASE (DO_NlogN)
               CALL fmm_store_test(fmm_test_LFF)
            CASE (DO_FMM)
               CALL fmm_store_test(fmm_test_LFF)
            CASE DEFAULT
               CALL fmm_quit ('unable to initialise T_pair_tests')
         END SELECT
      END IF

      init_tests = 'initialised'

   END SUBROUTINE fmm_init_T_pair_tests

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_close_T_pair_tests

      IMPLICIT NONE
      IF (init_tests /= 'initialised') STOP 'must initialise pair_tests!'
      init_tests = ' '

   END SUBROUTINE fmm_close_T_pair_tests

!-------------------------------------------------------------------------------
! routine to perform T-pair test on two raw or boxed moments, and if they
! pass the pair is sent for evaluation.

   SUBROUTINE fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)

      USE fmm_T_buffer, ONLY: fmm_add_to_T_buffer

      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      INTEGER(INTK),      INTENT(IN) :: weight
      TYPE(T_pair_single) :: T_pair
      LOGICAL, EXTERNAL   :: fmm_included_pair
      EXTERNAL fmm_stored_t_pair_mould

      IF (fmm_included_pair(LHS,RHS,id)) THEN
        ! pour all relevant info into the T-pair mould
         CALL fmm_stored_t_pair_mould(LHS,RHS,id,weight,T_pair)
        ! pass single T_pair entity to contractors via buffer
         CALL fmm_add_to_T_buffer(T_pair)
      END IF

   END SUBROUTINE fmm_test_and_buffer_T_pair

!-------------------------------------------------------------------------------
! fmm_test_pass is always TRUE, bypassing normal classical/nonclassical tests

   FUNCTION fmm_test_pass(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: fmm_test_pass
      fmm_test_pass = .TRUE.

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer_array(LHS%raw_paras%id)
         CALL Unused_integer_array(RHS%raw_paras%id)
         CALL Unused_integer(id%LHS)
      END IF
   END FUNCTION fmm_test_pass

!-------------------------------------------------------------------------------
! Logical function to test if two gaussian overlaps are separated by
! more than the sum of their extents.
! fmm_test_ext is TRUE if sufficiently well-separated.

   FUNCTION fmm_test_ext(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: fmm_test_ext
      REAL(REALK) :: r_ij(3), ext_ij, r_mod, r_zero

      ext_ij = RHS%raw_paras(id%RHS)%ext + LHS%raw_paras(id%LHS)%ext
      r_ij   = RHS%raw_paras(id%RHS)%cntr - LHS%raw_paras(id%LHS)%cntr
      r_mod  = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)
      ext_ij = ext_ij*ext_ij
      r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL
      fmm_test_ext = ((r_mod - ext_ij) > r_zero)

   END FUNCTION fmm_test_ext

!-------------------------------------------------------------------------------
! Interaction pair test routine for near-field phase.
! An interaction should be treated by multipoles if in NF boxes and
! sufficiently well-separated.

   FUNCTION fmm_test_NF_ext(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: fmm_test_NF_ext
      TYPE(box_mm_paras) :: LHStmp, RHStmp

      LHStmp%box = LHS%raw_paras(id%LHS)%box
      RHStmp%box = RHS%raw_paras(id%RHS)%box
      LHStmp%bra = LHS%raw_paras(id%LHS)%bra
      RHStmp%bra = RHS%raw_paras(id%RHS)%bra
      LHStmp%level = 0
      RHStmp%level = 0

      fmm_test_NF_ext = fmm_NF_boxes(LHStmp,RHStmp) .AND.   &
                        fmm_test_ext(LHS,RHS,id)

   END FUNCTION fmm_test_NF_ext

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Far Field"
! based on separation.  Used in simple full quadratic (FQ) algorithm
! when we wish to compare unboxed moments (c.f. fmm_test_FF function)

   FUNCTION fmm_test_raw_FF(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL            :: fmm_test_raw_FF
      TYPE(box_mm_paras) :: LHStmp, RHStmp

      LHStmp%box = LHS%raw_paras(id%LHS)%box
      RHStmp%box = RHS%raw_paras(id%RHS)%box
      LHStmp%bra = LHS%raw_paras(id%LHS)%bra
      RHStmp%bra = RHS%raw_paras(id%RHS)%bra
      LHStmp%level = 0
      RHStmp%level = 0

      fmm_test_raw_FF = .NOT. (fmm_NF_boxes(LHStmp,RHStmp))

   END FUNCTION fmm_test_raw_FF

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Far Field"
! based on separation.  Used by Boxed Quadratic (BQ) algorithm.

   FUNCTION fmm_test_FF(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: fmm_test_FF

      fmm_test_FF = .NOT.(fmm_NF_boxes(LHS%box_paras(id%LHS),   &
                                       RHS%box_paras(id%RHS)))

   END FUNCTION fmm_test_FF

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Local Far Field"
! based on separation;  used by FMM and NlogN algorithms;
! note that in the NlogN scheme boxes are at different levels, but the
! definition of LFF is such that we can then translate the box at the
! deeper level to the higher level and apply the same criteria as for
! two boxes at the deeper level (assuming branches join at higher levels)
! FIXME: how about branch scheme where branches don't join at higher levels??

   FUNCTION fmm_test_LFF(LHS_in,RHS_in,id)

      USE fmm_box_utils, ONLY: fmm_translate_to_common_grid

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS_in, RHS_in
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: fmm_test_LFF

      ! introduce these so we can translate to common grid if needed
      TYPE(box_mm_paras) :: LHS, RHS

      LHS = LHS_in%box_paras(id%LHS)
      RHS = RHS_in%box_paras(id%RHS)
      IF (LHS%level /= RHS%level) THEN
         ! All NF, LFF and RFF definitions are based on boxes in ONE grid.
         ! We must therefore translate the boxes into a common grid first.
         ! Of course, this can be avoided by passing the transformed paras.
         CALL fmm_translate_to_common_grid(LHS,RHS)
      END IF

      fmm_test_LFF = .FALSE.
      IF (fmm_NF_boxes(LHS,RHS)) RETURN
      IF (fmm_RFF_boxes(LHS,RHS)) RETURN
      fmm_test_LFF = .TRUE.

   END FUNCTION fmm_test_LFF

!-------------------------------------------------------------------------------

END MODULE fmm_T_pair_tests

