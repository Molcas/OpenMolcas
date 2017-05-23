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
MODULE fmm_T_pair_builder

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_T_pair_builder,  &
             fmm_close_T_pair_builder, &
             fmm_gen_local_T_pairs,    &
             fmm_gen_nonlocal_T_pairs

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_T_pair_builder(scheme,pair_type)

      USE fmm_T_pair_mould,  ONLY: fmm_init_T_pair_mould
      USE fmm_T_pair_tests,  ONLY: fmm_init_T_pair_tests
      USE fmm_T_buffer,      ONLY: fmm_open_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK), INTENT(IN) :: pair_type

      CALL fmm_init_T_pair_mould(scheme,pair_type)
      CALL fmm_init_T_pair_tests(scheme)
      CALL fmm_open_T_buffer(scheme)

   END SUBROUTINE fmm_init_T_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_close_T_pair_builder(scheme)

      USE fmm_T_pair_mould,  ONLY: fmm_close_T_pair_mould
      USE fmm_T_pair_tests,  ONLY: fmm_close_T_pair_tests
      USE fmm_T_buffer,      ONLY: fmm_close_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme

      CALL fmm_close_T_pair_mould
      CALL fmm_close_T_pair_tests
      CALL fmm_close_T_buffer(scheme)

   END SUBROUTINE fmm_close_T_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_gen_local_T_pairs(LHS,RHS,pair_type)

      USE fmm_T_pair_tests, ONLY: fmm_test_and_buffer_T_pair
      USE fmm_local_search, ONLY: fmm_get_local_paras

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(INOUT) :: LHS, RHS
      INTEGER(INTK),      INTENT(IN)    :: pair_type
      TYPE(gen_mm_paras) :: RHS_local
      TYPE(LHS_RHS_type) :: id
      INTEGER(INTK)      :: weight, i,j, ndim

      weight = 1
      NULLIFY (RHS_local%box_paras, RHS_local%raw_paras)

      SELECT CASE (pair_type)
         CASE (LHS_raw_RHS_raw)
            DO j = 1, SIZE(LHS%raw_paras)
               id%LHS = j
               DO i = 1, SIZE(RHS%raw_paras)
                  id%RHS = i
                  CALL fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
               END DO
            END DO
         CASE (LHS_box_RHS_box)
            DO j = 1, SIZE(LHS%box_paras)
               CALL fmm_get_local_paras(j,RHS,pair_type,RHS_local,ndim)
               IF (ndim == 0) CYCLE
               id%LHS = j
               DO i = 1, SIZE(RHS_local%box_paras)
                  id%RHS = i
                  CALL fmm_test_and_buffer_T_pair(LHS,RHS_local,id,weight)
               END DO
               IF (ASSOCIATED(RHS_local%box_paras)) THEN
                  DEALLOCATE (RHS_local%box_paras)
                  NULLIFY (RHS_local%box_paras)
               END IF
            END DO
         CASE DEFAULT
            CALL fmm_quit ('cannot reconcile requested T_pair type!')
      END SELECT

   END SUBROUTINE fmm_gen_local_T_pairs

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)

      USE fmm_T_pair_tests, ONLY: fmm_test_and_buffer_T_pair

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(INOUT) :: LHS, RHS
      INTEGER(INTK),      INTENT(IN)    :: pair_type
      TYPE(LHS_RHS_type) :: id
      INTEGER(INTK)      :: weight, i,j

      weight = 1

      SELECT CASE (pair_type)
         CASE (LHS_raw_RHS_raw)
            DO j = 1, SIZE(LHS%raw_paras)
               id%LHS = j
               DO i = 1, SIZE(RHS%raw_paras)
                  id%RHS = i
                  CALL fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
               END DO
            END DO
         CASE (LHS_box_RHS_box)
            DO j = 1, SIZE(LHS%box_paras)
               id%LHS = j
               DO i = 1, SIZE(RHS%box_paras)
                  id%RHS = i
                  CALL fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
               END DO
            END DO
         CASE DEFAULT
            CALL fmm_quit ('cannot reconcile requested T_pair type!')
      END SELECT

   END SUBROUTINE fmm_gen_nonlocal_T_pairs

!-------------------------------------------------------------------------------

END MODULE fmm_T_pair_builder
