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
MODULE fmm_T_pair_mould

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_T_pair_mould,           &
             fmm_close_T_pair_mould

   INTEGER(INTK), SAVE :: LHS_LMAX, RHS_LMAX
    ! flag to test initialisation
   CHARACTER(LEN=11), SAVE :: fmm_init_mould

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_close_T_pair_mould

      IMPLICIT NONE
      IF (fmm_init_mould /= 'initialised') STOP 'mm_T_pair_mould init'
      fmm_init_mould = ' '
      LHS_LMAX = 0
      RHS_LMAX = 0

   END SUBROUTINE fmm_close_T_pair_mould

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_T_pair_basics(LHS,RHS,id,weight,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      INTEGER(INTK),       INTENT(IN)  :: weight
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%N_or_T = 'N' ! not used for T contractions, so null
      ! init. normalisation scalar for other options
      T_pair%paras%ratio = one
      ! here, T_pair%r_ab is the unnormalised vector

      ! weighting to account for double-counting in T-pair generation
      T_pair%paras%weight = weight

      T_pair%LMAX = MAX(T_pair%paras%LHS_LMAX,T_pair%paras%RHS_LMAX)
      T_pair%lm_max = (1+T_pair%LMAX)**2

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer_array(LHS%raw_paras%id)
         CALL Unused_integer_array(RHS%raw_paras%id)
         CALL Unused_integer(id%LHS)
      END IF
   END SUBROUTINE fmm_set_T_pair_basics

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_RR_paras(LHS,RHS,id,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      ! interaction vector for building T-matrix
      T_pair%r_ab = RHS%raw_paras(id%RHS)%cntr - LHS%raw_paras(id%LHS)%cntr
      ! indices to map back to actual moments in separate array
      T_pair%paras%LHS_id = LHS%raw_paras(id%LHS)%id
      T_pair%paras%RHS_id = RHS%raw_paras(id%RHS)%id
      ! check that paras:moments mapping was built
      IF (T_pair%paras%LHS_id == 0) CALL fmm_quit('LHS paras:moments mapping')
      IF (T_pair%paras%RHS_id == 0) CALL fmm_quit('RHS paras:moments mapping')

   END SUBROUTINE fmm_set_RR_paras

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_BB_paras(LHS,RHS,id,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      ! interaction vector for building T-matrix
      T_pair%r_ab = RHS%box_paras(id%RHS)%cntr - LHS%box_paras(id%LHS)%cntr
      ! indices to map back to actual moments in separate array
      T_pair%paras%LHS_id = LHS%box_paras(id%LHS)%id
      T_pair%paras%RHS_id = RHS%box_paras(id%RHS)%id
      ! check that paras:moments mapping was built
      IF (T_pair%paras%LHS_id == 0) CALL fmm_quit('LHS paras:moments mapping')
      IF (T_pair%paras%RHS_id == 0) CALL fmm_quit('RHS paras:moments mapping')

   END SUBROUTINE fmm_set_BB_paras

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_LHS_LMAX(X,Y,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: X        ! dummy variables
      TYPE(LHS_RHS_type),  INTENT(IN)  :: Y        ! dummy variables
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%LHS_LMAX = LHS_LMAX

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer_array(X%raw_paras%id)
         CALL Unused_integer(Y%LHS)
      END IF
   END SUBROUTINE fmm_set_LHS_LMAX

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_set_RHS_LMAX(X,Y,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: X        ! dummy variables
      TYPE(LHS_RHS_type),  INTENT(IN)  :: Y        ! dummy variables
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%RHS_LMAX = RHS_LMAX

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer_array(X%raw_paras%id)
         CALL Unused_integer(Y%LHS)
      END IF
   END SUBROUTINE fmm_set_RHS_LMAX

!-------------------------------------------------------------------------------
! This routine directs the saving of functions in fmm_proc_selector.c
! relevant to the formation of a T-pair entity.
! There are 4 parts to the making of a T-pair once raw data has been
! suppplied from the classical interaction search algorithm.
! These are all called consecutively from the C-code via
! fmm_stored_t_pair_mould.
! The 4 functions saved are selected here and depend on whether
! pure boxed moments, or unboxed moments, are being contracted.

   SUBROUTINE fmm_init_T_pair_mould(scheme,pair_type)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),      INTENT(IN) :: pair_type
      EXTERNAL fmm_store_t_pair_mould1   ! raw/box dependent variables
      EXTERNAL fmm_store_t_pair_mould2   ! set LMAX (LHS)
      EXTERNAL fmm_store_t_pair_mould3   ! set LMAX (RHS)
      EXTERNAL fmm_store_t_pair_mould4   ! common to all T-pair builds

      CALL fmm_store_t_pair_mould2(fmm_set_LHS_LMAX)
      CALL fmm_store_t_pair_mould3(fmm_set_RHS_LMAX)
      CALL fmm_store_t_pair_mould4(fmm_set_T_pair_basics)
      SELECT CASE (pair_type)
      CASE (LHS_raw_RHS_raw)
         LHS_LMAX = scheme%raw_LMAX
         RHS_LMAX = scheme%raw_LMAX
         CALL fmm_store_t_pair_mould1(fmm_set_RR_paras)
      CASE (LHS_box_RHS_box)
         LHS_LMAX = scheme%trans_LMAX
         RHS_LMAX = scheme%trans_LMAX
         CALL fmm_store_t_pair_mould1(fmm_set_BB_paras)
      CASE DEFAULT
         CALL fmm_quit ('cannot recognise T_pair type!')
      END SELECT

      fmm_init_mould = 'initialised'

   END SUBROUTINE fmm_init_T_pair_mould

!-------------------------------------------------------------------------------

END MODULE fmm_T_pair_mould
