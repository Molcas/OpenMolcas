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
MODULE fmm_Vff_driver

   USE fmm_stats
   USE fmm_global_paras
   USE fmm_T_contractors,  ONLY: fmm_set_T_con_ptrs,         &
                                 fmm_select_T_con,           &
                                 fmm_init_T_contractors,     &
                                 fmm_free_T_contractors
   USE fmm_W_contractors,  ONLY: fmm_init_W_contractors,     &
                                 fmm_free_W_contractors
   USE fmm_T_pair_builder, ONLY: fmm_init_T_pair_builder,    &
                                 fmm_close_T_pair_builder,   &
                                 fmm_gen_local_T_pairs,      &
                                 fmm_gen_nonlocal_T_pairs

   USE fmm_box_utils,      ONLY: fmm_deepest_level
   USE fmm_box_builder,    ONLY: fmm_get_box_qlm_at_level,   &
                                 fmm_get_box_paras_at_level
   USE fmm_W_pair_builder, ONLY: fmm_get_raw_Vff_from_boxed_Vff

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_Vff

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_verify_data_in(RHS_in)

      IMPLICIT NONE
      TYPE(raw_mm_data), INTENT(IN) :: RHS_in
      LOGICAL :: fail

      fail = .FALSE.
      IF (.NOT.ASSOCIATED(RHS_in%paras)) fail = .TRUE.
      IF (.NOT.ASSOCIATED(RHS_in%qlm_T)) fail = .TRUE.
      IF (.NOT.ASSOCIATED(RHS_in%qlm_W)) fail = .TRUE.
      IF (SIZE(RHS_in%paras) /= SIZE(RHS_in%qlm_T,2)) fail = .TRUE.
      IF (SIZE(RHS_in%paras) /= SIZE(RHS_in%qlm_W,2)) fail = .TRUE.
      IF (fail) CALL fmm_quit('mms pointers sent incorrectly to fmm_Vff_driver')

   END SUBROUTINE fmm_verify_data_in

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_Vff(scheme,LHS_paras,RHS_mms,Vff)

      USE fmm_box_builder,  ONLY: fmm_init_box_builder, fmm_free_box_builder
      USE fmm_local_search, ONLY: fmm_init_local_search
      USE fmm_local_search, ONLY: fmm_free_local_search

      IMPLICIT NONE
      TYPE(raw_mm_paras),  INTENT(INOUT) :: LHS_paras(:)
      TYPE(raw_mm_data),   INTENT(INOUT) :: RHS_mms
      TYPE(scheme_paras),  INTENT(INOUT) :: scheme
      REAL(REALK), TARGET, INTENT(OUT)   :: Vff(:,:)

      WRITE(LUPRI,*) "Computing multipole potential..."
      CALL fmm_verify_data_in(RHS_mms)

      CALL fmm_init_box_builder(LHS_paras,RHS_mms,scheme)
      IF (scheme%algorithm == DO_FMM) CALL fmm_init_local_search(scheme)

      scheme%phase = FAR_FIELD
      stat_NF_not_FF = .FALSE.
      SELECT CASE (scheme%algorithm)
         CASE (DO_NULL)
            CONTINUE
         CASE (DO_FQ)
            CALL fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)
         CASE (DO_BQ)
            CALL fmm_get_BQ_Vff(scheme,LHS_paras,Vff)
         CASE (DO_NlogN)
            CALL fmm_get_NlogN_Vff(scheme,LHS_paras,Vff)
         CASE (DO_FMM)
            CALL fmm_get_FMM_Vff(scheme,LHS_paras,Vff)
         CASE DEFAULT
            CALL fmm_quit('invalid algorithm requested!')
      END SELECT

      IF (scheme%include_near_field) THEN
         ! Also compute multipole interactions within NF boxes
         scheme%phase = NEAR_FIELD
         stat_NF_not_FF = .TRUE.
         CALL fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)
      END IF

      CALL fmm_free_box_builder
      CALL fmm_free_local_search

   END SUBROUTINE fmm_get_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_paras(:)
      TYPE(raw_mm_data),          INTENT(IN)    :: RHS_mms
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      INTEGER(INTK), PARAMETER :: pair_type = LHS_raw_RHS_raw
      TYPE(gen_mm_paras) :: LHS, RHS
      REAL(REALK) :: T0, fmm_second, TTOTFQ

      T0 = fmm_second()
      NULLIFY(LHS%raw_paras,LHS%box_paras)
      NULLIFY(RHS%raw_paras,RHS%box_paras)

      LHS%raw_paras => LHS_paras(:)
      RHS%raw_paras => RHS_mms%paras(:)

      ! select the T-contractor to be stored/called via C code;
      CALL fmm_select_T_con(scheme)
      ! set T_contractor pointers
      CALL fmm_set_T_con_ptrs(Vff,RHS_mms%qlm_T)
      CALL fmm_init_T_contractors(scheme)

      ! generate full multipole potential
      CALL fmm_init_T_pair_builder(scheme,pair_type)
      CALL fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
      CALL fmm_close_T_pair_builder(scheme)
      CALL fmm_free_T_contractors

      TTOTFQ = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_FQ_Vff', TTOTFQ, LUPRI)

   END SUBROUTINE fmm_get_FQ_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_BQ_Vff(scheme,LHS_raw_paras,Vff)

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)       :: LHS, RHS
      REAL(REALK), POINTER     :: qlm_T(:,:)
      REAL(REALK), POINTER     :: boxed_Vff(:,:)
      INTEGER(INTK)            :: level, LMAX, lm_dim, foo
      INTEGER(INTK), PARAMETER :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, fmm_second, TTOT

      T0 = fmm_second()
      NULLIFY(qlm_T,boxed_Vff)
      NULLIFY(LHS%raw_paras,LHS%box_paras)
      NULLIFY(RHS%raw_paras,RHS%box_paras)

      CALL fmm_init_W_contractors(scheme%trans_LMAX)

      level = fmm_deepest_level(scheme)
      IF (level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (1+LMAX)**2

      ! set up LHS and RHS boxed parameters
      CALL fmm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
      CALL fmm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
      CALL fmm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','free')

      ! allocate temporary boxed potential
      foo = SIZE(LHS%box_paras)
      ALLOCATE (boxed_Vff(lm_dim,foo))
      boxed_Vff = zero

      ! select the T-contractor to be stored/called via C code
      CALL fmm_select_T_con(scheme)
      ! set T_contractor pointers
      CALL fmm_set_T_con_ptrs(boxed_Vff,qlm_T)
      CALL fmm_init_T_contractors(scheme)

      ! generate FF Vff
      CALL fmm_init_T_pair_builder(scheme,pair_type)
      CALL fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
      CALL fmm_close_T_pair_builder(scheme)
      CALL fmm_free_T_contractors

      ! translate boxed potential to raw LHS centres
      CALL fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
      DEALLOCATE(boxed_Vff)

      TTOT = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_BQ_Vff', TTOT, LUPRI)

      CALL fmm_free_W_contractors

   END SUBROUTINE fmm_get_BQ_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_NlogN_Vff(scheme,LHS_raw_paras,Vff)

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)       :: LHS, RHS
      REAL(REALK), POINTER     :: qlm_T(:,:)
      REAL(REALK), POINTER     :: boxed_Vff(:,:)
      INTEGER(INTK)            :: LHS_level, RHS_level, deepest_level
      INTEGER(INTK)            :: LMAX, lm_dim, foo
      INTEGER(INTK), PARAMETER :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, fmm_second, TTOT

      T0 = fmm_second()
      NULLIFY(qlm_T,boxed_Vff)
      NULLIFY(LHS%raw_paras,LHS%box_paras)
      NULLIFY(RHS%raw_paras,RHS%box_paras)

      CALL fmm_init_W_contractors(scheme%trans_LMAX)

      deepest_level = fmm_deepest_level(scheme)
      IF (deepest_level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (1+LMAX)**2

      ! select the T-contractor to be stored/called via C code
      CALL fmm_select_T_con(scheme)
      CALL fmm_init_T_contractors(scheme)
      ! get LHS boxed MM parameters
      LHS_level = deepest_level
      CALL fmm_get_box_paras_at_level(LHS_level,scheme,LHS%box_paras,'LHS')
      ! allocate boxed potential based on number of LHS boxes
      foo = SIZE(LHS%box_paras)
      ALLOCATE (boxed_Vff(lm_dim,foo))
      boxed_Vff = zero

      ! generate far-field potential at deepest box centres
      DO RHS_level = deepest_level, TOP_LEVEL, -1
         ! get packed RHS MM paras and translated moments
         CALL fmm_get_box_paras_at_level(RHS_level,scheme,RHS%box_paras,'RHS')
         CALL fmm_get_box_qlm_at_level(RHS_level,scheme,qlm_T,'RHS','free')
         ! set T_contractor pointers
         CALL fmm_set_T_con_ptrs(boxed_Vff,qlm_T)
         ! generate LFF potential at prescribed centres
         CALL fmm_init_T_pair_builder(scheme,pair_type)
!FIXME: would be nice to add NlogN to local search algorithm
         CALL fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
         CALL fmm_close_T_pair_builder(scheme)
      END DO
      CALL fmm_free_T_contractors

      ! translate boxed potential to raw LHS centres
      CALL fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
      DEALLOCATE(boxed_Vff)

      TTOT = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_NlogN_Vff', TTOT, LUPRI)

      CALL fmm_free_W_contractors

   END SUBROUTINE fmm_get_NlogN_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_FMM_Vff(scheme,LHS_raw_paras,Vff)

      USE fmm_W_pair_builder, ONLY: fmm_translate_parents_Vff

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),                INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)          :: LHS, RHS
      TYPE(box_mm_paras), POINTER :: p_box_paras(:)   ! parents data
      REAL(REALK), POINTER        :: qlm_T(:,:)
      REAL(REALK), POINTER        :: Vff_tmp1(:,:)
      REAL(REALK), POINTER        :: Vff_tmp2(:,:)
      REAL(REALK), POINTER        :: Vff_ptr(:,:)
      INTEGER(INTK)               :: level, l_up, deepest_level
      INTEGER(INTK)               :: LMAX, lm_dim, foo
      INTEGER(INTK)               :: use_Vff_tmp
      INTEGER(INTK), PARAMETER    :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, fmm_second, TTOT

      T0 = fmm_second()
      NULLIFY(qlm_T,Vff_tmp1,Vff_tmp2,Vff_ptr)
      NULLIFY(LHS%raw_paras,LHS%box_paras)
      NULLIFY(RHS%raw_paras,RHS%box_paras)

      CALL fmm_init_W_contractors(scheme%trans_LMAX)

      use_Vff_tmp = 0
      deepest_level = fmm_deepest_level(scheme)
      IF (deepest_level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (LMAX+1)**2

      ! select the T-contractor to be stored/called via C code
      CALL fmm_select_T_con(scheme)
      CALL fmm_init_T_contractors(scheme)
      ! generate far-field potential using whole box hierarchy
      DO level = TOP_LEVEL, deepest_level
         l_up  = level-1
         ! initialise T-contractor and allocate T-matrix
         ! get LHS boxed MM parameters
         CALL fmm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
         ! get packed RHS MM paras and translated moments
         CALL fmm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
         CALL fmm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','keep')
         ! allocate temporary boxed potentials
         IF (use_Vff_tmp /= 2) THEN
            foo = SIZE(LHS%box_paras)
            ALLOCATE (Vff_tmp1(lm_dim,foo))
            Vff_tmp1 = zero
            Vff_ptr => Vff_tmp1(:,:)
         ELSE
            foo = SIZE(LHS%box_paras)
            ALLOCATE (Vff_tmp2(lm_dim,foo))
            Vff_tmp2 = zero
            Vff_ptr => Vff_tmp2(:,:)
         END IF
         ! set T_contractor pointers
         CALL fmm_set_T_con_ptrs(Vff_ptr,qlm_T)
         ! generate LFF contribution to Vff at this level
         CALL fmm_init_T_pair_builder(scheme,pair_type)
         CALL fmm_gen_local_T_pairs(LHS,RHS,pair_type)
         CALL fmm_close_T_pair_builder(scheme)
         ! get RFF contribution from parent level
         SELECT CASE (use_Vff_tmp)
            CASE (0)
               use_Vff_tmp = 2
            CASE (1)
               ! RFF contribution: parents' translated FF ( parent | child )
               CALL fmm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')
               CALL fmm_translate_parents_Vff(level,scheme,Vff_tmp2,Vff_tmp1, &
                                              LHS%box_paras)
               ! deallocate parent's space
               DEALLOCATE(Vff_tmp2)
               NULLIFY(Vff_tmp2)
               use_Vff_tmp = 2
            CASE (2)
               ! RFF contribution: parents' translated FF ( parent | child )
               CALL fmm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')
               CALL fmm_translate_parents_Vff(level,scheme,Vff_tmp1,Vff_tmp2, &
                                              LHS%box_paras)
               ! deallocate parent's space
               DEALLOCATE(Vff_tmp1)
               NULLIFY(Vff_tmp1)
               use_Vff_tmp = 1
         END SELECT
      END DO
      CALL fmm_free_T_contractors

      ! translate boxed potential to raw LHS centres for final contraction
      IF (use_Vff_tmp == 1) THEN  ! Vff_tmp2 holds final boxed potential
         CALL fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp2,Vff)
      ELSE  ! Vff_tmp1 holds final boxed potential
         CALL fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp1,Vff)
      END IF
      IF (ASSOCIATED(Vff_tmp1)) DEALLOCATE(Vff_tmp1)
      IF (ASSOCIATED(Vff_tmp2)) DEALLOCATE(Vff_tmp2)

      TTOT = fmm_second()-T0
      CALL TIMTXT('>>> TIME USED in fmm_get_FMM_Vff', TTOT, LUPRI)

      CALL fmm_free_W_contractors

   END SUBROUTINE fmm_get_FMM_Vff

!-------------------------------------------------------------------------------

END MODULE fmm_Vff_driver

