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
MODULE fmm_box_builder

   USE fmm_global_paras
   USE fmm_stats
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_box_builder,                          &
             fmm_free_box_builder,                          &
             fmm_get_box_paras_at_level,                    &
             fmm_get_box_qlm_at_level

   ! Pointers to unboxed LHS parameters and RHS moments & parameters
   TYPE(raw_mm_data),  POINTER, SAVE :: RHS_raw_mms
   TYPE(raw_mm_paras), POINTER, SAVE :: LHS_raw_paras(:)
   ! Packed LHS & RHS parameters and RHS moments at all levels
   TYPE(box_mm_data), POINTER, SAVE :: mms_at_lev(:)
   ! Deepest level of boxes in the hierarchy
   INTEGER(INTK), SAVE :: deepest_level

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_box_builder(LHS_paras,RHS_mms,scheme)

      USE fmm_box_utils, ONLY: fmm_deepest_level

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(INOUT) :: LHS_paras(:)
      TYPE(raw_mm_data),  TARGET, INTENT(INOUT) :: RHS_mms
      TYPE(scheme_paras),         INTENT(IN)    :: scheme

      deepest_level = fmm_deepest_level(scheme)
      stat_deepest_level = deepest_level

      ! initialise RHS pointers to raw moments and parameters
      RHS_raw_mms => RHS_mms
      ! initialise LHS pointers to raw parameters
      LHS_raw_paras => LHS_paras(:)
      ! allocate foundation of hierarchical moment data structure
      CALL preallocate_levels
      ! update raw_mm_paras with box and branch data
      CALL init_box_paras(LHS_paras,RHS_mms%paras,scheme)

   END SUBROUTINE fmm_init_box_builder

!-------------------------------------------------------------------------------

   SUBROUTINE init_box_paras(LHS,RHS,scheme)

      USE fmm_box_utils, ONLY: fmm_box, fmm_branch, fmm_grain, fmm_box_centre

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: LHS(:)
      TYPE(raw_mm_paras), INTENT(INOUT) :: RHS(:)
      TYPE(scheme_paras), INTENT(IN)    :: scheme

      INTEGER(INTK) :: i
      REAL(REALK)   :: grain, grain_inv

      grain     = fmm_grain(scheme,deepest_level)
      grain_inv = one/grain

      DO i = 1, SIZE(RHS)
         RHS(i)%box = fmm_box(RHS(i)%cntr,grain_inv)
         RHS(i)%bra = fmm_branch(RHS(i)%ext,grain_inv)
         RHS(i)%box_cntr = fmm_box_centre(RHS(i)%box,grain)
         RHS(i)%map_up = 0     ! not defined yet
      END DO

      DO i = 1, SIZE(LHS)
         LHS(i)%box = fmm_box(LHS(i)%cntr,grain_inv)
         LHS(i)%bra = fmm_branch(LHS(i)%ext,grain_inv)
         LHS(i)%box_cntr = fmm_box_centre(LHS(i)%box,grain)
         LHS(i)%map_up = 0     ! not defined yet
      END DO

   END SUBROUTINE init_box_paras

!-------------------------------------------------------------------------------

   SUBROUTINE build_paras_at_level(level,scheme)

      USE fmm_box_packer, ONLY: fmm_init_pkd_paras

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),      INTENT(IN) :: level

      TYPE(raw_mm_paras), POINTER :: ptr(:)

      IF ((level < TOP_LEVEL) .OR. (level > deepest_level)) THEN
         CALL fmm_quit ('cannot iterate paras to this level!')
      END IF

      ! build packed parameters at deepest level if not already done
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%RHS_paras)) THEN
         ptr => RHS_raw_mms%paras(:)
         CALL fmm_init_pkd_paras(deepest_level,scheme,ptr,   &
                                 mms_at_lev(deepest_level)%RHS_paras)
      END IF
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%LHS_paras)) THEN
         ptr => LHS_raw_paras(:)
         CALL fmm_init_pkd_paras(deepest_level,scheme,ptr,  &
                                 mms_at_lev(deepest_level)%LHS_paras)
      END IF

      IF (level < deepest_level) THEN
         ! iterate paras through the box hierarchy to the reqd level
         CALL iterate_paras_to_level(level,scheme,'RHS')
         CALL iterate_paras_to_level(level,scheme,'LHS')
      END IF

   END SUBROUTINE build_paras_at_level

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE iterate_paras_to_level(l,scheme,side)

      USE fmm_box_packer, ONLY: fmm_shift_and_pack_paras

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN) :: l         ! level in box hierarchy
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CHARACTER(3),       INTENT(IN) :: side

      TYPE(box_mm_paras), POINTER :: ptr(:)
      INTEGER(INTK) :: l_down

      l_down = l+1
      SELECT CASE (side)
      CASE ('RHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%RHS_paras)) THEN
            ! must build paras at deeper levels before higher levels
            CALL iterate_paras_to_level(l_down,scheme,'RHS')
         END IF
         ptr => mms_at_lev(l_down)%RHS_paras(:)
         ! build new paras from paras at previous level (including packing)
         CALL fmm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%RHS_paras)
      CASE ('LHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%LHS_paras)) THEN
            ! must build paras at deeper levels before higher levels
            CALL iterate_paras_to_level(l_down,scheme,'LHS')
         END IF
          ptr => mms_at_lev(l_down)%LHS_paras(:)
          ! build new paras from paras at previous level (including packing)
          CALL fmm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%LHS_paras)
      CASE DEFAULT
         CALL fmm_quit ('must build LHS or RHS paras!')
      END SELECT

   END SUBROUTINE iterate_paras_to_level

!-------------------------------------------------------------------------------

   SUBROUTINE build_qlm_at_level(level,scheme,memory)

      USE fmm_W_pair_builder, ONLY: fmm_translate_raw_moments

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),      INTENT(IN) :: level
      CHARACTER(4),       INTENT(IN) :: memory

      TYPE(raw_mm_data), POINTER :: ptr1
      TYPE(box_mm_data), POINTER :: ptr2
      INTEGER(INTK) :: mms_dim

      IF ((level < TOP_LEVEL) .OR. (level > deepest_level)) THEN
         CALL fmm_quit ('cannot iterate boxed moments to this level!')
      END IF

      ! build boxed moments at the deepest level if not already done
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%qlm_W)) THEN
         mms_dim = SIZE(mms_at_lev(deepest_level)%RHS_paras)
         CALL allocate_lm_at_level(deepest_level,mms_dim,scheme%trans_LMAX)
         IF (.NOT.ASSOCIATED(RHS_raw_mms))    &
               CALL fmm_quit('mm_box_builder not correctly initialised!')
         ptr1 => RHS_raw_mms
         ptr2 => mms_at_lev(deepest_level)
         CALL fmm_translate_raw_moments(scheme,ptr1,ptr2)
      END IF

      IF (level < deepest_level) THEN
         ! iterate RHS moments through the box hierarchy to the reqd level
         CALL iterate_qlm_to_level(level,scheme,memory)
      END IF

   END SUBROUTINE build_qlm_at_level

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE iterate_qlm_to_level(level,scheme,memory)

      USE fmm_W_pair_builder, ONLY: fmm_translate_boxed_moments

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN) :: level
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CHARACTER(4),       INTENT(IN) :: memory

      INTEGER(INTK) :: mms_dim, l_down
      TYPE(box_mm_data), POINTER :: ptr1, ptr2

      l_down = level+1
      ! note we always build qlm_W first, then scale to give qlm_T
      IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%qlm_W)) THEN
         ! must have qlm at deeper levels before higher levels
         CALL iterate_qlm_to_level(l_down,scheme,memory)
      END IF

      ! must build boxed paras before boxed moments
      IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%RHS_paras)) THEN
         CALL build_paras_at_level(l_down,scheme)
      END IF

      mms_dim = SIZE(mms_at_lev(level)%RHS_paras)
      CALL allocate_lm_at_level(level,mms_dim,scheme%trans_LMAX)
      ptr1 => mms_at_lev(l_down)
      ptr2 => mms_at_lev(level)
      CALL fmm_translate_boxed_moments(scheme,ptr1,ptr2)

   END SUBROUTINE iterate_qlm_to_level

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_box_paras_at_level(l,scheme,box_paras,side)

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN) :: l
      TYPE(scheme_paras), INTENT(IN) :: scheme
      TYPE(box_mm_paras), POINTER    :: box_paras(:)
      CHARACTER(3),       INTENT(IN) :: side

      IF (.NOT.ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should be allocated!'

      SELECT CASE (side)
      CASE ('RHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l)%RHS_paras)) THEN
            ! RHS paras not available at this level yet. So make them!
            CALL build_paras_at_level(l,scheme)
         END IF
         box_paras => mms_at_lev(l)%RHS_paras(:)
         stat_RHS_boxes(l) = SIZE(box_paras)
      CASE ('LHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l)%LHS_paras)) THEN
            ! LHS paras not available at this level yet. So make them!
            CALL build_paras_at_level(l,scheme)
         END IF
         box_paras => mms_at_lev(l)%LHS_paras(:)
         stat_LHS_boxes(l) = SIZE(box_paras)
      CASE DEFAULT
         CALL fmm_quit ('must select just LHS or RHS paras to use')
      END SELECT

   END SUBROUTINE fmm_get_box_paras_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_box_qlm_at_level(l,scheme,qlm_T,side,memory)

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)  :: l
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      REAL(REALK),        POINTER     :: qlm_T(:,:)
      CHARACTER(3),       INTENT(IN)  :: side
      CHARACTER(4),       INTENT(IN)  :: memory

      IF (.NOT.ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should be allocated!'
      IF (.NOT.ASSOCIATED(mms_at_lev(l)%qlm_T)) THEN
         ! qlm data not available at this level yet. So make them!
         CALL build_qlm_at_level(l,scheme,memory)
      END IF

      IF (side == 'LHS') CALL fmm_quit('currently no LHS boxed mms built!')
      IF (side == 'RHS') THEN
         qlm_T => mms_at_lev(l)%qlm_T(:,:)
      ELSE
         CALL fmm_quit("must select LHS or RHS boxed moments!")
      END IF

   END SUBROUTINE fmm_get_box_qlm_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE preallocate_levels

      IMPLICIT NONE
      INTEGER(INTK) :: i

      IF (deepest_level == 0) RETURN
      IF (ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should not be allocated!'
      IF (deepest_level < TOP_LEVEL) THEN
         CALL fmm_quit('error allocating levels in box hierarchy')
      END IF
      ! We allow for possibility of all levels being used, but sub-variables
      ! are only allocated if required, so this is not a big deal
      ALLOCATE (mms_at_lev(deepest_level))
      DO i = LBOUND(mms_at_lev,1), UBOUND(mms_at_lev,1)
         NULLIFY (mms_at_lev(i)%LHS_paras)
         NULLIFY (mms_at_lev(i)%RHS_paras)
         NULLIFY (mms_at_lev(i)%qlm_W)
         NULLIFY (mms_at_lev(i)%qlm_T)
      END DO

   END SUBROUTINE preallocate_levels

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_lm_at_level(l,mms_dim,LMAX)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: l, mms_dim, LMAX

      INTEGER(INTK) :: lm_dim
      LOGICAL :: fail

      lm_dim = (1+LMAX)**2
      IF (l > deepest_level) CALL fmm_quit('invalid level to allocate!')
      IF (l < TOP_LEVEL) CALL fmm_quit('invalid level to allocate!')

      ! must test if pointers already allocated (compiler may not notice)
      fail = .FALSE.
      IF (ASSOCIATED(mms_at_lev(l)%qlm_T)) fail = .TRUE.
      IF (ASSOCIATED(mms_at_lev(l)%qlm_W)) fail = .TRUE.
      IF (fail) CALL fmm_quit('box lm data already allocated!')
      ALLOCATE (mms_at_lev(l)%qlm_T(lm_dim,mms_dim))
      ALLOCATE (mms_at_lev(l)%qlm_W(lm_dim,mms_dim))
      ! must zero as they are built additively
      mms_at_lev(l)%qlm_T = zero
      mms_at_lev(l)%qlm_W = zero

   END SUBROUTINE allocate_lm_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_box_builder

      IMPLICIT NONE
      INTEGER(INTK) :: l

      NULLIFY (RHS_raw_mms)
!FIXME
      NULLIFY (LHS_raw_paras)

      IF (ASSOCIATED(mms_at_lev)) THEN
        DO l = LBOUND(mms_at_lev,1), UBOUND(mms_at_lev,1)
          IF (ASSOCIATED(mms_at_lev(l)%LHS_paras,mms_at_lev(l)%RHS_paras)) THEN
             ! LHS and RHS paras the same and only point the same space
             DEALLOCATE(mms_at_lev(l)%RHS_paras)
          ELSE
             IF (ASSOCIATED(mms_at_lev(l)%RHS_paras))    &
                 DEALLOCATE(mms_at_lev(l)%RHS_paras)
             IF (ASSOCIATED(mms_at_lev(l)%LHS_paras))    &
                 DEALLOCATE(mms_at_lev(l)%LHS_paras)
          END IF
          IF (ASSOCIATED(mms_at_lev(l)%qlm_W)) DEALLOCATE(mms_at_lev(l)%qlm_W)
          IF (ASSOCIATED(mms_at_lev(l)%qlm_T)) DEALLOCATE(mms_at_lev(l)%qlm_T)
          NULLIFY (mms_at_lev(l)%LHS_paras)
          NULLIFY (mms_at_lev(l)%RHS_paras)
          NULLIFY (mms_at_lev(l)%qlm_W)
          NULLIFY (mms_at_lev(l)%qlm_T)
        END DO
        DEALLOCATE (mms_at_lev)
      END IF
      deepest_level = 0

   END SUBROUTINE fmm_free_box_builder

!-------------------------------------------------------------------------------

END MODULE fmm_box_builder

