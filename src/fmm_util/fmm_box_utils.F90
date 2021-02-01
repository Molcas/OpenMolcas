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
MODULE fmm_box_utils

   USE fmm_global_paras, ONLY: INTK, REALK, scheme_paras, gen_mm_paras, LHS_RHS_TYPE, box_mm_paras, MAX_LEVEL, TOP_LEVEL, WS_MIN, &
                               Half, Two
   USE fmm_stats, ONLY: stat_level_saturation, stat_max_branch, stat_min_branch
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_box,                       &
             fmm_box_centre,                &
             fmm_branch,                    &
             fmm_grain,                     &
             fmm_parent_box,                &
             fmm_parent_bra,                &
             fmm_deepest_level,             &
             fmm_def_WS_NF,                 &
             fmm_def_WS_RFF,                &
             fmm_same_box,                  &
             fmm_RFF_boxes,                 &
             fmm_NF_boxes,                  &
             fmm_translate_to_common_grid

CONTAINS

!-------------------------------------------------------------------------------

   FUNCTION fmm_parent_box(box)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: box(3)
      INTEGER(INTK) :: fmm_parent_box(3)

      fmm_parent_box = 1+((box-1)/2)    ! get largest integer after /2

   END FUNCTION fmm_parent_box

!-------------------------------------------------------------------------------

   FUNCTION fmm_parent_bra(branch)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: branch
      INTEGER(INTK) :: fmm_parent_bra

!      fmm_parent_bra = (((branch-1)/2)/2)*2 +2
!      fmm_parent_bra = MAX(WS_MIN, fmm_parent_bra)
      fmm_parent_bra = WS_MIN

! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(branch)
   END FUNCTION fmm_parent_bra

!-------------------------------------------------------------------------------
! We use the coordinate shift to ensure all the box indices are positive

   FUNCTION fmm_box(centre,grain_inv)

      USE fmm_qlm_builder, ONLY: fmm_coord_shift

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN) :: centre(3)
      REAL(REALK),   INTENT(IN) :: grain_inv
      INTEGER(INTK) :: fmm_box(3)

!      fmm_box(:) = 1+ INT(centre(:)*grain_inv)
      fmm_box(:) = 1+ INT( (centre(:)-fmm_coord_shift(:)) *grain_inv )

   END FUNCTION fmm_box

!-------------------------------------------------------------------------------
! We take account of the coordinate shift used in fmm_box

   FUNCTION fmm_box_centre(box,grain)

      USE fmm_qlm_builder, ONLY: fmm_coord_shift

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: box(3)
      REAL(REALK),   INTENT(IN) :: grain
      REAL(REALK) :: fmm_box_centre(3)

!      fmm_box_centre(:) = grain*(box(:)-half)
      fmm_box_centre(:) = grain*(box(:)-half) + fmm_coord_shift(:)

   END FUNCTION fmm_box_centre

!-------------------------------------------------------------------------------

   FUNCTION fmm_branch(extent,grain_inv)

      IMPLICIT NONE
!      TYPE(scheme_paras), INTENT(IN) :: scheme
      REAL(REALK), INTENT(IN) :: extent, grain_inv
      INTEGER(INTK) :: fmm_branch

      ! Note branches are always even in CFMM algorithm
!      IF (scheme%branch_free) THEN
!         fmm_branch = WS_MIN
!      ELSE
!         fmm_branch = 2 * CEILING(extent*grain_inv)
!      END IF

!      fmm_branch = 2 * CEILING(extent*grain_inv)
!      fmm_branch = MAX(WS_MIN, fmm_branch)

      fmm_branch = WS_MIN

      stat_max_branch = MAX(stat_max_branch,fmm_branch)
      stat_min_branch = MIN(stat_min_branch,fmm_branch)

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_real(extent)
         CALL Unused_real(grain_inv)
      END IF
   END FUNCTION fmm_branch

!-------------------------------------------------------------------------------

   FUNCTION fmm_deepest_level(scheme)

      USE fmm_qlm_builder, ONLY: fmm_system_size
      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK) :: fmm_deepest_level
      REAL(REALK) :: x

      x = fmm_system_size/scheme%grain
      ! Note we round UP the integer
      fmm_deepest_level = MAX(TOP_LEVEL,(1+INT(LOG(x)/LOG(two))))

      IF (fmm_deepest_level > MAX_LEVEL) THEN
         stat_level_saturation = fmm_deepest_level
         fmm_deepest_level = MAX_LEVEL
      END IF

   END FUNCTION fmm_deepest_level

!-------------------------------------------------------------------------------
! Return the box size at a given level, based on the smallest box size.

   FUNCTION fmm_grain(scheme,level)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER(INTK),      INTENT(IN) :: level
      REAL(REALK) :: fmm_grain
      INTEGER(INTK) :: deepest_level

      deepest_level = fmm_deepest_level(scheme)
      fmm_grain = (scheme%grain)*(2**(deepest_level-level))

   END FUNCTION fmm_grain

!-------------------------------------------------------------------------------
! Routine to define the separation parameter for NF interactions
! i.e. boxes closer than this separation are treated as near field.
! only designed for cases where boxes are at the same level (same size)

   SUBROUTINE fmm_def_WS_NF(LHS,RHS,id,WS_para)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN)  :: id
      INTEGER(INTK),      INTENT(OUT) :: WS_para

      WS_para = (LHS%box_paras(id%LHS)%bra + RHS%box_paras(id%RHS)%bra) /2

   END SUBROUTINE fmm_def_WS_NF

!-------------------------------------------------------------------------------
! Routine to define local space that contains both LFF and NF space;
! i.e. boxes further apart than this separation are RFF boxes;
! note that this routine is generalised for cases where the boxes are
! at different levels (of different sizes) as is needed for the NlogN
! contraction scheme; in this case, we must choose which box size to express
! WS_para in terms of; by default we take the box size at highest level

   SUBROUTINE fmm_def_WS_RFF(LHS_in,RHS_in,id,WS_para)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN)  :: LHS_in, RHS_in
      TYPE(LHS_RHS_type), INTENT(IN)  :: id
      INTEGER(INTK),      INTENT(OUT) :: WS_para

      TYPE(box_mm_paras) :: LHS, RHS
      LHS = LHS_in%box_paras(id%LHS)
      RHS = RHS_in%box_paras(id%RHS)

      ! we define WS_para in terms of the highest common grid
      CALL fmm_translate_to_common_grid(LHS,RHS)

      ! LFF is related to the parent's NF space
      LHS%bra = fmm_parent_bra(LHS%bra)
      RHS%bra = fmm_parent_bra(RHS%bra)

      ! hence NF-space for _parents_ of common grid
      WS_para = (LHS%bra + RHS%bra) /2
      ! hence LFF-space in common grid itself is
      WS_para = 1 + 2*WS_para
      ! +1 because of asymmetry in LFF space, and we are conservative.

   END SUBROUTINE fmm_def_WS_RFF

!-------------------------------------------------------------------------------
! Are 2 boxes within their near-field space (hence too close to be
! treated by boxed multipole expansions exactly)?
! If number of boxes between A and B < ((Br(A)+Br(B))/2) they are NF

   FUNCTION fmm_NF_boxes(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: fmm_NF_boxes
      INTEGER(INTK) :: WS_NF
      INTEGER(INTK) :: space

      ! test only valid if box grids are at the same level depth
      IF (LHS%level /= RHS%level) CALL fmm_quit('levels not equal in NF_boxes')

      WS_NF = (LHS%bra + RHS%bra)/2

      fmm_NF_boxes = .FALSE.
      space = ABS(LHS%box(3) - RHS%box(3))
      IF (space <= WS_NF) THEN
         space = ABS(LHS%box(2) - RHS%box(2))
         IF (space <= WS_NF) THEN
            space = ABS(LHS%box(1) - RHS%box(1))
            IF (space <= WS_NF) THEN
               ! pair are in near field
               fmm_NF_boxes = .TRUE.
            END IF
         END IF
      END IF

   END FUNCTION fmm_NF_boxes

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Remote Far Field"
! based on separation.  Used by FMM and NlogN algorithms.

   FUNCTION fmm_RFF_boxes(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: fmm_RFF_boxes

      TYPE(box_mm_paras) :: LHS_up, RHS_up

      ! test only valid if box grids are at the same level depth
      IF (LHS%level /= RHS%level) CALL fmm_quit('levels in fmm_RFF_boxes')

      LHS_up = LHS
      RHS_up = RHS
      LHS_up%box = fmm_parent_box(LHS%box)
      LHS_up%bra = fmm_parent_bra(LHS%bra)
      RHS_up%box = fmm_parent_box(RHS%box)
      RHS_up%bra = fmm_parent_bra(RHS%bra)

      fmm_RFF_boxes = .NOT.(fmm_NF_boxes(LHS_up,RHS_up))

   END FUNCTION fmm_RFF_boxes

!-------------------------------------------------------------------------------

   FUNCTION fmm_same_box(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: fmm_same_box

      IF (LHS%level /= RHS%level) CALL fmm_quit('levels not equal in same_box')
      fmm_same_box = .FALSE.
      IF (LHS%box(1) == RHS%box(1)  .AND.   &
          LHS%box(2) == RHS%box(2)  .AND.   &
          LHS%box(3) == RHS%box(3))         &
      THEN
         fmm_same_box = .TRUE.
      END IF

   END FUNCTION fmm_same_box

!-------------------------------------------------------------------------------
! Routine to generate a temporary set of box parameters at a common grid
! level.  Parameters at the deeper level are translated to their parent's
! grid until the levels match.

   SUBROUTINE fmm_translate_to_common_grid(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: LHS, RHS

      IF (LHS%level == RHS%level) RETURN
      IF (LHS%level > RHS%level) THEN
         DO WHILE (LHS%level > RHS%level)
            LHS%box = fmm_parent_box(LHS%box)
            LHS%bra = fmm_parent_bra(LHS%bra)
            LHS%level = LHS%level -1
         END DO
      ELSE
         DO WHILE (RHS%level > LHS%level)
            RHS%box = fmm_parent_box(RHS%box)
            RHS%bra = fmm_parent_bra(RHS%bra)
            RHS%level = RHS%level -1
         END DO
      END IF

   END SUBROUTINE fmm_translate_to_common_grid

!-------------------------------------------------------------------------------

END MODULE fmm_box_utils

