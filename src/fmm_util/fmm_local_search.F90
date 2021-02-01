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
MODULE fmm_local_search

   USE fmm_global_paras, ONLY: INTK, id_list, id_node, scheme_paras, box_mm_data, box_mm_paras, gen_mm_paras, LHS_raw_RHS_raw, &
                               LHS_box_RHS_box, TOP_LEVEL, DO_FMM
   USE fmm_box_utils, ONLY: fmm_parent_box, fmm_same_box, fmm_RFF_boxes

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_local_search,   &
             fmm_free_local_search,   &
             fmm_get_local_paras

   TYPE fmm_local_map
      TYPE(id_list), POINTER :: box_list(:)
   END TYPE fmm_local_map

   TYPE(fmm_local_map), ALLOCATABLE, SAVE :: map_at_level(:)
   INTEGER(INTK), SAVE :: deepest_level

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_local_search(scheme)

      USE fmm_box_builder, ONLY: fmm_get_box_paras_at_level
      USE fmm_box_utils,   ONLY: fmm_deepest_level

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme

      TYPE(fmm_local_map), ALLOCATABLE :: RHS_map(:)
      TYPE(box_mm_data),   ALLOCATABLE :: mm_paras(:)
      TYPE(box_mm_paras),  POINTER :: ptr(:)
      INTEGER(INTK) :: level, j

      IF ( .NOT. scheme%branch_free )  &
         CALL fmm_quit('local search module may not work with CFMM!')
      IF ( .NOT. scheme%algorithm == DO_FMM )  &
         CALL fmm_quit('local search module is only for FMM algorithm')

      deepest_level = fmm_deepest_level(scheme)
      NULLIFY (ptr)
      ALLOCATE (mm_paras(deepest_level))
      ALLOCATE (RHS_map(deepest_level))

      ! Get multipole moment box information over whole hierarchy
      DO level = TOP_LEVEL, deepest_level
         NULLIFY (RHS_map(level)%box_list)
         NULLIFY (mm_paras(level)%LHS_paras)
         NULLIFY (mm_paras(level)%RHS_paras)
         NULLIFY (mm_paras(level)%qlm_T)
         NULLIFY (mm_paras(level)%qlm_W)
         CALL fmm_get_box_paras_at_level(level,scheme,ptr,'LHS')
         mm_paras(level)%LHS_paras => ptr
         CALL fmm_get_box_paras_at_level(level,scheme,ptr,'RHS')
         mm_paras(level)%RHS_paras => ptr
      END DO

      ! Build parent:child map of RHS parameters and use it to make
      ! local interaction map of LHS-RHS pairs
      CALL fmm_init_local_map(mm_paras)
      CALL fmm_build_RHS_map(mm_paras,RHS_map)
      CALL fmm_build_local_map(mm_paras,RHS_map)

      ! Deallocations
      DO level = TOP_LEVEL, deepest_level
         IF (ASSOCIATED(RHS_map(level)%box_list)) THEN
            DO j = 1, SIZE(RHS_map(level)%box_list)
               CALL free_linked_list(RHS_map(level)%box_list(j)%head)
            END DO
            DEALLOCATE (RHS_map(level)%box_list)
         END IF
         NULLIFY (RHS_map(level)%box_list)
         NULLIFY (mm_paras(level)%LHS_paras)
         NULLIFY (mm_paras(level)%RHS_paras)
      END DO
      DEALLOCATE (mm_paras, RHS_map)

   END SUBROUTINE fmm_init_local_search

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_local_map(mm_paras)

      IMPLICIT NONE
      TYPE(box_mm_data), INTENT(IN) :: mm_paras(:)

      INTEGER(INTK) :: i,j, ndim

      ALLOCATE (map_at_level(deepest_level))

      DO i = TOP_LEVEL, deepest_level
         ndim = SIZE(mm_paras(i)%LHS_paras)
         NULLIFY (map_at_level(i)%box_list)
         ALLOCATE (map_at_level(i)%box_list(ndim))
         DO j = 1, ndim
            NULLIFY (map_at_level(i)%box_list(j)%head)
            map_at_level(i)%box_list(j)%occ = 0
         END DO
      END DO

   END SUBROUTINE fmm_init_local_map

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_local_search

      IMPLICIT NONE
      INTEGER(INTK) :: i,j,ndim

      IF ( deepest_level == TOP_LEVEL ) RETURN

      DO i = TOP_LEVEL, deepest_level
         IF (ASSOCIATED(map_at_level(i)%box_list)) THEN
            ndim = SIZE(map_at_level(i)%box_list)
            DO j = 1, ndim
               CALL free_linked_list(map_at_level(i)%box_list(j)%head)
            END DO
            DEALLOCATE (map_at_level(i)%box_list)
         END IF
         NULLIFY (map_at_level(i)%box_list)
      END DO

      IF (ALLOCATED(map_at_level)) DEALLOCATE (map_at_level)

   END SUBROUTINE fmm_free_local_search

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE free_linked_list(occ_node)

      TYPE(id_node), POINTER :: occ_node

      IF (.NOT. ASSOCIATED(occ_node)) RETURN

      ! traverse linked-list freeing from bottom up recursively
      IF (ASSOCIATED(occ_node%next)) THEN
         CALL free_linked_list(occ_node%next)
         IF (ASSOCIATED(occ_node)) DEALLOCATE(occ_node)
         NULLIFY(occ_node)
      END IF
      IF (ASSOCIATED(occ_node)) DEALLOCATE(occ_node)
      NULLIFY(occ_node)

   END SUBROUTINE free_linked_list

!-------------------------------------------------------------------------------
! For each box, a list of the children boxes is generated.
! Box occupancy is only computed for levels _above_ the deepest
! level since deepest level occupancy corresponds to raw moments,
! not boxed moments, which is done separately.

   SUBROUTINE fmm_build_RHS_map(mm_paras,RHS_map)

      IMPLICIT NONE
      TYPE(box_mm_data),   INTENT(IN)  :: mm_paras(:)
      TYPE(fmm_local_map), INTENT(OUT) :: RHS_map(:)

      TYPE(box_mm_paras) :: box
      INTEGER(INTK) :: i,id, level, ndim

      DO level = deepest_level-1, TOP_LEVEL, -1

         ndim = SIZE(mm_paras(level)%RHS_paras)
         ALLOCATE (RHS_map(level)%box_list(ndim))

         ! Initialise
         DO i = 1, ndim
            NULLIFY (RHS_map(level)%box_list(i)%head)
            RHS_map(level)%box_list(i)%occ = 0
         END DO

         ! We drive the loop over the child box level
         DO i = 1, SIZE(mm_paras(level+1)%RHS_paras)
            ! 'map_up' indexes moments but, if in sync, then also the para set
            id = mm_paras(level+1)%RHS_paras(i)%map_up
            box%box = fmm_parent_box(mm_paras(level+1)%RHS_paras(i)%box)
            box%level = level
            ! 'id' should point to RHS parent box parameters
            IF ( .NOT. fmm_same_box(mm_paras(level)%RHS_paras(id),box))  &
               CALL fmm_quit('RHS paras map-up to parent broken')
            ! Append item to linked list
            CALL fmm_add_item(RHS_map(level)%box_list(id),i)
         END DO

      END DO

   END SUBROUTINE fmm_build_RHS_map

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_add_item(box_map,raw_id)

      IMPLICIT NONE
      TYPE(id_list), INTENT(INOUT) :: box_map
      INTEGER(INTK), INTENT(IN) :: raw_id
      TYPE(id_node), POINTER :: new_node

      IF ( box_map%occ == 0 ) THEN
         ! Start new linked list for these boxes
         box_map%occ = 1
         ALLOCATE (box_map%head)
         box_map%head%id = raw_id
         NULLIFY (box_map%head%next)
         RETURN
      END IF

      box_map%occ = box_map%occ + 1

      ALLOCATE(new_node)
      new_node%id = raw_id
      IF (ASSOCIATED(box_map%head%next)) THEN
         ! More than one entry in list (including head)
         ! so point new_node to old second entry
         new_node%next => box_map%head%next
         ! Point head to new_node
         NULLIFY(box_map%head%next)
         box_map%head%next => new_node
      ELSE
         ! Only head so far; make new_node our second entry
         box_map%head%next => new_node
         NULLIFY(new_node%next)   ! end of list
      END IF

      END SUBROUTINE fmm_add_item

!-------------------------------------------------------------------------------
! Here we generate a list for each LHS box of all
! the RHS boxes spatially local to the LHS box.
! For now 'local' includes all the parent's near-field boxes.

   SUBROUTINE fmm_build_local_map(mm_paras,RHS_map)

      IMPLICIT NONE
      TYPE(box_mm_data),   INTENT(IN) :: mm_paras(:)
      TYPE(fmm_local_map), INTENT(IN) :: RHS_map(:)

      TYPE(box_mm_paras) :: LHS_box, RHS_box
      INTEGER(INTK) :: i,j,id, lev
      TYPE(id_node), POINTER :: LHS_ptr, RHS_ptr

      lev = TOP_LEVEL
      DO i = 1, SIZE(mm_paras(lev)%LHS_paras)
         DO j = 1, SIZE(mm_paras(lev)%RHS_paras)
            LHS_box = mm_paras(lev)%LHS_paras(i)
            RHS_box = mm_paras(lev)%RHS_paras(j)
            IF ( .NOT. fmm_RFF_boxes(LHS_box,RHS_box) )  THEN
               CALL fmm_add_item(map_at_level(lev)%box_list(i),j)
            END IF
         END DO
      END DO

      DO lev = TOP_LEVEL+1, deepest_level
         DO i = 1, SIZE(mm_paras(lev)%LHS_paras)

            ! 'map_up' indexes moments but, if in sync, then also the para set
            id = mm_paras(lev)%LHS_paras(i)%map_up
            LHS_box%box = fmm_parent_box(mm_paras(lev)%LHS_paras(i)%box)
            LHS_box%level = lev-1
            ! 'id' should point to LHS parent box parameters
            IF ( .NOT. fmm_same_box(mm_paras(lev-1)%LHS_paras(id),LHS_box)) &
               CALL fmm_quit('LHS paras map-up to parent broken')

            IF ( map_at_level(lev-1)%box_list(id)%occ == 0 ) CYCLE
            LHS_ptr => map_at_level(lev-1)%box_list(id)%head
            ! Loop over parent's near-field boxes
            parent: DO
               IF (RHS_map(lev-1)%box_list(LHS_ptr%id)%occ == 0) THEN
                  IF (.NOT.ASSOCIATED(LHS_ptr%next)) EXIT parent
                  LHS_ptr => LHS_ptr%next
                  CYCLE parent
               END IF
               RHS_ptr => RHS_map(lev-1)%box_list(LHS_ptr%id)%head
               ! Get children and add to list of local boxes
               LHS_box = mm_paras(lev)%LHS_paras(i)
               child: DO
                  RHS_box = mm_paras(lev)%RHS_paras(RHS_ptr%id)
                  IF ( .NOT. fmm_RFF_boxes(LHS_box,RHS_box) )  THEN
                    CALL fmm_add_item(map_at_level(lev)%box_list(i),RHS_ptr%id)
                  END IF
                  IF (.NOT.ASSOCIATED(RHS_ptr%next)) EXIT child
                  RHS_ptr => RHS_ptr%next
               END DO child

               IF (.NOT.ASSOCIATED(LHS_ptr%next)) EXIT parent
               LHS_ptr => LHS_ptr%next
            END DO parent

         END DO
      END DO

   END SUBROUTINE fmm_build_local_map

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_local_paras(id,RHS_all,pair_type,RHS_local,ndim)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(OUT) :: RHS_local
      TYPE(gen_mm_paras), INTENT(IN)  :: RHS_all
      INTEGER(INTK),      INTENT(OUT) :: ndim
      INTEGER(INTK),      INTENT(IN)  :: pair_type
      INTEGER(INTK),      INTENT(IN)  :: id

      TYPE(id_node), POINTER :: ptr
      INTEGER(INTK) :: level, i

      SELECT CASE (pair_type)
         CASE (LHS_raw_RHS_raw)
            CALL fmm_quit ('local_paras: raw_raw NYI')

         CASE (LHS_box_RHS_box)

            IF (ASSOCIATED(RHS_local%box_paras)) CALL fmm_quit('RHS_local')
            level = RHS_all%box_paras(1)%level
            ndim = map_at_level(level)%box_list(id)%occ
            IF (ndim == 0) RETURN
            ALLOCATE (RHS_local%box_paras(ndim))

            ptr => map_at_level(level)%box_list(id)%head
            i = 0
            DO ! over linked-list of local boxes
               i = i+1
               RHS_local%box_paras(i) = RHS_all%box_paras(ptr%id)
               IF (.NOT.ASSOCIATED(ptr%next)) EXIT
               ptr => ptr%next
            END DO

         CASE DEFAULT
            CALL fmm_quit ('local_paras: requested T_pair type!')
      END SELECT

   END SUBROUTINE fmm_get_local_paras

!-------------------------------------------------------------------------------

END MODULE fmm_local_search
