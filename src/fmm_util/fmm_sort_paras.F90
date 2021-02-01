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
MODULE fmm_sort_paras

! Sorting module with a very simple O(N^2) insertion sort
! algorithm or using the average O(NlogN) Quicksort algorithm
!
! Sort multipole moment parameters wrt
!     (a) expansion centre component
!     (b) box component
!     (c) branch

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_quicksort_wrt_branches,      &
             fmm_quicksort_wrt_vector,        &
             fmm_quicksort_wrt_boxes

   ! criteria for switching from Quicksort to Insertion-Sort algorithms
   INTEGER(INTK), PARAMETER :: quicksort_CUTOFF = 10

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_branches(arr)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      TYPE(box_mm_paras) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%bra > tmp%bra) THEN
               arr(j+1) = arr(j)
            ELSE
               arr(j+1) = tmp
               CYCLE iloop
            END IF
         END DO
         arr(1) = tmp
      END DO iloop

   END SUBROUTINE insertion_sort_wrt_branches

!-------------------------------------------------------------------------------

   SUBROUTINE bra_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: i, j
      TYPE(box_mm_paras) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE bra_swap_elements

!-------------------------------------------------------------------------------
! Subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION bra_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: left, right
      INTEGER(INTK) :: bra_median_of_three
      INTEGER(INTK) :: cntr

      cntr = (left + right)/2
      IF(arr(left)%bra > arr(cntr)%bra)  CALL bra_swap_elements(arr,left,cntr)
      IF(arr(left)%bra > arr(right)%bra) CALL bra_swap_elements(arr,left,right)
      IF(arr(cntr)%bra > arr(right)%bra) CALL bra_swap_elements(arr,cntr,right)
      CALL bra_swap_elements(arr,cntr,right-1_INTK)
      bra_median_of_three = arr(right-1)%bra

   END FUNCTION bra_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_branches(arr)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)

      INTEGER(INTK) :: left, right
      INTEGER(INTK) :: i,j, pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_branches(arr)
         RETURN
      END IF

      pivot = bra_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%bra < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%bra > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL bra_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be
      ! swap pivot back to middle of array
      CALL bra_swap_elements(arr,i,right-1_INTK)

      ! now sort sub-arrays either side of pivot
      CALL fmm_quicksort_wrt_branches( arr( left:(i-1) ) )
      CALL fmm_quicksort_wrt_branches( arr( (i+1):right ) )

   END SUBROUTINE fmm_quicksort_wrt_branches

!-------------------------------------------------------------------------------
! Insertion sort for sorting of multipole moment centres
!-----------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_vector(arr,xyz)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: xyz

      TYPE(raw_mm_paras) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%cntr(xyz) > tmp%cntr(xyz)) THEN
               arr(j+1) = arr(j)
            ELSE
               arr(j+1) = tmp
               CYCLE iloop
            END IF
         END DO
         arr(1) = tmp
      END DO iloop

   END SUBROUTINE insertion_sort_wrt_vector

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt 1 component of T-vector
!------------------------------------------------------------------

   SUBROUTINE vector_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: i, j
      TYPE(raw_mm_paras) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE vector_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION vector_median_of_three(arr,left,right,xyz)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: left, right
      INTEGER(INTK),      INTENT(IN)    :: xyz

      REAL(REALK)   :: vector_median_of_three
      INTEGER(INTK) :: cntr

      cntr = (left + right)/2

      IF (arr(left)%cntr(xyz) > arr(cntr)%cntr(xyz))  &
         CALL vector_swap_elements(arr,left,cntr)
      IF (arr(left)%cntr(xyz) > arr(right)%cntr(xyz))  &
         CALL vector_swap_elements(arr,left,right)
      IF (arr(cntr)%cntr(xyz) > arr(right)%cntr(xyz))  &
         CALL vector_swap_elements(arr,cntr,right)

      CALL vector_swap_elements(arr,cntr,right-1_INTK)
      vector_median_of_three = arr(right-1)%cntr(xyz)

   END FUNCTION vector_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_vector(arr,xyz)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: xyz

      INTEGER(INTK) :: left, right
      INTEGER(INTK) :: i,j
      REAL(REALK)   :: pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_vector(arr,xyz)
         RETURN
      END IF

      pivot = vector_median_of_three(arr,left,right,xyz)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%cntr(xyz) < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%cntr(xyz) > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL vector_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be
      ! swap pivot back to middle of array
      CALL vector_swap_elements(arr,i,right-1_INTK)

      ! now sort sub-arrays either side of pivot
      CALL fmm_quicksort_wrt_vector( arr( left:(i-1) ),xyz )
      CALL fmm_quicksort_wrt_vector( arr( (i+1):right ),xyz )

   END SUBROUTINE fmm_quicksort_wrt_vector

!-------------------------------------------------------------------------------
! Insertion sort for sorting wrt one box index
!-----------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_boxes(arr,xyz)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: xyz

      TYPE(box_mm_paras) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%box(xyz) > tmp%box(xyz)) THEN
               arr(j+1) = arr(j)
            ELSE
               arr(j+1) = tmp
               CYCLE iloop
            END IF
         END DO
         arr(1) = tmp
      END DO iloop

   END SUBROUTINE insertion_sort_wrt_boxes

!-------------------------------------------------------------------------------
! Quicksort for sorting wrt one box index
!------------------------------------------------------------------

   SUBROUTINE boxes_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: i, j
      TYPE(box_mm_paras) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE boxes_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION boxes_median_of_three(arr,left,right,xyz)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: left, right
      INTEGER(INTK),      INTENT(IN)    :: xyz

      INTEGER(INTK) :: boxes_median_of_three
      INTEGER(INTK) :: box

      box = (left + right)/2

      IF (arr(left)%box(xyz) > arr(box)%box(xyz))  &
         CALL boxes_swap_elements(arr,left,box)
      IF (arr(left)%box(xyz) > arr(right)%box(xyz))  &
         CALL boxes_swap_elements(arr,left,right)
      IF (arr(box)%box(xyz) > arr(right)%box(xyz))  &
         CALL boxes_swap_elements(arr,box,right)

      CALL boxes_swap_elements(arr,box,right-1_INTK)
      boxes_median_of_three = arr(right-1)%box(xyz)

   END FUNCTION boxes_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_boxes(arr,xyz)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),      INTENT(IN)    :: xyz

      INTEGER(INTK) :: left, right
      INTEGER(INTK) :: i,j
      INTEGER(INTK) :: pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_boxes(arr,xyz)
         RETURN
      END IF

      pivot = boxes_median_of_three(arr,left,right,xyz)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%box(xyz) < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%box(xyz) > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL boxes_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be
      ! swap pivot back to middle of array
      CALL boxes_swap_elements(arr,i,right-1_INTK)

      ! now sort sub-arrays either side of pivot
      CALL fmm_quicksort_wrt_boxes( arr( left:(i-1) ),xyz )
      CALL fmm_quicksort_wrt_boxes( arr( (i+1):right ),xyz )

   END SUBROUTINE fmm_quicksort_wrt_boxes

!-------------------------------------------------------------------------------

END MODULE fmm_sort_paras

!==============================================================================
!
!MODULE fmm_sort_sh_pairs
!
!! Sorting module with a very simple O(N^2) insertion sort
!! algorithm or using the average O(NlogN) Quicksort algorithm
!!
!! Sort shell pair array wrt
!!     (a) box component
!
!   USE fmm_global_paras
!   IMPLICIT NONE
!   PRIVATE
!   ! Public procedures
!   PUBLIC :: fmm_quicksort_wrt_boxes
!
!   ! criteria for switching from Quicksort to Insertion-Sort algorithms
!   INTEGER(INTK), PARAMETER :: quicksort_CUTOFF = 10
!
!CONTAINS
!
!!-------------------------------------------------------------------------------
!! Insertion sort for sorting wrt one box index
!!-----------------------------------------------------------------------
!
!   SUBROUTINE insertion_sort_wrt_boxes(arr,xyz)
!
!      IMPLICIT NONE
!      TYPE(fmm_shell_pair_node), INTENT(INOUT) :: arr(:)
!      INTEGER(INTK),             INTENT(IN)    :: xyz
!
!      TYPE(fmm_shell_pair_node) :: tmp
!      INTEGER(INTK) :: i, j
!
!      iloop: DO i = 1, SIZE(arr)
!         tmp = arr(i)
!         DO j = (i-1), 1, -1
!            IF (arr(j)%box(xyz) > tmp%box(xyz)) THEN
!               arr(j+1) = arr(j)
!            ELSE
!               arr(j+1) = tmp
!               CYCLE iloop
!            END IF
!         END DO
!         arr(1) = tmp
!      END DO iloop
!
!   END SUBROUTINE insertion_sort_wrt_boxes
!
!!-------------------------------------------------------------------------------
!! Quicksort for sorting wrt one box index
!!------------------------------------------------------------------
!
!   SUBROUTINE boxes_swap_elements(arr,i,j)
!
!      IMPLICIT NONE
!      TYPE(fmm_shell_pair_node), INTENT(INOUT) :: arr(:)
!      INTEGER(INTK),             INTENT(IN)    :: i, j
!      TYPE(fmm_shell_pair_node) :: tmp
!
!      tmp = arr(i)
!      arr(i) = arr(j)
!      arr(j) = tmp
!
!   END SUBROUTINE boxes_swap_elements
!
!!-------------------------------------------------------------------------------
!! subroutine to return the median of first, last, and centre array elements.
!! order these and hide the pivot (median) as the second to last array element
!! (acts as a sentinel whilst the last is the largest and put in correct place.)
!
!   FUNCTION boxes_median_of_three(arr,left,right,xyz)
!
!      IMPLICIT NONE
!      TYPE(fmm_shell_pair_node), INTENT(INOUT) :: arr(:)
!      INTEGER(INTK),             INTENT(IN)    :: left, right
!      INTEGER(INTK),             INTENT(IN)    :: xyz
!
!      INTEGER(INTK) :: boxes_median_of_three
!      INTEGER(INTK) :: box
!
!      box = (left + right)/2
!
!      IF (arr(left)%box(xyz) > arr(box)%box(xyz))  &
!         CALL boxes_swap_elements(arr,left,box)
!      IF (arr(left)%box(xyz) > arr(right)%box(xyz))  &
!         CALL boxes_swap_elements(arr,left,right)
!      IF (arr(box)%box(xyz) > arr(right)%box(xyz))  &
!         CALL boxes_swap_elements(arr,box,right)
!
!      CALL boxes_swap_elements(arr,box,right-1_INTK)
!      boxes_median_of_three = arr(right-1)%box(xyz)
!
!   END FUNCTION boxes_median_of_three
!
!!-------------------------------------------------------------------------------
!
!   RECURSIVE SUBROUTINE fmm_quicksort_wrt_boxes(arr,xyz)
!
!      IMPLICIT NONE
!      TYPE(fmm_shell_pair_node), INTENT(INOUT) :: arr(:)
!      INTEGER(INTK),             INTENT(IN)    :: xyz
!
!      INTEGER(INTK) :: left, right
!      INTEGER(INTK) :: i,j
!      INTEGER(INTK) :: pivot
!
!      left = 1
!      right = SIZE(arr)
!
!      IF (left + quicksort_CUTOFF > right) THEN
!         ! array too small for Quicksort to be efficient, so don't use it
!         CALL insertion_sort_wrt_boxes(arr,xyz)
!         RETURN
!      END IF
!
!      pivot = boxes_median_of_three(arr,left,right,xyz)
!      ! partition array elements about chosen pivot
!      i = left
!      j = right-2
!      DO ! until left, right pointers cross
!         ! increment left until large element found
!         DO WHILE ( arr(i)%box(xyz) < pivot )
!            i = i+1
!         END DO
!         ! increment right until small element found
!         DO WHILE ( arr(j)%box(xyz) > pivot )
!            j = j-1
!         END DO
!         IF (i < j) THEN
!            ! swap items (i,j)
!            CALL boxes_swap_elements(arr,i,j)
!            ! now increment (i,j) to avoid infinite loop if both = pivot
!            i = i+1
!            j = j-1
!         ELSE
!            EXIT
!         END IF
!      END DO
!      ! i is now the position where the pivot should be
!      ! swap pivot back to middle of array
!      CALL boxes_swap_elements(arr,i,right-1_INTK)
!
!      ! now sort sub-arrays either side of pivot
!      CALL fmm_quicksort_wrt_boxes( arr( left:(i-1) ),xyz )
!      CALL fmm_quicksort_wrt_boxes( arr( (i+1):right ),xyz )
!
!   END SUBROUTINE fmm_quicksort_wrt_boxes
!
!!-------------------------------------------------------------------------------
!
!END MODULE fmm_sort_sh_pairs
!
