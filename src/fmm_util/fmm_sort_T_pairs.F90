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
MODULE fmm_sort_T_pairs

! Sorting module with a very simple O(N^2) insertion sort
! algorithm or using the average O(NlogN) Quicksort algorithm
!
! Sort array of interaction T-pairs using Quicksort wrt
!     (a) RHS index
!     (b) T-vector component
!     (c) T-vector modulus (ratio)

   USE fmm_global_paras, ONLY: INTK, REALK, T_pair_single
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_quicksort_wrt_RHS,           &
             fmm_quicksort_wrt_vector,        &
             fmm_quicksort_wrt_ratio

   ! criteria for switching from Quicksort to Insertion-Sort algorithms
   INTEGER(INTK), PARAMETER :: quicksort_CUTOFF = 10

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_RHS(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      TYPE(T_pair_single) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%paras%RHS_id > tmp%paras%RHS_id) THEN
               arr(j+1) = arr(j)
            ELSE
               arr(j+1) = tmp
               CYCLE iloop
            END IF
         END DO
         arr(1) = tmp
      END DO iloop

   END SUBROUTINE insertion_sort_wrt_RHS

!-------------------------------------------------------------------------------

   SUBROUTINE RHS_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE RHS_swap_elements

!-------------------------------------------------------------------------------
! Subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION RHS_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: left, right
      INTEGER(INTK) :: RHS_median_of_three
      INTEGER(INTK) :: cntr

      cntr = (left + right)/2

      IF (arr(left)%paras%RHS_id > arr(cntr)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,left,cntr)
      IF (arr(left)%paras%RHS_id > arr(right)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,left,right)
      IF (arr(cntr)%paras%RHS_id > arr(right)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,cntr,right)

      CALL RHS_swap_elements(arr,cntr,right-1_INTK)
      RHS_median_of_three = arr(right-1)%paras%RHS_id

   END FUNCTION RHS_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_RHS(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      INTEGER(INTK) :: left, right
      INTEGER(INTK) :: i,j, pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_RHS(arr)
         RETURN
      END IF

      pivot = RHS_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%paras%RHS_id < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%paras%RHS_id > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL RHS_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be
      ! swap pivot back to middle of array
      CALL RHS_swap_elements(arr,i,right-1_INTK)

      ! now sort sub-arrays either side of pivot
      CALL fmm_quicksort_wrt_RHS( arr( left:(i-1) ) )
      CALL fmm_quicksort_wrt_RHS( arr( (i+1):right ) )

   END SUBROUTINE fmm_quicksort_wrt_RHS

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt 1 component of T-vector
!-----------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_vector(arr,xyz)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: xyz

      TYPE(T_pair_single) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%r_ab(xyz) > tmp%r_ab(xyz)) THEN
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
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

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
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: left, right
      INTEGER(INTK),       INTENT(IN)    :: xyz

      REAL(REALK)   :: vector_median_of_three
      INTEGER(INTK) :: cntr

      cntr = (left + right)/2

      IF (arr(left)%r_ab(xyz) > arr(cntr)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,left,cntr)
      IF (arr(left)%r_ab(xyz) > arr(right)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,left,right)
      IF (arr(cntr)%r_ab(xyz) > arr(right)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,cntr,right)

      CALL vector_swap_elements(arr,cntr,right-1_INTK)
      vector_median_of_three = arr(right-1)%r_ab(xyz)

   END FUNCTION vector_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_vector(arr,xyz)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: xyz

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
         DO WHILE ( arr(i)%r_ab(xyz) < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%r_ab(xyz) > pivot )
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
! Insertion sort for sorting of T-pair batch wrt T-vector modulus (ratio)
!------------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_ratio(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      TYPE(T_pair_single) :: tmp
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr)
         tmp = arr(i)
         DO j = (i-1), 1, -1
            IF (arr(j)%paras%ratio > tmp%paras%ratio) THEN
               arr(j+1) = arr(j)
            ELSE
               arr(j+1) = tmp
               CYCLE iloop
            END IF
         END DO
         arr(1) = tmp
      END DO iloop

   END SUBROUTINE insertion_sort_wrt_ratio

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt T-vector modulus (ratio)
!-------------------------------------------------------------------

   SUBROUTINE ratio_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE ratio_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION ratio_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER(INTK),       INTENT(IN)    :: left, right

      REAL(REALK)   :: ratio_median_of_three
      INTEGER(INTK) :: cntr

      cntr = (left + right)/2

      IF (arr(left)%paras%ratio > arr(cntr)%paras%ratio)  &
         CALL ratio_swap_elements(arr,left,cntr)
      IF (arr(left)%paras%ratio > arr(right)%paras%ratio)  &
         CALL ratio_swap_elements(arr,left,right)
      IF (arr(cntr)%paras%ratio > arr(right)%paras%ratio)  &
         CALL ratio_swap_elements(arr,cntr,right)

      CALL ratio_swap_elements(arr,cntr,right-1_INTK)
      ratio_median_of_three = arr(right-1)%paras%ratio

   END FUNCTION ratio_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_quicksort_wrt_ratio(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      INTEGER(INTK) :: left, right
      INTEGER(INTK) :: i,j
      REAL(REALK)   :: pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_ratio(arr)
         RETURN
      END IF

      pivot = ratio_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%paras%ratio < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%paras%ratio > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL ratio_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be
      ! swap pivot back to middle of array
      CALL ratio_swap_elements(arr,i,right-1_INTK)

      ! now sort sub-arrays either side of pivot
      CALL fmm_quicksort_wrt_ratio( arr( left:(i-1) ) )
      CALL fmm_quicksort_wrt_ratio( arr( (i+1):right ) )

   END SUBROUTINE fmm_quicksort_wrt_ratio

!-------------------------------------------------------------------------------

END MODULE fmm_sort_T_pairs

