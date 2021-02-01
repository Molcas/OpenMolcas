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
MODULE fmm_scale_T_buffer

   USE fmm_global_paras, ONLY: INTK, REALK, T_pair_batch, T_pair_single, Zero, One
   USE fmm_stats, ONLY: stat_tpack_chunks, stat_tpack_unique, stat_tpack_total

   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: fmm_init_scale_T_buffer,    &
             fmm_free_scale_T_buffer,    &
             fmm_scale_T_buffer_add

   INTEGER(INTK), PARAMETER :: BUFFER_SIZE = 500000
   ! module wide variables
   TYPE(T_pair_batch), SAVE :: T_pair_buffer

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_scale_T_buffer

      IMPLICIT NONE

      NULLIFY (T_pair_buffer%items)
      ALLOCATE (T_pair_buffer%items(BUFFER_SIZE))
      T_pair_buffer%ndim = 0

   END SUBROUTINE fmm_init_scale_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_scale_T_buffer(T_contractor)

      IMPLICIT NONE
      EXTERNAL T_contractor

      IF (.NOT.ASSOCIATED(T_pair_buffer%items))  &
         CALL fmm_quit('T_pair_buffer not alloc.')
      IF ( T_pair_buffer%ndim /= 0 ) THEN
         CALL expunge_scale_buffer(T_contractor)
         T_pair_buffer%ndim = 0
      END IF
      DEALLOCATE (T_pair_buffer%items)
      NULLIFY (T_pair_buffer%items)

   END SUBROUTINE fmm_free_scale_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_scale_T_buffer_add(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor
      REAL(REALK) :: ratio

      stat_tpack_total = stat_tpack_total + one
      T_pair_buffer%ndim = T_pair_buffer%ndim +1
      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

      ! normalise T-vectors held in buffer for sorting purposes
      ratio = SQRT(SUM(T_pair%r_ab*T_pair%r_ab))
      ! we sort wrt x-axis first, and we want opposite vectors to be together
      IF (T_pair%r_ab(1) < zero) ratio = -ratio
      T_pair_buffer%items(T_pair_buffer%ndim)%paras%ratio = ratio
      T_pair_buffer%items(T_pair_buffer%ndim)%r_ab = T_pair%r_ab/ratio

      IF ( T_pair_buffer%ndim == BUFFER_SIZE ) THEN
         ! sort the buffer and pass all T-pairs to contractor
         CALL expunge_scale_buffer(T_contractor)
      END IF

   END SUBROUTINE fmm_scale_T_buffer_add

!-------------------------------------------------------------------------------

   SUBROUTINE expunge_scale_buffer(T_contractor)

      USE fmm_sort_T_pairs, ONLY: fmm_quicksort_wrt_vector,     &
                                  fmm_quicksort_wrt_ratio

      IMPLICIT NONE
      EXTERNAL T_contractor

      TYPE(T_pair_batch) :: ptr, ptr2
      INTEGER(INTK) :: i, lo, hi
      REAL(REALK)   :: q1,q2

      ptr%ndim = MIN(BUFFER_SIZE,T_pair_buffer%ndim)
      ptr%items => T_pair_buffer%items(1:ptr%ndim)

      ! recursively sort wrt T-vectors, starting with x-component
      CALL sort_wrt_axis(1_INTK,ptr%items)

      ! expunge in batches of the same T-vector direction
      lo = 1
      DO i = 2, ptr%ndim
         q1 = ptr%items(i)%r_ab(1)
         q2 = ptr%items(i-1)%r_ab(1)
         IF (q1 == q2) THEN
            q1 = ptr%items(i)%r_ab(2)
            q2 = ptr%items(i-1)%r_ab(2)
            IF (q1 == q2) THEN
               q1 = ptr%items(i)%r_ab(3)
               q2 = ptr%items(i-1)%r_ab(3)
               IF (q1 == q2) CYCLE
            END IF
         END IF
         hi = i-1
         ptr2%ndim = hi-lo+1
         ptr2%items => ptr%items(lo:hi)
         stat_tpack_unique = stat_tpack_unique + one
         CALL T_contractor(ptr2)
         lo = i
      END DO

      ! finally do last batch
      hi = ptr%ndim
      ptr2%ndim = hi-lo+1
      ptr2%items => ptr%items(lo:hi)
      stat_tpack_unique = stat_tpack_unique + one
      CALL T_contractor(ptr2)

      T_pair_buffer%ndim = 0
      stat_tpack_chunks = stat_tpack_chunks + one

   CONTAINS

!-------------------------------------------------------------------------------

      RECURSIVE SUBROUTINE sort_wrt_axis(xyz,items)

         IMPLICIT NONE
         INTEGER(INTK),       INTENT(IN)    :: xyz
         TYPE(T_pair_single), INTENT(INOUT) :: items(:)

         INTEGER(INTK) :: i, lo, hi
         REAL(REALK)   :: q1,q2

         IF (SIZE(items) == 1) RETURN

         ! sort only if needed
         q1 = items(1)%r_ab(xyz)
         DO i = 2, SIZE(items)
            q2 = items(i)%r_ab(xyz)
            IF ( q2 < q1 ) THEN
               CALL fmm_quicksort_wrt_vector(items,xyz)
               EXIT
            END IF
            q1 = q2
         END DO

         ! sub-sort next T-vector component
         lo = 1
         DO i = 2, SIZE(items)
            q1 = items(i-1)%r_ab(xyz)
            q2 = items(i)%r_ab(xyz)
            IF ( q2 /= q1 ) THEN
               hi = i-1
               IF (xyz == 3) THEN
                  CALL fmm_quicksort_wrt_ratio(items(lo:hi))
                  RETURN
               ELSE
                  CALL sort_wrt_axis(xyz+1_INTK,items(lo:hi))
               END IF
               lo = i
            END IF
         END DO

         ! do last batch
         hi = SIZE(items)
         IF (xyz == 3) THEN
            CALL fmm_quicksort_wrt_ratio(items(lo:hi))
            RETURN
         ELSE
            CALL sort_wrt_axis(xyz+1_INTK,items(lo:hi))
         END IF

      END SUBROUTINE sort_wrt_axis

!-------------------------------------------------------------------------------

   END SUBROUTINE expunge_scale_buffer

!-------------------------------------------------------------------------------

END MODULE fmm_scale_T_buffer
