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
MODULE fmm_multi_T_buffer

   USE fmm_stats
   USE fmm_global_paras

   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: fmm_init_multi_T_buffer,    &
             fmm_free_multi_T_buffer,    &
             fmm_multi_T_buffer_add

   INTEGER(INTK), PARAMETER :: BUFFER_SIZE = 1000
   ! module wide variables
   INTEGER(INTK),      SAVE :: ndim_max
   TYPE(T_pair_batch), SAVE :: T_pair_buffer

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_multi_T_buffer(ndim_max_in)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: ndim_max_in

!      ndim_max = ndim_max_in*2  ! we multiply by two for 'paired" algorithm
      ndim_max = ndim_max_in
      IF (ndim_max < 1) CALL fmm_quit('invalid multiple T-matrix dimension!')
      NULLIFY (T_pair_buffer%items)
      ALLOCATE (T_pair_buffer%items(BUFFER_SIZE))
      T_pair_buffer%ndim = 0

   END SUBROUTINE fmm_init_multi_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_multi_T_buffer(T_contractor)

      IMPLICIT NONE
      EXTERNAL T_contractor

      IF (.NOT.ASSOCIATED(T_pair_buffer%items))  &
         CALL fmm_quit('T_pair_buffer not alloc.')
      IF ( T_pair_buffer%ndim /= 0 ) THEN
         CALL expunge_multi_buffer(T_contractor)
         T_pair_buffer%ndim = 0
      END IF
      DEALLOCATE (T_pair_buffer%items)
      NULLIFY (T_pair_buffer%items)

   END SUBROUTINE fmm_free_multi_T_buffer

!-------------------------------------------------------------------------------
!
!   SUBROUTINE fmm_multi_T_buffer_add(T_contractor,T_pair)
!
!      IMPLICIT NONE
!      TYPE(T_pair_single), INTENT(IN) :: T_pair
!      EXTERNAL T_contractor
!
!      INTEGER(INTK), SAVE :: iRHS_last = 0
!      INTEGER(INTK) :: iRHS
!
!      iRHS = T_pair%paras%RHS_id
!
!!!!!      IF ( BTEST(T_pair_buffer%ndim+1,0) ) THEN
!         ! number of buffer entries is even; try to expunge
!         IF ( (T_pair_buffer%ndim == ndim_max)   &
!              .OR. ((iRHS /= iRHS_last) .AND. (iRHS_last /= 0)) ) THEN
!            ! expunge buffer and build all the T-matrices at once
!            CALL T_contractor(T_pair_buffer)
!            T_pair_buffer%ndim = 0
!         END IF
!         iRHS_last = T_pair%paras%RHS_id
!!!!!      END IF
!
!      T_pair_buffer%ndim = T_pair_buffer%ndim +1
!      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair
!
!   END SUBROUTINE fmm_multi_T_buffer_add
!
!-------------------------------------------------------------------------------

   SUBROUTINE fmm_multi_T_buffer_add(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor

      IF ( T_pair_buffer%ndim == BUFFER_SIZE ) THEN
         ! expunge buffer and build all the T-matrices at once
         CALL expunge_multi_buffer(T_contractor)
      END IF

      stat_tpack_total = stat_tpack_total + one
      T_pair_buffer%ndim = T_pair_buffer%ndim +1
      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

   END SUBROUTINE fmm_multi_T_buffer_add

!-------------------------------------------------------------------------------

   SUBROUTINE expunge_multi_buffer(T_contractor)

      USE fmm_sort_T_pairs, ONLY: fmm_quicksort_wrt_RHS

      IMPLICIT NONE
      EXTERNAL T_contractor

      INTEGER(INTK) :: i, lo
      INTEGER(INTK) :: iRHS, iRHS_next, item_max
      TYPE(T_pair_batch) :: ptr

      lo = 1
      item_max = MIN((BUFFER_SIZE-1),(T_pair_buffer%ndim-1))

      ! sort only if needed
      iRHS = T_pair_buffer%items(1)%paras%RHS_id
      DO i = 2, item_max
         iRHS_next = T_pair_buffer%items(i)%paras%RHS_id
         IF ( iRHS_next < iRHS ) THEN
            CALL fmm_quicksort_wrt_RHS(T_pair_buffer%items(1:item_max))
            EXIT
         END IF
         iRHS = iRHS_next
      END DO

      DO i = 1, item_max
         iRHS = T_pair_buffer%items(i)%paras%RHS_id
         iRHS_next = T_pair_buffer%items(i+1)%paras%RHS_id
         ptr%ndim = i-lo+1
         IF ((iRHS /= iRHS_next) .OR. (ptr%ndim == ndim_max)) THEN
            ptr%items => T_pair_buffer%items(lo:i)
            CALL T_contractor(ptr)
            lo = i+1
         END IF
      END DO

      ptr%ndim = (item_max+1)-lo+1
      ptr%items => T_pair_buffer%items(lo:(item_max+1))
      CALL T_contractor(ptr)

      T_pair_buffer%ndim = 0
      stat_tpack_chunks = stat_tpack_chunks + one

   END SUBROUTINE expunge_multi_buffer

!-------------------------------------------------------------------------------

END MODULE fmm_multi_T_buffer
