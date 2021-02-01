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
MODULE fmm_2darray_sort

! Sorting module with a very simple O(N^2) insertion sort algorithm
! for integer 2-d arrays like (3,:) according to input key 1-3
!
   USE fmm_global_paras, ONLY: INTK
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_insertion_sort

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_insertion_sort(arr,xyz)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(INOUT) :: arr(:,:)
      INTEGER(INTK), INTENT(IN)    :: xyz

      INTEGER(INTK) :: tmp(3)
      INTEGER(INTK) :: i, j

      iloop: DO i = 1, SIZE(arr,2)
         tmp(:) = arr(:,i)
         DO j = (i-1), 1, -1
            IF (arr(xyz,j) > tmp(xyz)) THEN
               arr(:,j+1) = arr(:,j)
            ELSE
               arr(:,j+1) = tmp(:)
               CYCLE iloop
            END IF
         END DO
         arr(:,1) = tmp(:)
      END DO iloop

   END SUBROUTINE fmm_insertion_sort

!-------------------------------------------------------------------------------

END MODULE fmm_2darray_sort
