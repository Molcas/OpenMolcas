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
MODULE fmm_boundary

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, raw_mm_paras, fmm_planes, scheme_paras, Zero, One

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_opt_near_field

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE verify_planes(grid,planes)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: grid(:)
      TYPE(fmm_planes),   INTENT(INOUT) :: planes

      REAL(REALK), PARAMETER :: THR = 1.0D-15
      INTEGER(INTK) :: i

      DO i = 1, SIZE(grid)
         IF ( ABS(grid(i)%cntr(1) - planes%xmin) > THR .AND.     &
              ABS(grid(i)%cntr(1) - planes%xmax) > THR .AND.     &
              ABS(grid(i)%cntr(2) - planes%ymin) > THR .AND.     &
              ABS(grid(i)%cntr(2) - planes%ymax) > THR .AND.     &
              ABS(grid(i)%cntr(3) - planes%zmin) > THR .AND.     &
              ABS(grid(i)%cntr(3) - planes%zmax) > THR ) THEN
            CALL fmm_quit('boundary planes not // to coordinate axes!')
         END IF
      END DO

   END SUBROUTINE verify_planes

!-------------------------------------------------------------------------------

   SUBROUTINE get_boundary_planes(grid,planes)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: grid(:)
      TYPE(fmm_planes),   INTENT(INOUT) :: planes

      INTEGER(INTK) :: i

      planes%xmin = zero
      planes%xmax = zero
      planes%ymin = zero
      planes%ymax = zero
      planes%zmin = zero
      planes%zmax = zero
      DO i = 1, SIZE(grid)
         planes%xmin = MIN(planes%xmin,grid(i)%cntr(1))
         planes%xmax = MAX(planes%xmax,grid(i)%cntr(1))
         planes%ymin = MIN(planes%ymin,grid(i)%cntr(2))
         planes%ymax = MAX(planes%ymax,grid(i)%cntr(2))
         planes%zmin = MIN(planes%zmin,grid(i)%cntr(3))
         planes%zmax = MAX(planes%zmax,grid(i)%cntr(3))
      END DO

   END SUBROUTINE get_boundary_planes

!-------------------------------------------------------------------------------

   SUBROUTINE get_closest_approach(RHS,planes,gap)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)  :: RHS(:)
      TYPE(fmm_planes),   INTENT(IN)  :: planes
      REAL(REALK),        INTENT(OUT) :: gap

      INTEGER(INTK) :: i

      gap = 1.0d10
      DO i = 1, SIZE(RHS)
         gap = MIN(gap, ABS(RHS(i)%cntr(1) - planes%xmin) )
         gap = MIN(gap, ABS(RHS(i)%cntr(1) - planes%xmax) )
         gap = MIN(gap, ABS(RHS(i)%cntr(2) - planes%ymin) )
         gap = MIN(gap, ABS(RHS(i)%cntr(2) - planes%ymax) )
         gap = MIN(gap, ABS(RHS(i)%cntr(3) - planes%zmin) )
         gap = MIN(gap, ABS(RHS(i)%cntr(3) - planes%zmax) )
      END DO

   END SUBROUTINE get_closest_approach

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_opt_near_field(scheme,LHS,RHS)

      USE fmm_box_utils, ONLY: fmm_deepest_level, fmm_branch, fmm_grain

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: LHS(:)
      TYPE(raw_mm_paras), INTENT(IN)    :: RHS(:)
      TYPE(scheme_paras), INTENT(INOUT) :: scheme

      TYPE(fmm_planes) :: planes
      REAL(REALK)      :: gap, dummy
      REAL(REALK)      :: grain
      INTEGER(INTK)    :: branch

      ! We have only considered branch-free algorithm
      IF (.NOT. scheme%branch_free) RETURN

      CALL get_boundary_planes(LHS,planes)
      CALL verify_planes(LHS,planes)
      CALL get_closest_approach(RHS,planes,gap)
      WRITE(LUPRI,'(A,E15.7)') ' Minimum distance to boundary =', gap

      IF (gap < scheme%extent_min) THEN
         CALL fmm_quit('conflict between branch-free radius and boundary gap!')
      END IF

      grain  = fmm_grain(scheme,fmm_deepest_level(scheme))
      branch = fmm_branch(dummy,one/grain)

      ! we choose this condition to skip the NF conservatively
      IF (gap > (branch+2)*grain) THEN
         WRITE(LUPRI,*) 'There are no near-field interactions!'
         scheme%include_near_field = .FALSE.
      END IF

   END SUBROUTINE fmm_opt_near_field

!-------------------------------------------------------------------------------

END MODULE fmm_boundary
