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
MODULE fmm_box_packer

   USE fmm_global_paras, ONLY: INTK, REALK, scheme_paras, raw_mm_paras, box_mm_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_pkd_paras,              &
             fmm_shift_and_pack_paras

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_pkd_paras(deepest_level,scheme,raw_paras,pkd_paras)

      USE fmm_box_utils, ONLY: fmm_box_centre, fmm_parent_box, fmm_grain

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)    :: deepest_level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_paras), INTENT(INOUT) :: raw_paras(:)
      TYPE(box_mm_paras), POINTER       :: pkd_paras(:)

      TYPE(box_mm_paras) :: tmp_paras(SIZE(raw_paras))
      INTEGER(INTK)      :: i, l_up, box(3), foo
      INTEGER(INTK)      :: tmp_map(SIZE(raw_paras))
      REAL(REALK)        :: grain

      l_up = deepest_level-1
      grain = fmm_grain(scheme,l_up)
      DO i = 1, SIZE(raw_paras)
         tmp_paras(i)%box(:) = raw_paras(i)%box(:)
         tmp_paras(i)%cntr(:) = raw_paras(i)%box_cntr(:)
         tmp_paras(i)%bra = raw_paras(i)%bra
         tmp_paras(i)%level = deepest_level
         tmp_paras(i)%id = i
         box(:) = fmm_parent_box(tmp_paras(i)%box(:))
         tmp_paras(i)%cntr_up(:) = fmm_box_centre(box,grain)
         tmp_paras(i)%map_up = 0   ! not defined yet
      END DO

      IF (ASSOCIATED(pkd_paras)) CALL fmm_quit ('paras should be nullified!')
      IF (.FALSE.) THEN
         ! just build unpacked paras
         foo = SIZE(raw_paras)
         ALLOCATE (pkd_paras(foo))
         pkd_paras(:) = tmp_paras(:)
         tmp_map = (/ (i,i=1,SIZE(raw_paras)) /)
      ELSE
         ! combine paras in same box and branch (i.e. packed form)
         CALL pack_boxed_paras(tmp_paras,pkd_paras,tmp_map)
      END IF
      ! store map of packing indices (raw:boxed)
      raw_paras(:)%map_up = tmp_map(:)

   END SUBROUTINE fmm_init_pkd_paras

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_shift_and_pack_paras(level,scheme,paras_in,paras_out)

      USE fmm_box_utils, ONLY: fmm_box_centre, fmm_grain,              &
                               fmm_parent_box, fmm_parent_bra

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)    :: level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(box_mm_paras), INTENT(INOUT) :: paras_in(:)
      TYPE(box_mm_paras), POINTER       :: paras_out(:)

      TYPE(box_mm_paras) :: tmp_paras(SIZE(paras_in))
      INTEGER(INTK)      :: i, l_up, box(3), foo
      INTEGER(INTK)      :: tmp_map(SIZE(paras_in))
      REAL(REALK)        :: grain, grain_up

      ! build tmp array for unpacked paras at next level up
      l_up = level-1
      grain = fmm_grain(scheme,level)
      grain_up = fmm_grain(scheme,l_up)
      DO i = 1, SIZE(paras_in)
         tmp_paras(i)%box = fmm_parent_box(paras_in(i)%box)
         tmp_paras(i)%cntr(:) = fmm_box_centre(tmp_paras(i)%box,grain)
         tmp_paras(i)%bra = fmm_parent_bra(paras_in(i)%bra)
         tmp_paras(i)%level = level
         tmp_paras(i)%id = i
         box(:) = fmm_parent_box(tmp_paras(i)%box(:))
         tmp_paras(i)%cntr_up(:) = fmm_box_centre(box,grain_up)
         tmp_paras(i)%map_up = 0   ! not defined yet
      END DO

      IF (.FALSE.) THEN
         ! just build unpacked paras
         foo = SIZE(paras_in)
         ALLOCATE (paras_out(foo))
         paras_out(:) = tmp_paras(:)
         tmp_map = (/ (i,i=1,SIZE(paras_in)) /)
      ELSE
         ! combine paras in same box and branch (i.e. packed form)
         CALL pack_boxed_paras(tmp_paras,paras_out,tmp_map)
      END IF
      ! store map of packing indices (raw:boxed)
      paras_in(:)%map_up = tmp_map(:)

   END SUBROUTINE fmm_shift_and_pack_paras

!-------------------------------------------------------------------------------

   SUBROUTINE pack_boxed_paras(paras_in,paras_out,map)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: paras_in(:)
      TYPE(box_mm_paras), POINTER       :: paras_out(:)
      INTEGER(INTK),      INTENT(OUT)   :: map(:)

      TYPE(box_mm_paras) :: tmp_paras(SIZE(paras_in))
      INTEGER(INTK) :: i, k

      CALL fmm_sort_wrt_boxes_and_branches(1_INTK,paras_in)

      map(:) = 0
      tmp_paras(1) = paras_in(1)
      tmp_paras(1)%id = 1
      map(paras_in(1)%id) = 1

      k = 1
      DO i = 2, SIZE(paras_in)

         IF (paras_in(i)%box(3) /= paras_in(i-1)%box(3)) THEN
            k = k+1
         ELSE IF (paras_in(i)%box(2) /= paras_in(i-1)%box(2)) THEN
            k = k+1
         ELSE IF (paras_in(i)%box(1) /= paras_in(i-1)%box(1)) THEN
            k = k+1
         ELSE IF (paras_in(i)%bra /= paras_in(i-1)%bra) THEN
            k = k+1
         END IF

         tmp_paras(k) = paras_in(i)
         tmp_paras(k)%id = k
         map(paras_in(i)%id) = k

      END DO

      ! store packed parameters in exactly allocated array
      ALLOCATE (paras_out(k))
      paras_out = tmp_paras(1:k)

   END SUBROUTINE pack_boxed_paras

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_sort_wrt_boxes_and_branches(xyz,paras)

      USE fmm_sort_paras, ONLY: fmm_quicksort_wrt_boxes,      &
                                fmm_quicksort_wrt_branches

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)    :: xyz
      TYPE(box_mm_paras), INTENT(INOUT) :: paras(:)

      INTEGER(INTK) :: i, lo, hi
      REAL(REALK)   :: q1,q2

      IF (SIZE(paras) == 1) RETURN

      ! sort only if needed
      q1 = paras(1)%box(xyz)
      DO i = 2, SIZE(paras)
         q2 = paras(i)%box(xyz)
         IF ( q2 < q1 ) THEN
            CALL fmm_quicksort_wrt_boxes(paras,xyz)
            EXIT
         END IF
         q1 = q2
      END DO

      ! sub-sort next box component
      lo = 1
      DO i = 2, SIZE(paras)
         q1 = paras(i-1)%box(xyz)
         q2 = paras(i)%box(xyz)
         IF ( q2 /= q1 ) THEN
            hi = i-1
            IF (xyz == 3) THEN
               CALL fmm_quicksort_wrt_branches(paras(lo:hi))
            ELSE
               CALL fmm_sort_wrt_boxes_and_branches(xyz+1_INTK,paras(lo:hi))
            END IF
            lo = i
         END IF
      END DO

      ! do last batch
      hi = SIZE(paras)
      IF (xyz == 3) THEN
         CALL fmm_quicksort_wrt_branches(paras(lo:hi))
      ELSE
         CALL fmm_sort_wrt_boxes_and_branches(xyz+1_INTK,paras(lo:hi))
      END IF

   END SUBROUTINE fmm_sort_wrt_boxes_and_branches

!-------------------------------------------------------------------------------

END MODULE fmm_box_packer
