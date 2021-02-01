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
MODULE fmm_qlm_utils

   USE fmm_global_paras, ONLY: INTK, REALK, raw_mm_paras, raw_mm_data, id_list, id_node, One, Two, Half
   USE fmm_stats, ONLY: stat_pkd_moms_LHS, stat_pkd_moms_rHS, stat_screened_moms_RHS
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_renormalise_qlm,             &
             fmm_sort_paras_wrt_centre,       &
             fmm_assign_batches,              &
             fmm_factor_in_dens,              &
             fmm_get_T_sym_qlm,               &
             fmm_pack_raw_moments,            &
             fmm_pack_raw_parameters

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_renormalise_qlm(LMAX,qlm)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)    :: LMAX
      REAL(REALK),   INTENT(INOUT) :: qlm(:,:)

      INTEGER(INTK) :: i,L,M,p,pp
      REAL(REALK)   :: pref

      ! prefactor to symmetrize T-matrix
      DO i = 1, SIZE(qlm,2)
         DO L = 0, LMAX
            pp = L*(L+1) +1
            DO M = -L, -1
               pref = -one/(SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m)))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
            pref = one/FACTORIAL(L)
            p = pp  ! M=0
            qlm(p,i) = pref*qlm(p,i)
            DO M = 1, L
               pref = ((-1)**m)/SQRT(two*FACTORIAL(l-m)*FACTORIAL(l+m))
               p = pp+M
               qlm(p,i) = pref*qlm(p,i)
            END DO
         END DO
      END DO

   CONTAINS

      REAL(REALK) FUNCTION FACTORIAL(n)
         IMPLICIT NONE
         INTEGER(INTK), INTENT(IN) :: n
         INTEGER(INTK) :: i
         FACTORIAL = 1
         DO i = n, 2, -1
            FACTORIAL = FACTORIAL*i
         END DO
      END FUNCTION FACTORIAL

   END SUBROUTINE fmm_renormalise_qlm

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE fmm_sort_paras_wrt_centre(xyz,paras)

      USE fmm_sort_paras, ONLY: fmm_quicksort_wrt_vector

      IMPLICIT NONE
      INTEGER(INTK),      INTENT(IN)    :: xyz
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:)

      INTEGER(INTK) :: i, lo, hi
      REAL(REALK)   :: q1,q2

      IF (SIZE(paras) == 1) RETURN

      ! sort only if needed
      q1 = paras(1)%cntr(xyz)
      DO i = 2, SIZE(paras)
         q2 = paras(i)%cntr(xyz)
         IF ( q2 < q1 ) THEN
            CALL fmm_quicksort_wrt_vector(paras,xyz)
            EXIT
         END IF
         q1 = q2
      END DO

      ! sub-sort next cartesian component
      lo = 1
      DO i = 2, SIZE(paras)
         q1 = paras(i-1)%cntr(xyz)
         q2 = paras(i)%cntr(xyz)
         IF ( q2 /= q1 ) THEN
            hi = i-1
            IF (xyz == 3) THEN
               RETURN
            ELSE
               CALL fmm_sort_paras_wrt_centre(xyz+1_INTK,paras(lo:hi))
            END IF
            lo = i
         END IF
      END DO

      ! do last batch
      hi = SIZE(paras)
      IF (xyz == 3) THEN
         RETURN
      ELSE
         CALL fmm_sort_paras_wrt_centre(xyz+1_INTK,paras(lo:hi))
      END IF

   END SUBROUTINE fmm_sort_paras_wrt_centre

!-------------------------------------------------------------------------------
! Identify batches of moments with same centre assuming sorted order and
! assign batch index to parameter list

   SUBROUTINE fmm_assign_batches(paras)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:)

      INTEGER(INTK) :: i, batch
!FIXME
      REAL(REALK), PARAMETER :: TOLERANCE = 1.0D-20

      batch = 1
      paras(1)%batch = batch
      DO i = 2, SIZE(paras)
         IF (paras(i)%cntr(3) - paras(i-1)%cntr(3) > TOLERANCE) THEN
            batch = batch + 1
         ELSE IF (paras(i)%cntr(2) - paras(i-1)%cntr(2) > TOLERANCE) THEN
            batch = batch + 1
         ELSE IF (paras(i)%cntr(1) - paras(i-1)%cntr(1) > TOLERANCE) THEN
            batch = batch + 1
         END IF
         paras(i)%batch = batch
      END DO

   END SUBROUTINE fmm_assign_batches

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_factor_in_dens(dens,qlm)

      IMPLICIT NONE
      REAL(REALK), INTENT(IN)    :: dens(:)
      REAL(REALK), INTENT(INOUT) :: qlm(:,:)

      INTEGER(INTK) :: i

      DO i = 1, SIZE(qlm,2)
         qlm(:,i) = qlm(:,i)*dens(i)
      END DO

   END SUBROUTINE fmm_factor_in_dens

!-------------------------------------------------------------------------------
! Prefactorising moments to symmetrize modified T-matrix

   SUBROUTINE fmm_get_T_sym_qlm(LMAX,qlm_in,qlm_out)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: qlm_in(:,:)
      REAL(REALK),   INTENT(OUT) :: qlm_out(:,:)

      INTEGER(INTK) :: i,L,u, hi,lo
      REAL(REALK)   :: pref

      DO i = 1, SIZE(qlm_in,2)
         DO L = 0, LMAX
            u = L*(L+1) +1     ! m=0
            hi = u+L
            lo = u-L
            pref = two*((-1)**L)
            qlm_out(lo:hi,i) = pref*qlm_in(lo:hi,i)
            qlm_out(u,i) = half*pref*qlm_in(u,i)
         END DO
      END DO

   END SUBROUTINE fmm_get_T_sym_qlm

!-------------------------------------------------------------------------------
! Get number of unique batches;
! Also check raw_paras are sorted by batch;
! Do this by checking the batch ID is always increasing

   SUBROUTINE get_nbatch(paras,nbatch)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)  :: paras(:)
      INTEGER(INTK),      INTENT(OUT) :: nbatch
      INTEGER(INTK) :: i, ndim

      ndim = SIZE(paras)
      nbatch = 1
      DO i = 2, ndim
         IF ( paras(i)%batch < paras(i-1)%batch ) THEN
            CALL fmm_quit('batches of packed moments not sorted!')
         END IF
         ! some batch numbers are skipped so need to take care of this
         IF ( paras(i)%batch /= paras(i-1)%batch ) THEN
            nbatch = nbatch + 1
         END IF
      END DO

   END SUBROUTINE get_nbatch

!-------------------------------------------------------------------------------

   SUBROUTINE get_pkd_data(add_dens,raw_mms,pkd_paras,pkd_qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)  :: add_dens
      TYPE(raw_mm_data),  INTENT(IN)  :: raw_mms
      TYPE(raw_mm_paras), INTENT(OUT) :: pkd_paras(:)
      REAL(REALK),        INTENT(OUT) :: pkd_qlm(:,:)
      INTEGER(INTK) :: i,j,k, last_batch

      j = 0
      last_batch = -1
      DO i = 1, SIZE(raw_mms%paras)
         ! Use mapping of raw parameters to raw moments
         k = raw_mms%paras(i)%id
         IF ( raw_mms%paras(i)%batch == last_batch ) THEN
            ! Element same batch as previous
            IF (add_dens) THEN
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,k)*raw_mms%dens(k)
            ELSE
               pkd_qlm(:,j) = pkd_qlm(:,j) + raw_mms%qlm(:,k)
            END IF
         ELSE
            ! Element in new batch
            j = j + 1
            pkd_paras(j) = raw_mms%paras(i)
            IF (add_dens) THEN
               pkd_qlm(:,j) = raw_mms%qlm(:,k)*raw_mms%dens(k)
            ELSE
               pkd_qlm(:,j) = raw_mms%qlm(:,k)
            END IF
         END IF
         last_batch = raw_mms%paras(i)%batch
      END DO

   END SUBROUTINE get_pkd_data

!-------------------------------------------------------------------------------

   SUBROUTINE get_screened_nmom(screen_thr,qlm,skip,nskip)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN)  :: screen_thr
      REAL(REALK),   INTENT(IN)  :: qlm(:,:)
      LOGICAL,       INTENT(OUT) :: skip(:)
      INTEGER(INTK), INTENT(OUT) :: nskip
      INTEGER(INTK) :: i,j

      nskip = 0
      batches: DO i = 1, SIZE(qlm,2)
         skip(i) = .TRUE.
         nskip = nskip + 1
         lm_loop: DO j = 1, SIZE(qlm,1)
            IF ( ABS(qlm(j,i)) > screen_thr ) THEN
               skip(i) = .FALSE.
               nskip = nskip - 1
               EXIT lm_loop
            END IF
         END DO lm_loop
      END DO batches

   END SUBROUTINE get_screened_nmom

!-------------------------------------------------------------------------------
! Here we squeeze all the significant batches of moments to be sequential
! at the top of the array, with all the insignificant moments overwritten;
! we use skip(:) as a logical mask to direct the skipping.

   SUBROUTINE squeeze_significant_batches(skip,paras,qlm)

      IMPLICIT NONE
      LOGICAL,            INTENT(IN)    :: skip(:)
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:)
      REAL(REALK),        INTENT(INOUT) :: qlm(:,:)
      INTEGER(INTK) :: i, j

      IF (SIZE(paras) /= SIZE(qlm,2)) STOP 'paras and qlm should be same size!'
      IF (SIZE(paras) /= SIZE(skip)) STOP 'paras and skip should be same size!'

      j = 0
      DO i = 1, SIZE(paras)
         IF (skip(i)) CYCLE
         j = j + 1
         paras(j) = paras(i)
         qlm(:,j) = qlm(:,i)
      END DO

   END SUBROUTINE squeeze_significant_batches

!-------------------------------------------------------------------------------
! Routine to pack a set of raw moments by batches,
! where members of a batch share a common centre and extent,
! including density factoring and screening if requested.
! Note that no record is kept of the packing for later "unpacking".

   SUBROUTINE fmm_pack_raw_moments(mm_data,dens,dens_thr)

      IMPLICIT NONE
      TYPE(raw_mm_data),  INTENT(INOUT) :: mm_data
      LOGICAL,            INTENT(IN)    :: dens
      REAL(REALK),        INTENT(IN)    :: dens_thr

      TYPE(raw_mm_paras), ALLOCATABLE :: pkd_paras(:)
      REAL(REALK), ALLOCATABLE :: pkd_qlm(:,:)
      LOGICAL, ALLOCATABLE :: skip(:)
      INTEGER(INTK) :: nbatch, nskip, foo, ndim

      ! Get number of unique batches
      CALL get_nbatch(mm_data%paras,nbatch)

      ! Now build packed data by summing raw data in the same batch
      ALLOCATE (pkd_paras(nbatch))
      foo = SIZE(mm_data%qlm,1)
      ALLOCATE (pkd_qlm(foo,nbatch))
      CALL get_pkd_data(dens,mm_data,pkd_paras,pkd_qlm)
      ndim = nbatch

      IF ( dens ) THEN
         ! Perform density-based screening;
         ! Note we do not reallocate whole array, but push all the
         ! significant terms to the top, and only this array section
         ! of significant moments is then pointed to.
         ALLOCATE (skip(nbatch))
         CALL get_screened_nmom(dens_thr,pkd_qlm,skip,nskip)
         ndim = nbatch - nskip  ! number of significant moments
         CALL squeeze_significant_batches(skip,pkd_paras,pkd_qlm)
         DEALLOCATE (skip)
      END IF

      ! This is just for statistics
      stat_pkd_moms_RHS = nbatch
      stat_screened_moms_RHS = ndim

      ! Reallocate the original moments and copy across the packed ones
      DEALLOCATE (mm_data%paras, mm_data%qlm)
      NULLIFY (mm_data%paras, mm_data%qlm)
      ALLOCATE (mm_data%paras(ndim), mm_data%qlm(foo,ndim))
      mm_data%paras(:) = pkd_paras(:ndim)
      mm_data%qlm(:,:) = pkd_qlm(:,:ndim)

      DEALLOCATE (pkd_paras, pkd_qlm)

   END SUBROUTINE fmm_pack_raw_moments

!-------------------------------------------------------------------------------
! Routine to drive the packing of a set of raw mm parameters by batches,
! where members of a batch share a common centre and extent,
! including the build of a mapping between the packed paras and the
! original raw paras.  This map can be used for later "unpacking".

   SUBROUTINE fmm_pack_raw_parameters(mm_data)

      IMPLICIT NONE
      TYPE(raw_mm_data), INTENT(INOUT) :: mm_data

      TYPE(raw_mm_paras), ALLOCATABLE :: pkd_paras(:)
      INTEGER(INTK) :: i,j, nbatch, last_batch

      ! Get number of unique batches;
      CALL get_nbatch(mm_data%paras,nbatch)
      stat_pkd_moms_LHS = nbatch

      ! Initialise packed paras and batch map
      ALLOCATE(pkd_paras(nbatch))
      ALLOCATE(mm_data%batch_map(nbatch))
      DO i = 1, nbatch
         mm_data%batch_map(i)%occ = 0
         NULLIFY( mm_data%batch_map(i)%head )
      END DO

      ! Now build packed paras by compressing raw paras in same batch
      j = 0
      last_batch = -1
      DO i = 1, SIZE(mm_data%paras)
         IF ( mm_data%paras(i)%batch == last_batch ) THEN
            ! Add raw parameter mapping to existing linked-list for this batch
            CALL add_batch_item(mm_data%batch_map(j),mm_data%paras(i)%id)
         ELSE
            ! Element in new batch
            j = j + 1
            pkd_paras(j) = mm_data%paras(i)
            ! Linked-list for this batch is empty, so start one
            mm_data%batch_map(j)%occ = 1
            ALLOCATE( mm_data%batch_map(j)%head )
            ! Maintain mapping between parameters and moments
            mm_data%batch_map(j)%head%id = mm_data%paras(i)%id
            NULLIFY( mm_data%batch_map(j)%head%next ) ! rest of list is empty
         END IF
         last_batch = mm_data%paras(i)%batch
      END DO

      ! Reallocate the original parameters and copy across the packed ones
      DEALLOCATE (mm_data%paras)
      NULLIFY (mm_data%paras)
      ALLOCATE (mm_data%paras(nbatch))
      mm_data%paras(:) = pkd_paras(:)

      DEALLOCATE (pkd_paras)

   CONTAINS

      SUBROUTINE add_batch_item(batch_list,raw_id)

         IMPLICIT NONE
         TYPE(id_list), INTENT(INOUT) :: batch_list
         INTEGER(INTK), INTENT(IN)    :: raw_id
         TYPE(id_node), POINTER :: new_node

         batch_list%occ = batch_list%occ + 1
         ALLOCATE( new_node )
         new_node%id = raw_id
         IF (ASSOCIATED(batch_list%head%next)) THEN
            ! More than one entry in list (including head)
            ! so point new_node to old second entry
            new_node%next => batch_list%head%next
            ! Point head to new_node
            NULLIFY(batch_list%head%next)
            batch_list%head%next => new_node
         ELSE
            ! Only head so far; make new_node our second entry
            batch_list%head%next => new_node
            NULLIFY(new_node%next)   ! end of list
         END IF

      END SUBROUTINE add_batch_item

   END SUBROUTINE fmm_pack_raw_parameters

!-------------------------------------------------------------------------------

END MODULE fmm_qlm_utils
