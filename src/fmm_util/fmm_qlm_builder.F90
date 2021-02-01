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
MODULE fmm_qlm_builder

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, LUINTM, fmm_counters, scheme_paras, raw_mm_data, id_node, raw_mm_paras, &
                               ELECTRONIC_ONLY, NUCLEAR_ONLY, Zero, One
   USE fmm_stats, ONLY: stat_n_basis, stat_raw_moms_LHS, stat_raw_moms_RHS
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_raw_qlm,             &
             fmm_deallocate_qlm

   ! Public variables
   PUBLIC :: fmm_system_size, fmm_coord_shift

   ! Coordinate translation required to ensure all primitive (x,y,z) positive
   REAL(REALK), SAVE :: fmm_coord_shift(3)

   ! Number of multipole moments written to interface file
   TYPE(fmm_counters), SAVE :: n_mms
   ! Number of AO (contracted) basis functions
   INTEGER(INTK), SAVE :: nbas
   ! (minimum) dimension of a cube which encloses all moments
   REAL(REALK), SAVE :: fmm_system_size

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_raw_qlm(scheme,dens,LHS,RHS)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      REAL(REALK),        INTENT(IN)  :: dens(:,:)
      TYPE(raw_mm_data),  INTENT(OUT) :: LHS, RHS
      TYPE(raw_mm_data) :: mm_data
      INTEGER(INTK) :: LMAX

      LMAX = scheme%raw_LMAX

      ! Get all multipole data from program interface
      CALL fmm_get_n_mms_from_file(LMAX)
      CALL fmm_allocate_mms_arrays(LMAX,n_mms%tot,mm_data)
      CALL fmm_read_in_raw_data(dens,mm_data)
      CALL fmm_get_system_size_and_shift(mm_data%paras)
      IF (scheme%branch_free) CALL fmm_make_branch_free_extents(scheme,mm_data)

      ! Assign LHS range multipole data
      CALL fmm_distribute_LHS_RHS_data(LMAX,scheme%LHS_mm_range,mm_data,LHS)
      stat_raw_moms_LHS = SIZE(LHS%paras)
      ! Assign RHS range multipole data
      CALL fmm_distribute_LHS_RHS_data(LMAX,scheme%RHS_mm_range,mm_data,RHS)
      stat_raw_moms_RHS = SIZE(RHS%paras)

      CALL fmm_deallocate_mms_arrays(mm_data)

   END SUBROUTINE fmm_get_raw_qlm

!-------------------------------------------------------------------------------
! Branch-free scheme basically uses CFMM code but we simplify all the
! extents such that they are all the same and only _very_ close interactions
! are avoided (for numerical stability)

   SUBROUTINE fmm_make_branch_free_extents(scheme,mm_data)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(INOUT) :: mm_data

      INTEGER(INTK) :: i

      DO i = 1, n_mms%elec
         mm_data%paras(i)%ext = scheme%extent_min
      END DO

   END SUBROUTINE fmm_make_branch_free_extents

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_distribute_LHS_RHS_data(LMAX,mm_range,all_data,sub_data)

      IMPLICIT NONE
      INTEGER(INTK),     INTENT(IN)  :: LMAX, mm_range
      TYPE(raw_mm_data), INTENT(IN)  :: all_data
      TYPE(raw_mm_data), INTENT(OUT) :: sub_data
      INTEGER(INTK) :: foo, i, hi, lo

      lo = 1
      hi = n_mms%tot
      IF (mm_range == NUCLEAR_ONLY)    lo = 1+n_mms%elec
      IF (mm_range == ELECTRONIC_ONLY) hi = n_mms%elec
      foo = hi - lo + 1

      CALL fmm_allocate_mms_arrays(LMAX,foo,sub_data)

      sub_data%qlm(:,:)     = all_data%qlm(:,lo:hi)
      sub_data%dens(:)      = all_data%dens(lo:hi)
      sub_data%paras(:)     = all_data%paras(lo:hi)
      sub_data%J_indices(:) = all_data%J_indices(lo:hi)
      NULLIFY(sub_data%qlm_T)
      NULLIFY(sub_data%qlm_W)

      ! Initialise parameter:moments mapping,
      ! and batch numbers for unique centres
      DO i = 1, foo
         sub_data%paras(i)%id = i
         sub_data%paras(i)%batch = i
      END DO

   END SUBROUTINE fmm_distribute_LHS_RHS_data

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_deallocate_qlm(LHS,RHS)

      IMPLICIT NONE
      TYPE(raw_mm_data), INTENT(OUT) :: LHS, RHS

      CALL fmm_deallocate_mms_arrays(LHS)
      CALL fmm_deallocate_mms_arrays(RHS)

   END SUBROUTINE fmm_deallocate_qlm

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_n_mms_from_file(LMAX_in)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX_in
      INTEGER(INTK) :: LMAX !, ndim
      INTEGER(INTK), EXTERNAL :: IsFreeUnit

      ! Read number of electronic moments
      LUINTM = IsFreeUnit(LUINTM)
      OPEN (UNIT=LUINTM,FILE='multipoles.fmm1header',  &
            STATUS='OLD',ACTION='READ',FORM='UNFORMATTED')
      REWIND (LUINTM)
      READ (LUINTM) LMAX, nbas, n_mms%elec
      CLOSE(UNIT=LUINTM,STATUS='KEEP')

      IF (LMAX /= LMAX_in) THEN
         WRITE(LUPRI,*) LMAX, LMAX_in
         CALL fmm_quit('LMAX inconsistency in MM interface!')
      END IF
      IF (n_mms%elec < 1) CALL fmm_quit('No moments generated!')
!      This test only works for moments build over contracted AO pairs
!      ndim = nbas*(nbas+1)/2
!      IF ( n_mms%elec > ndim ) THEN
!         WRITE(LUPRI,*) LMAX, nbas, ndim, n_mms%elec
!         CALL fmm_quit('Too many moments generated, based on AO number!')
!      END IF

      ! Read number of nuclear moments or potential grid points
      LUINTM = IsFreeUnit(LUINTM)
      OPEN (UNIT=LUINTM,FILE='multipoles.fmm2header',  &
            STATUS='OLD',ACTION='READ',FORM='UNFORMATTED')
      REWIND (LUINTM)
      READ (LUINTM) n_mms%nuc
      CLOSE(UNIT=LUINTM,STATUS='KEEP')

      n_mms%tot = n_mms%elec + n_mms%nuc
      stat_n_basis = nbas

   END SUBROUTINE fmm_get_n_mms_from_file

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_allocate_mms_arrays(LMAX,ndim,mm_data)

      IMPLICIT NONE
      INTEGER(INTK),     INTENT(IN)  :: LMAX, ndim
      TYPE(raw_mm_data), INTENT(OUT) :: mm_data
      INTEGER(INTK) :: i, foo, alloc_error

      NULLIFY (mm_data%paras, mm_data%dens, mm_data%batch_map)
      NULLIFY (mm_data%qlm, mm_data%qlm_W, mm_data%qlm_T)

      ALLOCATE (mm_data%paras(ndim))
      ALLOCATE (mm_data%J_indices(ndim))

      ! Initialise parameters
      DO i = 1, ndim
         mm_data%paras(i)%cntr = zero
         mm_data%paras(i)%ext = zero
         mm_data%paras(i)%id = 0
         mm_data%paras(i)%batch = 0
         mm_data%paras(i)%map_up = 0
         mm_data%paras(i)%box = 0
         mm_data%paras(i)%bra = 0
         mm_data%paras(i)%box_cntr = zero
         mm_data%J_indices(i)%i_indx = 0
         mm_data%J_indices(i)%j_indx = 0
      END DO

      ALLOCATE (mm_data%dens(ndim))
      foo = (LMAX+1)**2
      WRITE(LUPRI,*) 'mms_arrays: Attempting to allocate',  &
                     MAX(1,foo*ndim*8/1000000), 'MB of memory...'
      ALLOCATE (mm_data%qlm(foo,ndim), STAT=alloc_error)
      IF (alloc_error /= 0) WRITE(LUPRI,*) '... Failed!'

      mm_data%qlm(:,:) = zero  ! only non-zero written explicitly

   END SUBROUTINE fmm_allocate_mms_arrays

!------------------------------------------------------------------------------

   SUBROUTINE fmm_deallocate_mms_arrays(mm_data)

      IMPLICIT NONE
      TYPE(raw_mm_data), INTENT(OUT) :: mm_data
      INTEGER(INTK) :: i

      IF (ASSOCIATED(mm_data%paras)) DEALLOCATE (mm_data%paras)
      IF (ASSOCIATED(mm_data%J_indices)) DEALLOCATE (mm_data%J_indices)
      IF (ASSOCIATED(mm_data%dens)) DEALLOCATE (mm_data%dens)
      IF (ASSOCIATED(mm_data%qlm)) DEALLOCATE (mm_data%qlm)
      IF (ASSOCIATED(mm_data%qlm_T)) DEALLOCATE (mm_data%qlm_T)
      IF (ASSOCIATED(mm_data%qlm_W)) DEALLOCATE (mm_data%qlm_W)

      IF (ASSOCIATED(mm_data%batch_map)) THEN
         DO i = 1, SIZE(mm_data%batch_map)
            CALL free_batch_map(mm_data%batch_map(i)%head)
         END DO
      END IF
      IF (ASSOCIATED(mm_data%batch_map)) DEALLOCATE(mm_data%batch_map)

      NULLIFY (mm_data%paras, mm_data%dens, mm_data%batch_map)
      NULLIFY (mm_data%qlm, mm_data%qlm_W, mm_data%qlm_T)
      NULLIFY (mm_data%J_indices)

   CONTAINS

      RECURSIVE SUBROUTINE free_batch_map(node)

         IMPLICIT NONE
         TYPE(id_node), POINTER :: node
         IF (ASSOCIATED(node%next)) THEN
            CALL free_batch_map(node%next)
         END IF
         DEALLOCATE(node)
         NULLIFY(node)

      END SUBROUTINE free_batch_map

   END SUBROUTINE fmm_deallocate_mms_arrays

!-------------------------------------------------------------------------------
! Read in multipole moment data from interface file.
! In all this MM code we assume the order of moments is:
!   (0),(-1,0,1),(-2,-1,0,1,2)... wrt (L,M)

   SUBROUTINE fmm_read_in_raw_data(dens,mm_data)

      IMPLICIT NONE
      REAL(REALK),       INTENT(IN)  :: dens(:,:)
      TYPE(raw_mm_data), INTENT(OUT) :: mm_data
      REAL(REALK)   :: PX,PY,PZ, SPH
      INTEGER(INTK) :: I,J,L,M, A,B, LM, X
      INTEGER(INTK), EXTERNAL :: IsFreeUnit

      ! Read electronic multipole moments into core
      LUINTM = IsFreeUnit(LUINTM)
      OPEN (UNIT=LUINTM,FILE='multipoles.fmm1',STATUS='OLD',  &
            ACTION='READ',FORM='UNFORMATTED')
      REWIND (LUINTM)

      readloop: DO

         READ (LUINTM) I,L,M, A,B, PX,PY,PZ, SPH
         ! EOF marked by record with negative angular momentum
         IF (L < 0) EXIT readloop

!         IF ((L == 0) .AND. (ABS(SPH) > 1.0d-12) )   &
!         WRITE (LUPRI,'(5I4,1X,3F8.4,2E13.4)') I,L,M,A,B,PX,PY,PZ,SPH,dens(A,B)

         IF (A > nbas) CALL fmm_quit('interface file error 0')
         IF (B > nbas) CALL fmm_quit('interface file error 00')
         IF (I < 1) CALL fmm_quit('interface file error 1')
         IF (I > SIZE(mm_data%qlm,2)) CALL fmm_quit('interface error 11')
         IF (((L+1)**2) > SIZE(mm_data%qlm,1)) CALL fmm_quit('interface 156')

          ! Indices to map moments to orbitals in J-matrix
         mm_data%J_indices(I)%i_indx = A
         mm_data%J_indices(I)%j_indx = B
          ! Multipole expansion centre
         mm_data%paras(I)%cntr = (/ PX, PY, PZ /)
          ! See defn of 'extent' in p.424 MEST Helgaker et al
         mm_data%paras(I)%ext  = 0
         LM = L*(L+1) +M +1
         mm_data%dens(I) = dens(A,B)
          ! Components (l,m) of MM without density factorised in
         mm_data%qlm(LM,I) = SPH

      END DO readloop

      CLOSE(UNIT=LUINTM,STATUS='KEEP')

      ! Next read nuclei data: charge and location
      ! This is also used for passing the grid points when
      ! computing an arbitrary potential
      IF (n_mms%nuc == 0) RETURN
      LUINTM = IsFreeUnit(LUINTM)
      OPEN (UNIT=LUINTM,FILE='multipoles.fmm2',STATUS='OLD',  &
            ACTION='READ',FORM='UNFORMATTED')
      REWIND (LUINTM)

      DO J = 1, n_mms%nuc
         I = n_mms%elec + J
         READ (LUINTM) X,L,M, A,B, PX,PY,PZ, SPH
!         WRITE (LUPRI,'(I4,1X,3F8.4,2E13.4)') I, PX,PY,PZ, SPH
         mm_data%qlm(1,I) = SPH
         mm_data%paras(I)%cntr = (/ PX, PY, PZ /)
      END DO

      mm_data%dens((n_mms%elec+1):) = one
      mm_data%paras((n_mms%elec+1):)%ext = zero         ! point charges
      mm_data%J_indices((n_mms%elec+1):)%i_indx = 0     ! not relevant
      mm_data%J_indices((n_mms%elec+1):)%j_indx = 0     ! not relevant

      CLOSE(UNIT=LUINTM,STATUS='KEEP')

#ifdef _WARNING_WORKAROUND_
      IF (.FALSE.) CALL Unused_integer(X)
#endif
   END SUBROUTINE fmm_read_in_raw_data

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_system_size_and_shift(paras)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: paras(:)
      REAL(REALK)   :: sys_min(3), sys_max(3)
      INTEGER(INTK) :: i

      sys_min = paras(1)%cntr
      sys_max = paras(1)%cntr
      DO i = 1, SIZE(paras)
         sys_min(:) = MIN(sys_min(:),paras(i)%cntr(:))
         sys_max(:) = MAX(sys_max(:),paras(i)%cntr(:))
      END DO
      fmm_system_size = MAXVAL(sys_max - sys_min)
      IF (fmm_system_size == zero) CALL fmm_quit('zero system size!')

      fmm_coord_shift = sys_min
!      DO i = 1, SIZE(paras)
!         paras(i)%cntr(:) = paras(i)%cntr(:) - sys_min(:)
!      END DO

   END SUBROUTINE fmm_get_system_size_and_shift

!-------------------------------------------------------------------------------

END MODULE fmm_qlm_builder

