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

module fmm_qlm_builder

use fmm_global_paras, only: INTK, REALK, LUPRI, LUINTM, fmm_counters, scheme_paras, raw_mm_data, id_node, raw_mm_paras, &
                            ELECTRONIC_ONLY,NUCLEAR_ONLY,Zero,One
use fmm_stats, only: stat_n_basis, stat_raw_moms_LHS, stat_raw_moms_RHS
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_get_raw_qlm, &
          fmm_deallocate_qlm

! Public variables
public :: fmm_system_size, fmm_coord_shift

! Coordinate translation required to ensure all primitive (x,y,z) positive
real(REALK), save :: fmm_coord_shift(3)

! Number of multipole moments written to interface file
type(fmm_counters), save :: n_mms
! Number of AO (contracted) basis functions
integer(INTK), save :: nbas
! (minimum) dimension of a cube which encloses all moments
real(REALK), save :: fmm_system_size

contains

!-------------------------------------------------------------------------------

subroutine fmm_get_raw_qlm(scheme,dens,LHS,RHS)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  real(REALK), intent(in)        :: dens(:,:)
  type(raw_mm_data), intent(out) :: LHS, RHS
  type(raw_mm_data) :: mm_data
  integer(INTK) :: LMAX

  LMAX = scheme%raw_LMAX

  ! Get all multipole data from program interface
  call fmm_get_n_mms_from_file(LMAX)
  call fmm_allocate_mms_arrays(LMAX,n_mms%tot,mm_data)
  call fmm_read_in_raw_data(dens,mm_data)
  call fmm_get_system_size_and_shift(mm_data%paras)
  if (scheme%branch_free) call fmm_make_branch_free_extents(scheme,mm_data)

  ! Assign LHS range multipole data
  call fmm_distribute_LHS_RHS_data(LMAX,scheme%LHS_mm_range,mm_data,LHS)
  stat_raw_moms_LHS = size(LHS%paras)
  ! Assign RHS range multipole data
  call fmm_distribute_LHS_RHS_data(LMAX,scheme%RHS_mm_range,mm_data,RHS)
  stat_raw_moms_RHS = size(RHS%paras)

  call fmm_deallocate_mms_arrays(mm_data)

end subroutine fmm_get_raw_qlm

!-------------------------------------------------------------------------------
! Branch-free scheme basically uses CFMM code but we simplify all the
! extents such that they are all the same and only _very_ close interactions
! are avoided (for numerical stability)

subroutine fmm_make_branch_free_extents(scheme,mm_data)

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  type(raw_mm_data), intent(inout) :: mm_data

  integer(INTK) :: i

  do i=1,n_mms%elec
    mm_data%paras(i)%ext = scheme%extent_min
  end do

end subroutine fmm_make_branch_free_extents

!-------------------------------------------------------------------------------

subroutine fmm_distribute_LHS_RHS_data(LMAX,mm_range,all_data,sub_data)

  implicit none
  integer(INTK), intent(in)      :: LMAX, mm_range
  type(raw_mm_data), intent(in)  :: all_data
  type(raw_mm_data), intent(out) :: sub_data
  integer(INTK) :: foo, i, hi, lo

  lo = 1
  hi = n_mms%tot
  if (mm_range == NUCLEAR_ONLY) lo = 1+n_mms%elec
  if (mm_range == ELECTRONIC_ONLY) hi = n_mms%elec
  foo = hi-lo+1

  call fmm_allocate_mms_arrays(LMAX,foo,sub_data)

  sub_data%qlm(:,:) = all_data%qlm(:,lo:hi)
  sub_data%dens(:) = all_data%dens(lo:hi)
  sub_data%paras(:) = all_data%paras(lo:hi)
  sub_data%J_indices(:) = all_data%J_indices(lo:hi)
  nullify(sub_data%qlm_T)
  nullify(sub_data%qlm_W)

  ! Initialise parameter:moments mapping,
  ! and batch numbers for unique centres
  do i=1,foo
    sub_data%paras(i)%id = i
    sub_data%paras(i)%batch = i
  end do

end subroutine fmm_distribute_LHS_RHS_data

!-------------------------------------------------------------------------------

subroutine fmm_deallocate_qlm(LHS,RHS)

  implicit none
  type(raw_mm_data), intent(out) :: LHS, RHS

  call fmm_deallocate_mms_arrays(LHS)
  call fmm_deallocate_mms_arrays(RHS)

end subroutine fmm_deallocate_qlm

!-------------------------------------------------------------------------------

subroutine fmm_get_n_mms_from_file(LMAX_in)

  implicit none
  integer(INTK), intent(in) :: LMAX_in
  integer(INTK) :: LMAX !, ndim
  integer(INTK), external :: IsFreeUnit

  ! Read number of electronic moments
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file='multipoles.fmm1header',status='OLD',action='READ',form='UNFORMATTED')
  rewind(LUINTM)
  read(LUINTM) LMAX,nbas,n_mms%elec
  close(unit=LUINTM,status='KEEP')

  if (LMAX /= LMAX_in) then
    write(LUPRI,*) LMAX,LMAX_in
    call fmm_quit('LMAX inconsistency in MM interface!')
  end if
  if (n_mms%elec < 1) call fmm_quit('No moments generated!')
  ! This test only works for moments build over contracted AO pairs
  !ndim = nbas*(nbas+1)/2
  !if (n_mms%elec > ndim) then
  !  write(LUPRI,*) LMAX,nbas,ndim,n_mms%elec
  !  call fmm_quit('Too many moments generated, based on AO number!')
  !end if

  ! Read number of nuclear moments or potential grid points
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file='multipoles.fmm2header',status='OLD',action='READ',form='UNFORMATTED')
  rewind(LUINTM)
  read(LUINTM) n_mms%nuc
  close(unit=LUINTM,status='KEEP')

  n_mms%tot = n_mms%elec+n_mms%nuc
  stat_n_basis = nbas

end subroutine fmm_get_n_mms_from_file

!-------------------------------------------------------------------------------

subroutine fmm_allocate_mms_arrays(LMAX,ndim,mm_data)

  implicit none
  integer(INTK), intent(in)      :: LMAX, ndim
  type(raw_mm_data), inteNt(out) :: mm_data
  integer(INTK) :: i, foo, alloc_error

  nullify(mm_data%paras,mm_data%dens,mm_data%batch_map)
  nullify(mm_data%qlm,mm_data%qlm_W,mm_data%qlm_T)

  allocate(mm_data%paras(ndim))
  allocate(mm_data%J_indices(ndim))

  ! Initialise parameters
  do i=1,ndim
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
  end do

  allocate(mm_data%dens(ndim))
  foo = (LMAX+1)**2
  write(LUPRI,*) 'mms_arrays: Attempting to allocate',max(1,foo*ndim*8/1000000),'MB of memory...'
  allocate(mm_data%qlm(foo,ndim),stat=alloc_error)
  if (alloc_error /= 0) write(LUPRI,*) '... Failed!'

  mm_data%qlm(:,:) = zero  ! only non-zero written explicitly

end subroutine fmm_allocate_mms_arrays

!------------------------------------------------------------------------------

subroutine fmm_deallocate_mms_arrays(mm_data)

  implicit none
  type(raw_mm_data), intent(out) :: mm_data
  integer(INTK) :: i

  if (associated(mm_data%paras)) deallocate(mm_data%paras)
  if (associated(mm_data%J_indices)) deallocate(mm_data%J_indices)
  if (associated(mm_data%dens)) deallocate(mm_data%dens)
  if (associated(mm_data%qlm)) deallocate(mm_data%qlm)
  if (associated(mm_data%qlm_T)) deallocate(mm_data%qlm_T)
  if (associated(mm_data%qlm_W)) deallocate(mm_data%qlm_W)

  if (associated(mm_data%batch_map)) then
    do i=1,size(mm_data%batch_map)
      call free_batch_map(mm_data%batch_map(i)%head)
    end do
  end if
  if (associated(mm_data%batch_map)) deallocate(mm_data%batch_map)

  nullify(mm_data%paras,mm_data%dens,mm_data%batch_map)
  nullify(mm_data%qlm,mm_data%qlm_W,mm_data%qlm_T)
  nullify(mm_data%J_indices)

  contains

  recursive subroutine free_batch_map(node)

    implicit none
    type(id_node), pointer :: node
    if (associated(node%next)) then
      call free_batch_map(node%next)
    end if
    deallocate(node)
    nullify(node)

  end subroutine free_batch_map

end subroutine fmm_deallocate_mms_arrays

!-------------------------------------------------------------------------------
! Read in multipole moment data from interface file.
! In all this MM code we assume the order of moments is:
!   (0),(-1,0,1),(-2,-1,0,1,2)... wrt (L,M)

subroutine fmm_read_in_raw_data(dens,mm_data)

# include "macros.fh"

  implicit none
  real(REALK), intent(in)        :: dens(:,:)
  type(raw_mm_data), intent(out) :: mm_data
  real(REALK)   :: PX, PY, PZ, SPH
  integer(INTK) :: I, J, L, M, A, B, LM, X
  integer(INTK), external :: IsFreeUnit

  ! Read electronic multipole moments into core
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file='multipoles.fmm1',status='OLD',action='READ',form='UNFORMATTED')
  rewind(LUINTM)

  readloop: do

    read(LUINTM) I,L,M,A,B,PX,PY,PZ,SPH
    ! EOF marked by record with negative angular momentum
    if (L < 0) exit readloop

    !if ((L == 0) .and. (ABS(SPH) > 1.0e-12_REALK)) write(LUPRI,'(5I4,1X,3F8.4,2E13.4)') I,L,M,A,B,PX,PY,PZ,SPH,dens(A,B)

    if (A > nbas) call fmm_quit('interface file error 0')
    if (B > nbas) call fmm_quit('interface file error 00')
    if (I < 1) call fmm_quit('interface file error 1')
    if (I > size(mm_data%qlm,2)) call fmm_quit('interface error 11')
    if (((L+1)**2) > size(mm_data%qlm,1)) call fmm_quit('interface 156')

    ! Indices to map moments to orbitals in J-matrix
    mm_data%J_indices(I)%i_indx = A
    mm_data%J_indices(I)%j_indx = B
    ! Multipole expansion centre
    mm_data%paras(I)%cntr = [PX,PY,PZ]
    ! See defn of 'extent' in p.424 MEST Helgaker et al
    mm_data%paras(I)%ext = 0
    LM = L*(L+1)+M+1
    mm_data%dens(I) = dens(A,B)
    ! Components (l,m) of MM without density factorised in
    mm_data%qlm(LM,I) = SPH

  end do readloop

  close(unit=LUINTM,status='KEEP')

  ! Next read nuclei data: charge and location
  ! This is also used for passing the grid points when
  ! computing an arbitrary potential
  if (n_mms%nuc == 0) return
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file='multipoles.fmm2',status='OLD',action='READ',form='UNFORMATTED')
  rewind(LUINTM)

  do J=1,n_mms%nuc
    I = n_mms%elec+J
    read(LUINTM) X,L,M,A,B,PX,PY,PZ,SPH
    !write(LUPRI,'(I4,1X,3F8.4,2E13.4)') I,PX,PY,PZ,SPH
    mm_data%qlm(1,I) = SPH
    mm_data%paras(I)%cntr = [PX,PY,PZ]
  end do

  unused_var(X)

  mm_data%dens((n_mms%elec+1):) = one
  mm_data%paras((n_mms%elec+1):)%ext = zero         ! point charges
  mm_data%J_indices((n_mms%elec+1):)%i_indx = 0     ! not relevant
  mm_data%J_indices((n_mms%elec+1):)%j_indx = 0     ! not relevant

  close(unit=LUINTM,status='KEEP')

end subroutine fmm_read_in_raw_data

!-------------------------------------------------------------------------------

subroutine fmm_get_system_size_and_shift(paras)

  implicit none
  type(raw_mm_paras), intent(inout) :: paras(:)
  real(REALK)   :: sys_min(3), sys_max(3)
  integer(INTK) :: i

  sys_min = paras(1)%cntr
  sys_max = paras(1)%cntr
  do i=1,size(paras)
    sys_min(:) = min(sys_min(:),paras(i)%cntr(:))
    sys_max(:) = max(sys_max(:),paras(i)%cntr(:))
  end do
  fmm_system_size = maxval(sys_max-sys_min)
  if (fmm_system_size == zero) call fmm_quit('zero system size!')

  fmm_coord_shift = sys_min
  !do i=1,size(paras)
  !  paras(i)%cntr(:) = paras(i)%cntr(:)-sys_min(:)
  !end do

end subroutine fmm_get_system_size_and_shift

!-------------------------------------------------------------------------------

end module fmm_qlm_builder
