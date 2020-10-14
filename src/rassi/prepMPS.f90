!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2016-2017, Stefan Knecht                               *
!***********************************************************************
  subroutine prepMPS(            &
                     trorb,      &
                     istate,     &
                     lsym,       &
                     mplet,      &
                     mspro,      &
                     nacte,      &
                     tra,        &
                     ntra,       &
                     nish,       &
                     nash,       &
                     nosh,       &
                     nsym,       &
                     lupri,      &
! Leon 8/12/2016 -- "real" state index across all JobIphs, needed when we read checkpoint file names
                     istatereal, &
                     job,        &
                     ist         &
                     )

  !> module dependencies
#ifdef _DMRG_
  use qcmaquis_interface_cfg
  use qcmaquis_info
  use qcmaquis_interface_mpssi
#endif

  implicit none
!-------------------------------------------------------------------------------
!
!    driver routine for the MPS rotation wrt the orbital transformation matrix
!    given in TRA.
!    --> rotate MPS state ISTATE
!    NOTE: TRA contains square matrices, one per symmetry
!
!-------------------------------------------------------------------------------

  integer,intent(in)    :: istate
  integer,intent(in)    :: istatereal
  integer,intent(in)    :: job
  integer,intent(in)    :: ist
  integer,intent(in)    :: lsym
  integer,intent(in)    :: mplet
  integer,intent(in)    :: mspro
  integer,intent(in)    :: nacte
  integer,intent(in)    :: nsym
  integer,intent(in)    :: ntra
  integer,intent(in)    :: lupri
  integer,intent(in)    :: nish(nsym)
  integer,intent(in)    :: nash(nsym)
  integer,intent(in)    :: nosh(nsym)
  real*8, intent(inout) :: tra(ntra)
  logical,intent(in)    :: trorb
!-------------------------------------------------------------------------------
#ifdef _DMRG_
  integer               :: i, isym, no, ii, ista, jorb, ni
  real*8                :: fac(1,1), ckk
  logical               :: debug_dmrg_rassi_code = .false.
  real*8, allocatable   :: tmat(:,:) ! active-active rotation matrix


  if(.not.trorb)then
    write(lupri,'(a,a)') ' prepMPS: no MPS rotation requested for state ', &
    qcm_group_names(job)%states(ist)
    return
  else
    write(lupri,'(a,a)') ' prepMPS:    MPS rotation requested for state ', qcm_group_names(job)%states(ist)
  end if

  dmrg_orbital_space%nash(1:nsym) = nash(1:nsym)
  dmrg_symmetry%nirrep            = nsym
  dmrg_state%nactel               = nacte
  dmrg_state%ms2                  = mplet
  dmrg_state%irefsm               = lsym

  if(debug_dmrg_rassi_code)then
    write(lupri,*)' Entering prepMPS. TRA='
    write(lupri,'(1x,5f16.8)')(TRA(I),I=1,NTRA)
  end if

  !> find first the scaling factor to transform wrt the inactive orbitals
  fac(1,1) = 1.0d0
  ista     = 1
  do isym = 1, nsym
    no = nosh(isym)
    do i = 1, nish(isym)
      ii  = ista + (no+1)*(i-1)
      ckk = tra(ii)
      fac(1,1) = fac(1,1)*ckk
    end do
    ista = ista + no**2
  end do
  fac(1,1) = fac(1,1)**2
  if(debug_dmrg_rassi_code)then
    write(lupri,*) ' scaling factor for MPS (inactive orbital rotations)',fac
  end if

  dmrg_state%ms2                  = mspro

  if (nsym.gt.1) stop "MPS rotation not supported with symmetry"

  allocate(tmat(nash(1),nash(1)))
  tmat = 0.0d0

  ista = 1
  jorb = 0
  do isym = 1, nsym
    ni = nish(isym)
    no = nosh(isym)
    do i = 1, nash(isym)
      ! copy the active-active part of the rotation matrix into tmat
      tmat(i,:) = tra(1+(no+i)*ni+nash(isym)*(i-1):(no+i)*ni+nash(isym)*i)
    end do
    ista = ista + no**2
  end do

  ! rotate MPS
  call qcmaquis_mpssi_rotate(qcm_prefixes(job), &
                             ist,               &
                             tmat, &
                             nash(1)**2, &
                             fac(1,1), &
                             mspro)

  if (allocated(tmat)) deallocate(tmat)

  ! Avoid unused variable warnings
  if (.false.) then
    call unused_integer(istatereal)
  end if
#else
  write(lupri,*) ' calling prepMPS w/o DMRG interface - foolish!'
  write(lupri,*) ' ... no actual task is performed.            '
! Avoid unused variable warnings if DMRG is disabled
  if (.false.) then
    call unused_real_array(tra)
    call unused_integer(istate)
    call unused_integer(istatereal)
    call unused_integer(job)
    call unused_integer(ist)
    call unused_integer(lsym)
    call unused_integer(mplet)
    call unused_integer(mspro)
    call unused_integer(nacte)
    call unused_integer(nsym)
    call unused_integer(ntra)
    call unused_integer(lupri)
    call unused_integer_array(nish)
    call unused_integer_array(nash)
    call unused_integer_array(nosh)
    call unused_logical(trorb)
  end if
#endif
  end subroutine prepMPS
