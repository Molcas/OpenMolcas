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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine read_prep()

  use rhodyn_data
  use rhodyn_utils, only: dashes
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5
  implicit none
  real(8),allocatable,dimension(:,:,:) :: DIPR, DIPI
  real(8),dimension(:,:),allocatable :: tmpr, tmpi

  if (mh5_is_hdf5('RDPREP')) then
    call dashes()
    write(*,*) 'reading RDPREP file, opening'
    call dashes()
    prep_id=mh5_open_file_r('RDPREP')
    write(*,*) 'RDPREP file opened'
  else
    write(*,*) 'RDPREP file is not of hdf5 format or not found'
    call abend()
  endif

  call mma_allocate(tmpr,Nstate,Nstate)
  call mma_allocate(tmpi,Nstate,Nstate)

  write(*,*)'hamiltonian extraction'
  if (mh5_exists_dset(prep_id,'FULL_H_R').and. &
      mh5_exists_dset(prep_id,'FULL_H_I')) then
    call mh5_fetch_dset_array_real(prep_id,'FULL_H_R',tmpr)
    call mh5_fetch_dset_array_real(prep_id,'FULL_H_I',tmpi)
    HTOT_CSF=dcmplx(tmpr,tmpi)
  else
    write(*,*) 'RDPREP does not contain Hamiltonian'
    call abend()
  endif

  write(*,*)'dm extraction'
  if (mh5_exists_dset(prep_id,'DM0_R').and. &
        mh5_exists_dset(prep_id,'DM0_I')) then
    call mh5_fetch_dset_array_real(prep_id, 'DM0_R',tmpr)
    call mh5_fetch_dset_array_real(prep_id,'DM0_I',tmpi)
    DM0=dcmplx(tmpr,tmpi)
  else
    write(*,*) 'RDPREP does not contain density matrix'
    call abend()
  endif

  write(*,*)'csfso extraction'
  if (mh5_exists_dset(prep_id,'CSF2SO_R').and. &
        mh5_exists_dset(prep_id,'CSF2SO_I')) then
    call mh5_fetch_dset_array_real(prep_id, 'CSF2SO_R',tmpr)
    call mh5_fetch_dset_array_real(prep_id,'CSF2SO_I',tmpi)
    CSF2SO=dcmplx(tmpr,tmpi)
  else
    write(*,*) 'RDPREP does not contain CSF2SO matrix'
    write(*,*) 'Propagation will be performed in the given basis!'
    flag_test=.True.
    flag_pulse=.False.
  endif
	  
	write(*,*) 'dipole real'
  if (mh5_exists_dset(prep_id,'SOS_EDIPMOM_REAL').and. &
      mh5_exists_dset(prep_id,'SOS_EDIPMOM_IMAG')) then
	  call mma_allocate (DIPR,Nstate,Nstate,3)
    call mma_allocate (DIPI,Nstate,Nstate,3)
    call mh5_fetch_dset_array_real(prep_id,'SOS_EDIPMOM_REAL',DIPR)
    call mh5_fetch_dset_array_real(prep_id,'SOS_EDIPMOM_IMAG',DIPI)
	  dipole = dcmplx(DIPR,DIPI)
	  call mma_deallocate (DIPR)
    call mma_deallocate (DIPI)
	  flag_pulse=.True.
  endif

  write(*,*)'u_ci extraction'
  if (mh5_exists_dset(prep_id,'U_CI')) then
    call mh5_fetch_dset_array_real(prep_id,'U_CI',U_CI)
  endif

  call mh5_close_file(prep_id)

  call mma_deallocate(tmpr)
  call mma_deallocate(tmpi)

  call dashes()
  write(*,*) 'reading RDPREP has been finished'
  call dashes()

end
