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

use rhodyn_data, only: CSF2SO, dipole, DM0, flag_pulse, HTOT_CSF, Nstate, prep_id, U_CI
use rhodyn_utils, only: dashes
use mh5, only: mh5_close_file, mh5_exists_dset, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, u6

implicit none
real(kind=wp), allocatable :: DIPI(:,:,:), DIPR(:,:,:), tmpi(:,:), tmpr(:,:)

if (mh5_is_hdf5('RDPREP')) then
  call dashes()
  write(u6,*) 'reading RDPREP file, opening'
  call dashes()
  prep_id = mh5_open_file_r('RDPREP')
  write(u6,*) 'RDPREP file opened'
else
  write(u6,*) 'RDPREP file is not of hdf5 format or not found'
  call abend()
end if

call mma_allocate(tmpr,Nstate,Nstate)
call mma_allocate(tmpi,Nstate,Nstate)

write(u6,*) 'hamiltonian extraction'
if (mh5_exists_dset(prep_id,'FULL_H_R') .and. mh5_exists_dset(prep_id,'FULL_H_I')) then
  call mh5_fetch_dset(prep_id,'FULL_H_R',tmpr)
  call mh5_fetch_dset(prep_id,'FULL_H_I',tmpi)
  HTOT_CSF(:,:) = cmplx(tmpr,tmpi,kind=wp)
else
  write(u6,*) 'RDPREP does not contain Hamiltonian'
  call abend()
end if

write(u6,*) 'dm extraction'
if (mh5_exists_dset(prep_id,'DM0_R') .and. mh5_exists_dset(prep_id,'DM0_I')) then
  call mh5_fetch_dset(prep_id,'DM0_R',tmpr)
  call mh5_fetch_dset(prep_id,'DM0_I',tmpi)
  DM0(:,:) = cmplx(tmpr,tmpi,kind=wp)
else
  write(u6,*) 'RDPREP does not contain density matrix'
  call abend()
end if

write(u6,*) 'csfso extraction'
if (mh5_exists_dset(prep_id,'CSF2SO_R') .and. mh5_exists_dset(prep_id,'CSF2SO_I')) then
  call mh5_fetch_dset(prep_id,'CSF2SO_R',tmpr)
  call mh5_fetch_dset(prep_id,'CSF2SO_I',tmpi)
  CSF2SO(:,:) = cmplx(tmpr,tmpi,kind=wp)
else
  write(u6,*) 'RDPREP does not contain CSF2SO matrix'
  write(u6,*) 'Propagation will be performed in the given basis!'
  flag_pulse = .false.
end if

write(u6,*) 'dipole real'
if (mh5_exists_dset(prep_id,'SOS_EDIPMOM_REAL') .and. mh5_exists_dset(prep_id,'SOS_EDIPMOM_IMAG')) then
  call mma_allocate(DIPR,Nstate,Nstate,3)
  call mma_allocate(DIPI,Nstate,Nstate,3)
  call mh5_fetch_dset(prep_id,'SOS_EDIPMOM_REAL',DIPR)
  call mh5_fetch_dset(prep_id,'SOS_EDIPMOM_IMAG',DIPI)
  dipole(:,:,:) = cmplx(DIPR,DIPI,kind=wp)
  call mma_deallocate(DIPR)
  call mma_deallocate(DIPI)
  flag_pulse = .true.
end if

write(u6,*) 'u_ci extraction'
if (mh5_exists_dset(prep_id,'U_CI')) then
  call mh5_fetch_dset(prep_id,'U_CI',U_CI)
end if

call mh5_close_file(prep_id)

call mma_deallocate(tmpr)
call mma_deallocate(tmpi)

call dashes()
write(u6,*) 'reading RDPREP has been finished'
call dashes()

end subroutine read_prep
