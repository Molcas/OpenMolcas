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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine read_rassisd()
!***********************************************************************
! Purpose :  In case of charge migration read in the complex
!            hamiltonian from the MOLCAS output rassisd file (SO)
!***********************************************************************

use rhodyn_data, only: dipole, dysamp, E_SF, n_sf, E_SO, flag_dyson, flag_so, HSOCX, ipglob, lrootstot, Nstate, &
                       runmode, SO_CI, V_SO, basis
use mh5, only: mh5_close_file, mh5_exists_dset, mh5_fetch_dset, mh5_open_file_r
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, fileid
real(kind=wp), allocatable :: DIPI(:,:,:), DIPR(:,:,:), tmpe(:), tmpi(:,:), tmpr(:,:)

call StatusLine('RhoDyn:','Read RASSI H5 file')

fileid = mh5_open_file_r('RASSISD')

if (runmode == 4) then
  if (flag_so) then
    ! reading complex Hamiltonian (already with SOC included)
    ! so far needed only for charge migration case
    ! take care of dimensions!
    call mma_allocate(tmpr,Nstate,Nstate,label='tmpr')
    call mma_allocate(tmpi,Nstate,Nstate,label='tmpi')
    if (mh5_exists_dset(fileid,'CH_SO_REAL') .and. mh5_exists_dset(fileid,'CH_SO_IMAG')) then
      call mh5_fetch_dset(fileid,'CH_SO_REAL',tmpr)
      call mh5_fetch_dset(fileid,'CH_SO_IMAG',tmpi)
    else if (mh5_exists_dset(fileid,'HSO_MATRIX_REAL') .and. mh5_exists_dset(fileid,'HSO_MATRIX_IMAG')) then
      ! if using standard rassi dsets in rassi.h5
      call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',tmpr)
      call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',tmpi)
    else
      write(u6,*) 'Error in reading RASSI file, no CH_SO_REAL matrix,'
      write(u6,*) 'nor HSO_MATRIX_REAL/IMAG datasets'
      call abend()
    end if
    HSOCX(:,:) = cmplx(tmpr,tmpi,kind=wp)
    if (allocated(tmpr)) call mma_deallocate(tmpr)
    if (allocated(tmpi)) call mma_deallocate(tmpi)
    ! filling in energies
    do i=1,Nstate
      E_SO(i) = HSOCX(i,i)
    end do
  else ! if flag_so is off
    write(u6,*) 'Reading SF energies SFS_ENERGIES and construct H'
    if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
      call mma_allocate(tmpe,Nstate)
      call mh5_fetch_dset(fileid,'SFS_ENERGIES',tmpe)
      E_SF(:) = tmpe
      call mma_deallocate(tmpe)
    else
      write(u6,*) 'Error in reading RASSI file, no SFS_ENERGIES'
      call abend()
    end if
    do i=1,Nstate
      HSOCX(i,i) = E_SF(i)
    end do
  end if
end if

! reading pure SOC matrix in SF basis
if ((runmode /= 4) .and. flag_so) then
  call mma_allocate(tmpr,lrootstot,lrootstot)
  call mma_allocate(tmpi,lrootstot,lrootstot)
  if (mh5_exists_dset(fileid,'V_SO_REAL') .and. mh5_exists_dset(fileid,'V_SO_IMAG')) then
    call mh5_fetch_dset(fileid,'V_SO_REAL',tmpr)
    call mh5_fetch_dset(fileid,'V_SO_IMAG',tmpi)
  else
    write(u6,*) 'Error in reading RASSISD file, no V_SO matrix'
    call abend()
  end if
  V_SO(:,:) = cmplx(tmpr,tmpi,kind=wp)
  if (allocated(tmpr)) call mma_deallocate(tmpr)
  if (allocated(tmpi)) call mma_deallocate(tmpi)
end if

! reading SOC coefficients
if (flag_so) then
  call mma_allocate(tmpr,lrootstot,lrootstot,label='tmpr')
  call mma_allocate(tmpi,lrootstot,lrootstot,label='tmpi')
  if (mh5_exists_dset(fileid,'SOCOEFF_REAL') .and. mh5_exists_dset(fileid,'SOCOEFF_IMAG')) then
    call mh5_fetch_dset(fileid,'SOCOEFF_REAL',tmpr)
    call mh5_fetch_dset(fileid,'SOCOEFF_IMAG',tmpi)
  else if (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL') .and. mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
    ! if using standard rassi dsets of rassi.h5
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',tmpr)
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',tmpi)
  else
    write(u6,*) 'Error in reading RASSI file, no SOCOEFF matrix'
    call abend()
  end if
  SO_CI(:,:) = cmplx(tmpr,tmpi,kind=wp)
  if (allocated(tmpr)) call mma_deallocate(tmpr)
  if (allocated(tmpi)) call mma_deallocate(tmpi)
end if

! matrices of dipole moments
call mma_allocate(DIPR,lrootstot,lrootstot,3,label='DIPR')
call mma_allocate(DIPI,lrootstot,lrootstot,3,label='DIPI')
if (flag_so) then
  if (mh5_exists_dset(fileid,'SOS_EDIPMOM_REAL') .and. mh5_exists_dset(fileid,'SOS_EDIPMOM_IMAG')) then
    call mh5_fetch_dset(fileid,'SOS_EDIPMOM_REAL',DIPR)
    call mh5_fetch_dset(fileid,'SOS_EDIPMOM_IMAG',DIPI)
  else
    write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SO basis'
    call abend()
  end if
else ! to read SFS_EDIPMOM (flag_so = off)
  if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
    call mh5_fetch_dset(fileid,'SFS_EDIPMOM',DIPR)
    DIPI = Zero
  else
    write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SF basis'
    call abend()
  end if
end if
!write(u6,*) 'dysorb read'
if (mh5_exists_dset(fileid,'DYSORB') .and. flag_dyson) then
  call mh5_fetch_dset(fileid,'DYSORB',dysamp)
else if (mh5_exists_dset(fileid,'DYSAMP') .and. flag_dyson) then
  call mh5_fetch_dset(fileid,'DYSAMP',dysamp)
else
  if (ipglob > 2) then
    write(u6,*) 'Ionization is not taken into account (set flag DYSO)'
    write(u6,*) ' and/or RASSI file does not contain Dyson amplitudes'
    flag_dyson = .false.
  end if
end if
!write(u6,*) 'dysorb has been read'

if (basis == 'SPH') then
  if (ipglob > 2) write(u6,*) 'Reading SF energies SFS_ENERGIES'
  if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
    call mma_allocate(tmpe,n_sf)
    call mh5_fetch_dset(fileid,'SFS_ENERGIES',tmpe)
    E_SF(:) = tmpe
    call mma_deallocate(tmpe)
  else
    write(u6,*) 'Error in reading RASSI file, no SFS_ENERGIES'
    call abend()
  end if
  ! sf dipole needed
  call mma_deallocate(dipole)
  call mma_deallocate(DIPR)
  call mma_deallocate(DIPI)
  call mma_allocate(dipole,n_sf,n_sf,3)
  call mma_allocate(DIPR,n_sf,n_sf,3)
  call mma_allocate(DIPI,n_sf,n_sf,3)
  if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
    call mh5_fetch_dset(fileid,'SFS_EDIPMOM',DIPR)
    DIPI = Zero
  else
    write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SF basis'
    call abend()
  end if
end if

dipole(:,:,:) = cmplx(DIPR,DIPI,kind=wp)
if (allocated(DIPR)) call mma_deallocate(DIPR)
if (allocated(DIPI)) call mma_deallocate(DIPI)
call mh5_close_file(fileid)

end subroutine read_rassisd
