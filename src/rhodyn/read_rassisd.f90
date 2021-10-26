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
subroutine read_rassisd
!***********************************************************************
!
! Purpose :  In case of charge migration read in the complex
!            hamiltonian from the MOLCAS output rassisd file (SO)
!
!***********************************************************************
  use rhodyn_data, only: Nstate, lrootstot, ipglob, HSOCX, E_SO, i, &
                         SO_CI, preparation, V_SO, dipole, flag_so, &
                         flag_dyson, dysamp, E_SF
  use rhodyn_utils, only: dashes
  use definitions, only: wp, iwp, u6
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5, only: mh5_open_file_r, mh5_exists_dset, mh5_fetch_dset, &
                 mh5_close_file
  implicit none
  real(kind=wp),dimension(:), allocatable :: tmpe
  real(kind=wp),dimension(:,:),allocatable :: tmpr, tmpi
  real(kind=wp),allocatable,dimension(:,:,:) :: DIPR, DIPI
  integer(kind=iwp):: fileid

  if (ipglob>3) then
    call dashes()
    write(u6,*) 'Reading RASSI file'
    call dashes()
  endif

  fileid = mh5_open_file_r('RASSISD')

  if (preparation==4) then
    if (flag_so) then
    ! reading complex Hamiltonian (already with SOC included)
    ! so far needed only for charge migration case
    ! take care of dimensions!
    call mma_allocate(tmpr,Nstate,Nstate)
    call mma_allocate(tmpi,Nstate,Nstate)
    if (mh5_exists_dset(fileid,'CH_SO_REAL').and. &
        mh5_exists_dset(fileid,'CH_SO_IMAG')) then
      call mh5_fetch_dset(fileid,'CH_SO_REAL',tmpr)
      call mh5_fetch_dset(fileid,'CH_SO_IMAG',tmpi)
    else if &! if using standard rassi dsets in rassi.h5
        (mh5_exists_dset(fileid,'HSO_MATRIX_REAL').and.&
         mh5_exists_dset(fileid,'HSO_MATRIX_IMAG')) then
      call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',tmpr)
      call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',tmpi)
    else
      write(u6,*) 'Error in reading RASSI file, no CH_SO_REAL matrix,'
      write(u6,*) 'nor HSO_MATRIX_REAL/IMAG datasets'
      call abend()
    endif
    HSOCX = dcmplx(tmpr,tmpi)
    if (allocated (tmpr)) call mma_deallocate(tmpr)
    if (allocated (tmpi)) call mma_deallocate(tmpi)
    ! filling in energies
    do i=1,Nstate
      E_SO(i) = HSOCX(i,i)
    enddo
    else ! if flag_so is off
      write(u6,*) 'Reading SF energies SFS_ENERGIES and construct H'
      if (mh5_exists_dset(fileid,'SFS_ENERGIES')) then
        call mma_allocate(tmpe,Nstate)
        call mh5_fetch_dset(fileid,'SFS_ENERGIES',tmpe)
        E_SF = tmpe
        call mma_deallocate(tmpe)
      else
        write(u6,*) 'Error in reading RASSI file, no SFS_ENERGIES'
        call abend()
      endif
      do i=1,Nstate
        HSOCX(i,i) = E_SF(i)
      enddo
    endif
  endif

! reading pure SOC matrix in SF basis
  if (preparation/=4.and.flag_so) then
    call mma_allocate(tmpr,lrootstot,lrootstot)
    call mma_allocate(tmpi,lrootstot,lrootstot)
    if (mh5_exists_dset(fileid,'V_SO_REAL').and. &
        mh5_exists_dset(fileid,'V_SO_IMAG')) then
      call mh5_fetch_dset(fileid,'V_SO_REAL',tmpr)
      call mh5_fetch_dset(fileid,'V_SO_IMAG',tmpi)
    else
      write(u6,*) 'Error in reading RASSISD file, no V_SO matrix'
      call abend()
    endif
    V_SO = dcmplx(tmpr,tmpi)
    if (allocated (tmpr)) call mma_deallocate(tmpr)
    if (allocated (tmpi)) call mma_deallocate(tmpi)
  endif

! reading SOC coefficients
  if (flag_so) then
  call mma_allocate(tmpr,lrootstot,lrootstot)
  call mma_allocate(tmpi,lrootstot,lrootstot)
  if (mh5_exists_dset(fileid,'SOCOEFF_REAL').and. &
      mh5_exists_dset(fileid,'SOCOEFF_IMAG')) then
    call mh5_fetch_dset(fileid,'SOCOEFF_REAL',tmpr)
    call mh5_fetch_dset(fileid,'SOCOEFF_IMAG',tmpi)
  else if &! if using standard rassi dsets of rassi.h5
      (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL').and.&
       mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',tmpr)
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',tmpi)
  else
    write(u6,*) 'Error in reading RASSI file, no SOCOEFF matrix'
    call abend()
  endif
  SO_CI = dcmplx(tmpr,tmpi)
  if (allocated (tmpr)) call mma_deallocate(tmpr)
  if (allocated (tmpi)) call mma_deallocate(tmpi)
  endif

! matrices of dipole moments
  call mma_allocate(DIPR,lrootstot,lrootstot,3)
  call mma_allocate(DIPI,lrootstot,lrootstot,3)
  if (flag_so) then
    if (mh5_exists_dset(fileid,'SOS_EDIPMOM_REAL').and. &
      mh5_exists_dset(fileid,'SOS_EDIPMOM_IMAG')) then
      call mh5_fetch_dset(fileid,'SOS_EDIPMOM_REAL',DIPR)
      call mh5_fetch_dset(fileid,'SOS_EDIPMOM_IMAG',DIPI)
    else
      write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SO basis'
      call abend()
    endif
  else ! to read SFS_EDIPMOM (flag_so = off)
    if (mh5_exists_dset(fileid,'SFS_EDIPMOM')) then
      call mh5_fetch_dset(fileid,'SFS_EDIPMOM',DIPR)
      DIPI=0d0
    else
      write(u6,*) 'Error in reading RASSISD file, no dipole matrix in SF basis'
      call abend()
    endif
  endif
  !write(u6,*) 'dysorb read'
  if (mh5_exists_dset(fileid,'DYSORB').and.flag_dyson) then
    call mh5_fetch_dset(fileid,'DYSORB',dysamp)
  else if (mh5_exists_dset(fileid,'DYSAMP').and.flag_dyson) then
    call mh5_fetch_dset(fileid,'DYSAMP',dysamp)
  else
    if (ipglob>2) then
      write(u6,*)'Ionization is not taken into account (set flag DYSO)'
      write(u6,*)' and/or RASSI file does not contain Dyson amplitudes'
      flag_dyson=.False.
    endif
  endif
  !write(u6,*) 'dysorb has been read'
  dipole = dcmplx(DIPR,DIPI)
  if (allocated(DIPR)) call mma_deallocate(DIPR)
  if (allocated(DIPI)) call mma_deallocate(DIPI)

  call mh5_close_file(fileid)

end
