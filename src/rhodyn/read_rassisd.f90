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
                         SO_CI
  use rhodyn_utils, only: dashes
  use definitions, only: wp, iwp, u6
  use mh5, only: mh5_open_file_r, mh5_exists_dset, mh5_fetch_dset, &
                 mh5_close_file
  implicit none
  !? dimension is suspicious
  real(kind=wp),dimension(Nstate,lrootstot) :: SO_CI_R,  SO_CI_I
  real(kind=wp),dimension(Nstate,Nstate) :: tmpr, tmpi
  integer(kind=iwp):: fileid

  if (ipglob>3) then
    call dashes()
    write(u6,*) 'Reading RASSISD file'
    call dashes()
  endif

  fileid = mh5_open_file_r('RASSISD')

  if (mh5_exists_dset(fileid,'CH_SO_REAL').and. &
      mh5_exists_dset(fileid,'CH_SO_IMAG')) then
    call mh5_fetch_dset(fileid,'CH_SO_REAL',tmpr)
    call mh5_fetch_dset(fileid,'CH_SO_IMAG',tmpi)
  else if &! if using standard rassi dsets of rassi.h5
      (mh5_exists_dset(fileid,'HSO_MATRIX_REAL').and.&
       mh5_exists_dset(fileid,'HSO_MATRIX_IMAG')) then
    call mh5_fetch_dset(fileid,'HSO_MATRIX_REAL',tmpr)
    call mh5_fetch_dset(fileid,'HSO_MATRIX_IMAG',tmpi)
  else
    write(u6,*) 'Error in reading RASSISD file, no CH_SO_REAL matrix,'
    write(u6,*) 'nor HSO_MATRIX_REAL/IMAG datasets'
    call abend()
  endif
  HSOCX = dcmplx(tmpr,tmpi)

! reading SOC coefficients
  if (mh5_exists_dset(fileid,'SOCOEFF_REAL').and. &
      mh5_exists_dset(fileid,'SOCOEFF_IMAG')) then
    call mh5_fetch_dset(fileid,'SOCOEFF_REAL',SO_CI_R)
    call mh5_fetch_dset(fileid,'SOCOEFF_IMAG',SO_CI_I)
  else if &! if using standard rassi dsets of rassi.h5
      (mh5_exists_dset(fileid,'SOS_COEFFICIENTS_REAL').and.&
       mh5_exists_dset(fileid,'SOS_COEFFICIENTS_IMAG')) then
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_REAL',SO_CI_R)
    call mh5_fetch_dset(fileid,'SOS_COEFFICIENTS_IMAG',SO_CI_I)
  else
    write(u6,*) 'Error in reading RASSISD file, no SOCOEFF matrix'
    call abend()
  endif

  SO_CI = dcmplx(SO_CI_R,SO_CI_I)

  do i=1,Nstate
    E_SO(i) = HSOCX(i,i)
  enddo

  call mh5_close_file(fileid)

end
