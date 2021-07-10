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
  use rhodyn_data
  use rhodyn_utils, only: dashes
  use definitions, only: wp, iwp, u6
  use mh5
  implicit none
  !? dimension is suspicious
  real(kind=wp),dimension(Nstate,lrootstot) :: SO_CI_R,  SO_CI_I
  real(kind=wp),dimension(Nstate,Nstate) :: tmpr, tmpi
  integer(kind=iwp):: fileid

!  call qEnter('read_rassisd')

  if (ipglob>3) then
    call dashes()
    write(u6,*) 'Reading RASSISD file'
    call dashes()
  endif

  fileid = mh5_open_file_r('RASSISD')

  if (mh5_exists_dset(fileid,'CH_SO_REAL')) then
    call mh5_fetch_dset(fileid,'CH_SO_REAL',tmpr)
  else
    write(u6,*) 'Error in reading RASSISD file, no CH_SO_REAL matrix'
    call abend()
  endif
  if (mh5_exists_dset(fileid,'CH_SO_IMAG')) then
    call mh5_fetch_dset(fileid,'CH_SO_IMAG',tmpi)
  else
    write(u6,*) 'Error in reading RASSISD file, no CH_SO_IMAG matrix'
    call abend()
  endif
  HSOCX = dcmplx(tmpr,tmpi)

  call mh5_fetch_dset(fileid,'SOCOEFF_REAL',SO_CI_R)
  call mh5_fetch_dset(fileid,'SOCOEFF_IMAG',SO_CI_I)
  SO_CI = dcmplx(SO_CI_R,SO_CI_I)

  do i=1,Nstate
    E_SO(i) = HSOCX(i,i)
  enddo

  call mh5_close_file(fileid)
!  call qExit('read_rassisd')
end
