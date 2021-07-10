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
subroutine read_rassd(nfile)

  use rhodyn_data
  use definitions, only: iwp, u6
  use mh5
  implicit none
  integer(kind=iwp),intent(in):: nfile
  integer(kind=iwp) :: fileid
  character(len=16) :: molcas_module

  write(u6,*) 'read_rassd entered'
!  if (ipglob>3) then
!  call dashes(72)
!  write(u6,*) 'Reading file: ', rassd_list(nfile)
!  endif

  fileid = mh5_open_file_r(rassd_list(nfile))

  i_rasscf = 0

! reading ci vectors and energies
  if (mh5_exists_attr(fileid,'MOLCAS_MODULE')) then
    call mh5_fetch_attr(fileid,'MOLCAS_MODULE', molcas_module)
! rasscf input file:
    if (trim(molcas_module(1:6))=='RASSCF') then
      i_rasscf=2
      write(u6,*) 'reading CI_VECTORS'
      call mh5_fetch_dset(fileid,'CI_VECTORS', &
                  CI(1:nconf(nfile),1:lroots(nfile),nfile), &
                  [nconf(nfile),lroots(nfile)],[0,0])
      write(u6,*) 'reading ROOT_ENERGIES'
      call mh5_fetch_dset(fileid,'ROOT_ENERGIES', &
                  E(1:lroots(nfile),nfile))
      write(u6,*) 'reading finished'
! caspt2 input file:
    elseif (trim(molcas_module(1:6))=='CASPT2') then
      i_rasscf=3
!       call mh5_fetch_dset(fileid,'H_EFF',
!     &              H_CSF(1:nconf(nfile),1:nconf(nfile),nfile))
      call mh5_fetch_dset(fileid,'CI_VECTORS', &
                  CI(1:nconf(nfile),1:lroots(nfile),nfile))
      call mh5_fetch_dset(fileid,'STATE_PT2_ENERGIES', &
                  E(1:lroots(nfile),nfile))
    endif
  endif

! reading Hamiltonian (at the moment molcas doesnot write rasscf
! Hamiltonian to HDF5 file)
  if (mh5_exists_dset(fileid,'HCSF')) then
    call mh5_fetch_dset(fileid,'HCSF', &
                  H_CSF(1:nconf(nfile),1:nconf(nfile),nfile))
    i_rasscf = 1
  endif
  if (mh5_exists_dset(fileid,'DTOC')) then
    call mh5_fetch_dset(fileid,'DTOC', &
                  DTOC(1:NDET(nfile),1:nconf(nfile),nfile))
  endif

  if (i_rasscf==0) then
    write(u6,*) 'Error in reading RASSCF outputs'
    call abend()
  endif

  call mh5_close_file(fileid)

end
