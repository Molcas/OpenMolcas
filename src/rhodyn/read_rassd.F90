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
! Copyright (C) 2021,2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine read_rassd(nfile)
!***********************************************************************
! reads input RASSCF file RASSDX, where X = nfile
!***********************************************************************

use rhodyn_data, only: CI, DTOC, E, H_CSF, i_rasscf, ipglob, lroots, nconf, NDET, rassd_list, runmode
use mh5, only: mh5_close_file, mh5_exists_attr, mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset, mh5_open_file_r
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nfile
integer(kind=iwp) :: fileid
character(len=16) :: molcas_module
real(kind=wp), allocatable :: tmp(:,:)

fileid = mh5_open_file_r(rassd_list(nfile))

i_rasscf = 0
! if i_rasscf = 0 at the end of the file, that is sign of error
!               1 - spin-free H obtained (dset 'HCSF' is found)
!               2 - CIs and Es obtained (dsets 'CI_VECTORS' and 'ROOT_ENERGIES' in RASSCF)
!               3 - CIs and Es (dsets 'CI_VECTORS' and 'STATE_PT2_ENERGIES' from CASPT2)

! reading ci vectors and energies
if (mh5_exists_attr(fileid,'MOLCAS_MODULE')) then
  call mma_allocate(tmp,nconf(nfile),lroots(nfile),label='tmp')
  call mh5_fetch_attr(fileid,'MOLCAS_MODULE',molcas_module)
  if (molcas_module(1:6) == 'RASSCF') then
    ! rasscf input file:
    i_rasscf = 2
    if (ipglob > 2) write(u6,*) 'reading CI_VECTORS'
    call mh5_fetch_dset(fileid,'CI_VECTORS',tmp,[nconf(nfile),lroots(nfile)],[0,0])
    CI(1:nconf(nfile),1:lroots(nfile),nfile) = tmp(:,:)
    if (runmode /= 4) then
      if (ipglob > 2) write(u6,*) 'reading ROOT_ENERGIES'
      call mh5_fetch_dset(fileid,'ROOT_ENERGIES',E(1:lroots(nfile),nfile))
    end if
  else if (molcas_module(1:6) == 'CASPT2') then
    ! caspt2 input file:
    i_rasscf = 3
    !call mh5_fetch_dset(fileid,'H_EFF',H_CSF(1:nconf(nfile),1:nconf(nfile),nfile))
    if (mh5_exists_dset(fileid,'CI_VECTORS')) call mh5_fetch_dset(fileid,'CI_VECTORS',tmp)
    CI(1:nconf(nfile),1:lroots(nfile),nfile) = tmp(:,:)
    if (runmode /= 4) then
      if (mh5_exists_dset(fileid,'STATE_PT2_ENERGIES')) call mh5_fetch_dset(fileid,'STATE_PT2_ENERGIES',E(1:lroots(nfile),nfile))
    end if
  end if
  call mma_deallocate(tmp)
end if

! reading Hamiltonian (at the moment molcas does not write Hamiltonian in det basis to HDF5 file)
if (mh5_exists_dset(fileid,'HCSF')) then
  call mma_allocate(tmp,nconf(nfile),nconf(nfile),label='tmp')
  call mh5_fetch_dset(fileid,'HCSF',tmp)
  H_CSF(1:nconf(nfile),1:nconf(nfile),nfile) = tmp(:,:)
  call mma_deallocate(tmp)
  i_rasscf = 1
end if
if (mh5_exists_dset(fileid,'DTOC')) then
  call mma_allocate(tmp,NDET(nfile),nconf(nfile),label='tmp')
  call mh5_fetch_dset(fileid,'DTOC',tmp)
  DTOC(1:NDET(nfile),1:nconf(nfile),nfile) = tmp(:,:)
  call mma_deallocate(tmp)
end if

call mh5_close_file(fileid)

if (i_rasscf == 0) then
  write(u6,*) 'Error in reading RASSCF file ',rassd_list(nfile)
  write(u6,*) 'Required dsets has not been found in RASSD file'
  call abend()
end if

end subroutine read_rassd
