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

module construct_lpdft_ham

implicit none
private

public :: construct_lpdft_ham_from_tdm

contains

subroutine construct_lpdft_ham_from_tdm()

  use rasscf_global, only: lRoots, NAC, NACPR2
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp, u6
#ifdef _HDF5_
  use mh5, only: mh5_close_file, mh5_exists_dset, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
#endif

  implicit none
#include "warnings.h"
  integer(kind=iwp) :: nPairs
  real(kind=wp), allocatable :: TDM1(:,:,:), TDM2(:,:)
#ifdef _HDF5_
  character(len=1024) :: currDir, fileName, project
  integer(kind=iwp) :: mh5id
  integer(kind=iwp), allocatable :: rootPairs(:,:)
#endif

  if (lRoots <= 1) then
    call WarningMessage(1,'LPDF MS-PDFT requested, but only one root is available; no transition densities will be read.')
    return
  end if

  nPairs = (lRoots*(lRoots-1))/2

#ifndef _HDF5_
  call WarningMessage(2,'LPDF MS-PDFT requested, but HDF5 support is disabled.')
  call Quit_OnUserError()
#else
  call get_environment_variable('Project',project)
  if (len_trim(project) == 0) project = 'Noname'
  call get_environment_variable('CurrDir',currDir)
  if (len_trim(currDir) > 0) then
    fileName = trim(currDir)//'/'//trim(project)//'_tdm_4_lpdft.h5'
  else
    fileName = trim(project)//'_tdm_4_lpdft.h5'
  end if

  if (.not. mh5_is_hdf5(trim(fileName))) then
    call WarningMessage(2,'LPDF transition density HDF5 file not found: '//trim(fileName))
    call Quit_OnUserError()
  end if

  mh5id = mh5_open_file_r(trim(fileName))
  if (.not. mh5_exists_dset(mh5id,'ACTIVE_1TDM')) then
    call WarningMessage(2,'LPDF transition density HDF5 file does not contain ACTIVE_1TDM.')
    call Quit_OnUserError()
  end if
  if (.not. mh5_exists_dset(mh5id,'ACTIVE_2TDM')) then
    call WarningMessage(2,'LPDF transition density HDF5 file does not contain ACTIVE_2TDM.')
    call Quit_OnUserError()
  end if

  call mma_allocate(TDM1,NAC,NAC,nPairs,Label='TDM1')
  call mma_allocate(TDM2,NACPR2,nPairs,Label='TDM2')
  call mh5_fetch_dset(mh5id,'ACTIVE_1TDM',TDM1)
  call mh5_fetch_dset(mh5id,'ACTIVE_2TDM',TDM2)

  if (mh5_exists_dset(mh5id,'ROOT_PAIRS')) then
    call mma_allocate(rootPairs,2,nPairs,Label='rootPairs')
    call mh5_fetch_dset(mh5id,'ROOT_PAIRS',rootPairs)
    call mma_deallocate(rootPairs)
  end if

  call mh5_close_file(mh5id)

  write(u6,'(6X,A)') 'LPDF: read transition density matrices from '//trim(fileName)
  write(u6,'(6X,A,I6,A)') 'LPDF: loaded transition density matrices for ',nPairs,' root pairs.'

  call mma_deallocate(TDM1)
  call mma_deallocate(TDM2)
#endif

end subroutine construct_lpdft_ham_from_tdm

end module construct_lpdft_ham
