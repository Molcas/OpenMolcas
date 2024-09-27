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

#include "compiler_features.h"
#ifdef _HDF5_

subroutine read_hdf5_init(file_h5,nstate,nss)

use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_close_file
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
character(len=180), intent(in) :: file_h5
integer(kind=iwp), intent(out) :: nstate, nss
integer(kind=iwp) :: fileid, i
logical(kind=iwp) :: Exists
character(len=256) :: sFile, tmp
character(len=180) :: tmp2
character(len=5) :: molcas_module_kind
integer(kind=iwp), allocatable :: spin_mult(:)

write(u6,'(A,A)') 'Read data from rassi.h5 file ',trim(file_h5)
NSS = 0
NSTATE = 0
Exists = .false.
! Check if it is a rassi-h5 file:

! Check the file exists in $WorkDir
call f_inquire(trim(file_h5),Exists)
if (Exists) write(u6,*) 'file ',trim(file_h5),' exists!!!'
! if not present, look for it in the submit directory:
if (.not. Exists) then
  call getenvf('MOLCAS_SUBMIT_DIR',tmp)
  if (tmp /= ' ') then
    i = index(tmp,' ')
    if (i > 0) then
      sFile = trim(tmp(1:i-1)//'/'//file_h5)
      call f_inquire(sFile,Exists)
    end if
  end if
  ! if still not present, warn the user and abort:
  if (.not. Exists) then
    call WarningMessage(2,'File '//trim(file_h5)//' is not found')
    call Quit_OnUserError()
  end if
end if
!----------------------------------------------------------------------|
! open the file
fileid = mh5_open_file_r(trim(file_h5))
! check if it is an HDF5 file produced by RASSI
call mh5_fetch_attr(fileid,'MOLCAS_MODULE',tmp2)
molcas_module_kind = trim(tmp2)
#ifdef _DEBUGPRINT_
write(u6,'(A,A)') 'read_hdf5::  molcas_module=',molcas_module_kind
#endif
if (molcas_module_kind(1:5) /= 'RASSI') then
  call WarningMessage(2,'Input HDF5 file '//trim(file_h5)//' is not produced by RASSI')
  call Quit_OnUserError()
end if

!----------------------------------------------------------------------|
! read number of spin free states
call mh5_fetch_attr(fileid,'NSTATE',nstate)
#ifdef _DEBUGPRINT_
write(u6,'(A,I6)') 'read_hdf5::  NSTATE=',NSTATE
#endif
call Put_iScalar('NSTATE_SINGLE   ',NSTATE)

! read spin multiplicity of each state:
call mma_allocate(spin_mult,nstate,'nstate')
call mh5_fetch_attr(fileid,'STATE_SPINMULT',spin_mult)
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'spin_mult'
write(u6,'(20I4)') (spin_mult(i),i=1,nstate)
#endif
! compute the number of spin-orbit states:
nss = sum(spin_mult(:))
#ifdef _DEBUGPRINT_
write(u6,'(A,I6)') 'read_hdf5::     NSS=',NSS
#endif
call Put_iScalar('NSS_SINGLE      ',NSS)
call mma_deallocate(spin_mult)

!----------------------------------------------------------------------|
! close the file
call mh5_close_file(fileid)

return

end subroutine read_hdf5_init

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(read_hdf5_init)

#endif
