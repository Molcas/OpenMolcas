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

subroutine Restart_Dynamix(File_H5res)

use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_open_attr, mh5_get_attr, mh5_close_attr, mh5_exists_dset, mh5_fetch_dset, &
               mh5_open_dset, mh5_get_dset_dims, mh5_close_dset, mh5_close_file
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
character(len=180), intent(in) :: File_H5Res
real(kind=wp), allocatable :: vel(:), NHC(:)
real(kind=wp) :: time, dt, e
integer(kind=iwp) :: attr_id, i, restart_fileid, natom, nsym, dset_id, nh(1)
character(len=256) :: tmp, sFile
logical(kind=iwp) :: Exists

write(u6,'(A)') 'Restarting dynamix from h5 file', file_h5res

! Check the file exists
sFile=File_H5Res
call f_inquire(sFile,Exists)
if (.not.Exists) then
  call getenvf('MOLCAS_SUBMIT_DIR',tmp)
  if (tmp /= ' ') then
    i=index(tmp,' ')
    if (i > 0) then
      sFile=tmp(1:i-1)//'/'//file_h5res
      call f_inquire(sFile,Exists)
    end if
  end if
  if (.not.Exists) then
    call WarningMessage(2,'File '//trim(sFile)//' is not found')
    call Quit_OnUserError()
  end if
end if

restart_fileid=mh5_open_file_r(sFile)

! read number of atoms
call mh5_fetch_attr(restart_fileid,'NSYM',nsym)
if (nsym > 1) then
  attr_id=mh5_open_attr(restart_fileid,'NATOMS_ALL')
else
  attr_id=mh5_open_attr(restart_fileid,'NATOMS_UNIQUE')
end if
call mh5_get_attr(attr_id,natom)
call mh5_close_attr(attr_id)

! read time and save in RunFile
call mh5_fetch_dset(restart_fileid,'TIME',time)
call Put_dScalar('MD_Time',time)

! read time step and save in RunFile
call mh5_fetch_dset(restart_fileid,'TIME_STEP',dt)
call Put_dScalar('Timestep',dt)

! read max hop and save in RunFile
if (mh5_exists_dset(restart_fileid,'MAX_HOP')) then
  call mh5_fetch_dset(restart_fileid,'MAX_HOP',i)
  call Put_iScalar('MaxHops',i)
end if

! read total energy and save in RunFile
call mh5_fetch_dset(restart_fileid,'ETOT',e)
call Put_dScalar('MD_Etot',e)

! read velocities and save in RunFile
call mma_allocate(vel,natom*3)
call mh5_fetch_dset(restart_fileid,'VELOCITIES',vel)
call Put_dArray('Velocities',vel,3*natom)
call mma_deallocate(vel)

! read thermostat info and save in RunFile
dset_id=mh5_open_dset(restart_fileid,'NOSEHOOVER')
call mh5_get_dset_dims(dset_id,nh)
call mh5_close_dset(dset_id)
call mma_allocate(NHC,nh(1))
call mh5_fetch_dset(restart_fileid,'NOSEHOOVER',NHC)
call Put_dArray('NOSEHOOVER',NHC,nh(1))
call mma_deallocate(NHC)

call mh5_close_file(restart_fileid)

return

end subroutine Restart_Dynamix

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(Restart_Dynamix)

#endif
