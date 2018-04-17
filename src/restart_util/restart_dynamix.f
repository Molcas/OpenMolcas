************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8
#ifdef _HDF5_
      SUBROUTINE Restart_Dynamix(file_h5res)
      IMPLICIT NONE
#include "stdalloc.fh"
#include "mh5.fh"
      REAL*8, ALLOCATABLE ::    vel(:)
      REAL*8              ::    time,dt,e
      INTEGER             ::    attr_id,i,restart_fileid,natom,nsym
      Character           ::    tmp*256,sFile*128
      CHARACTER*180       ::    File_H5Res
      LOGICAL             ::    Exist

      write(6,'(A)') 'Restarting dynamix from h5 file', file_h5res

C Check the file exists
      call f_inquire(sFile,Exist)
      if(.not.Exist) then
        call getenvf('MOLCAS_SUBMIT_DIR',tmp)
        if(tmp.ne.' ') then
          i=index(tmp,' ')
          if(i.gt.0) then
            sFile=tmp(1:i-1)//'/'//file_h5res
            call f_inquire(sFile,Exist)
          endif
        endif
        if(.not.Exist) then
          Call WarningMessage(2,'File '//
     &                        sFile//' is not found')
          call Quit_OnUserError()
        endif
      endif

      restart_fileid = mh5_open_file_r(sFile)

C read number of atoms
      call mh5_fetch_attr(restart_fileid,'NSYM',nsym)
      if (nsym .gt. 1) then
        attr_id = mh5_open_attr(restart_fileid,'NATOMS_ALL')
      else
        attr_id = mh5_open_attr(restart_fileid,'NATOMS_UNIQUE')
      endif
      call mh5_get_attr(attr_id,natom)

C read time and save in RunFile
      call mh5_fetch_dset_scalar_real(restart_fileid,'TIME',time)
      CALL Put_dScalar('MD_Time',time)

C read time step and save in RunFile
      call mh5_fetch_dset_scalar_real(restart_fileid,'TIME_STEP',dt)
      CALL Put_dScalar('Timestep',dt)

C read max hop and save in RunFile
      if (mh5_exists_dset(restart_fileid,'MAX_HOP')) then
        call mh5_fetch_dset_scalar_int(restart_fileid,'MAX_HOP',i)
        CALL Put_iScalar('MaxHops',i)
      endif

C read total energy and save in RunFile
      call mh5_fetch_dset_scalar_real(restart_fileid,'ETOT',e)
      CALL Put_dScalar('MD_Etot',e)

C read velocities and save in RunFile
      CALL mma_allocate(vel,natom*3)
      call mh5_fetch_dset_array_real(restart_fileid,'VELOCITIES',vel)
      call Put_dArray('Velocities',vel,3*natom)
      CALL mma_deallocate(vel)

      call mh5_close_file(restart_fileid)

      RETURN
      END
#elif defined (NAGFOR)
c Some compilers do not like empty files
      SUBROUTINE empty_Restart_Dynamix()
      END
#endif
