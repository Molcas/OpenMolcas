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
      SUBROUTINE Restart_SurfaceHop
      use Tully_variables
#ifdef _HDF5_
      use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_exists_attr,
     &               mh5_exists_dset, mh5_fetch_dset, mh5_close_file
#endif
      IMPLICIT NONE
#ifdef _HDF5_
#include "surfacehop.fh"
#include "stdalloc.fh"
      INTEGER             ::    i,nstates,nconfs,restart_fileid
      Character           ::    tmp*256,sFile*128
      LOGICAL             ::    Exist
      real*8              ::    dt
      real*8, allocatable ::    ener(:),ciarray(:),real_amatrix(:),
     &                          imag_amatrix(:)
      complex*16, allocatable :: amatrix(:)

      write(6,'(A)') 'Restarting surfacehop from h5 file', file_h5res

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

C Save on the RUNFILE information for DYNAMIX
      call restart_dynamix(file_h5res)

      call Get_dScalar('Timestep',DT)
      NSUBSTEPS=Int(200*DT/41.3)
      restart_fileid = mh5_open_file_r(sFile)
      call mh5_fetch_attr(restart_fileid,'NSTATES',nstates)
      call mh5_fetch_attr(restart_fileid,'NCONFS',nconfs)

C Save on the RUNFILE information for SURFACEHOP now

C read seed number and save in RunFile
      call mh5_fetch_dset_scalar_int(restart_fileid,'SEED',i)
      CALL Put_iScalar('Seed',i)

C read number of hops and save in RunFile
      if (mh5_exists_attr(restart_fileid,'NO. OF HOPS')) then
        call mh5_fetch_dset(restart_fileid,'NO. OF HOPS',i)
        CALL Put_iScalar('Number of Hops',i)
      endif

C read max hop for Tully and save in RunFile
      if (mh5_exists_dset(restart_fileid,'MAX_HOP_TULLY')) then
        call mh5_fetch_dset(restart_fileid,'MAX_HOP_TULLY',i)
        CALL Put_iScalar('MaxHopsTully',i)
      endif

C read relax root number and save in RunFile
      call mh5_fetch_dset(restart_fileid,'Relax CAS root',i)
      CALL Put_iScalar('Relax CASSCF root',i)

C read the energies of the previous step and save in RunFile
      CALL mma_allocate(ener,nstates)
      call mh5_fetch_dset(restart_fileid,'Energ Prev',ener)
      CALL Put_darray('VenergyP',ener,nstates)
      CALL mma_deallocate(ener)

C read the CI arrays of the previous step and save in RunFile
      CALL mma_allocate(ciarray,nstates*nconfs)
      call mh5_fetch_dset(restart_fileid,'CI Prev',ciarray)
      CALL Put_darray('AllCIP',ciarray,nstates*nconfs)
      CALL mma_deallocate(ciarray)

C read the CI arrays of the step before the previous step and save in RunFile
      CALL mma_allocate(ciarray,nstates*nconfs)
      call mh5_fetch_dset(restart_fileid,'CI PPrev',ciarray)
      CALL Put_darray('AllCIPP',ciarray,nstates*nconfs)
      CALL mma_deallocate(ciarray)

C read the AmatrixV and save in RunFile
      CALL mma_allocate(real_amatrix,nstates*nstates)
      CALL mma_allocate(imag_amatrix,nstates*nstates)
      CALL mma_allocate(amatrix,nstates*nstates)
      call mh5_fetch_dset(restart_fileid,'AmatrixV-R',real_amatrix)
      call mh5_fetch_dset(restart_fileid,'AmatrixV-I',imag_amatrix)
      amatrix(:) = DCMPLX(real_amatrix,imag_amatrix)
      CALL Put_zarray('AmatrixV',amatrix,nstates*nstates)
      CALL mma_deallocate(amatrix)
      CALL mma_deallocate(real_amatrix)
      CALL mma_deallocate(imag_amatrix)

      call mh5_close_file(restart_fileid)

#endif
      RETURN
      END
