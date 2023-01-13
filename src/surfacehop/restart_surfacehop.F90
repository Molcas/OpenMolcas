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

subroutine Restart_SurfaceHop

#ifdef _HDF5_
use Tully_variables, only: NSUBSTEPS
use stdalloc, only: mma_allocate, mma_deallocate
use mh5, only: mh5_open_file_r, mh5_fetch_attr, mh5_exists_attr, mh5_exists_dset, mh5_fetch_dset, mh5_close_file
use surfacehop_globals, only: File_H5Res
use Constants, only: auTofs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, nstates, nconfs, restart_fileid
character(len=256) :: tmp
character(len=128) :: sFile
logical(kind=iwp) :: Exists
real(kind=wp) :: dt
real(kind=wp), allocatable :: ener(:), ciarray(:), real_amatrix(:), imag_amatrix(:), overlap_save(:), oldphase(:)
complex(kind=wp), allocatable :: amatrix(:)

write(u6,'(A)') 'Restarting surfacehop from h5 file',file_h5res

! Check the file exists
call f_inquire(sFile,Exists)
if (.not. Exists) then
  call getenvf('MOLCAS_SUBMIT_DIR',tmp)
  if (tmp /= ' ') then
    i = index(tmp,' ')
    if (i > 0) then
      sFile = tmp(1:i-1)//'/'//File_H5Res
      call f_inquire(sFile,Exists)
    end if
  end if
  if (.not. Exists) then
    call WarningMessage(2,'File '//trim(sFile)//' is not found')
    call Quit_OnUserError()
  end if
end if

! Save on the RUNFILE information for DYNAMIX
call restart_dynamix(File_H5Res)

call Get_dScalar('Timestep',DT)
NSUBSTEPS = int(200*DT*auTofs)
restart_fileid = mh5_open_file_r(sFile)
call mh5_fetch_attr(restart_fileid,'NSTATES',nstates)
call mh5_fetch_attr(restart_fileid,'NCONFS',nconfs)

! Save on the RUNFILE information for SURFACEHOP now

! read seed number and save in RunFile
call mh5_fetch_dset(restart_fileid,'SEED',i)
call Put_iScalar('Seed',i)

! read number of hops and save in RunFile
if (mh5_exists_attr(restart_fileid,'NO. OF HOPS')) then
  call mh5_fetch_dset(restart_fileid,'NO. OF HOPS',i)
  call Put_iScalar('Number of Hops',i)
end if

! read max hop for Tully and save in RunFile
if (mh5_exists_dset(restart_fileid,'MAX_HOP_TULLY')) then
  call mh5_fetch_dset(restart_fileid,'MAX_HOP_TULLY',i)
  call Put_iScalar('MaxHopsTully',i)
end if

! read relax root number and save in RunFile
call mh5_fetch_dset(restart_fileid,'RELAX CAS ROOT',i)
call Put_iScalar('Relax CASSCF root',i)

! read the energies of the previous step and save in RunFile
call mma_allocate(ener,nstates)
call mh5_fetch_dset(restart_fileid,'ENERG PREV',ener)
call Put_darray('VenergyP',ener,nstates)
call mma_deallocate(ener)

! read the CI arrays of the previous step and save in RunFile
call mma_allocate(ciarray,nstates*nconfs)
call mh5_fetch_dset(restart_fileid,'CI PREV',ciarray)
call Put_darray('AllCIP',ciarray,nstates*nconfs)
call mma_deallocate(ciarray)

! read the CI arrays of the step before the previous step and save in RunFile
call mma_allocate(ciarray,nstates*nconfs)
call mh5_fetch_dset(restart_fileid,'CI PPREV',ciarray)
call Put_darray('AllCIPP',ciarray,nstates*nconfs)
call mma_deallocate(ciarray)

! read <t-2dt|t-dt> overlap and associated phase if exists and save in RunFile
if (mh5_exists_dset(restart_fileid,'RASSI_SAVE_OVLP')) then
  call mma_allocate(overlap_save,nstates*nstates)
  call mma_allocate(oldphase,nstates)
  call mh5_fetch_dset(restart_fileid,'RASSI_SAVE_OVLP',overlap_save)
  call mh5_fetch_dset(restart_fileid,'OLD_OVLP_PHASE',oldphase)
  call Put_darray('SH_Ovlp_Save',overlap_save,nstates*nstates)
  call Put_darray('Old_Phase',oldphase,nstates)
  call mma_deallocate(overlap_save)
  call mma_deallocate(oldphase)
end if

! read the AmatrixV and save in RunFile
call mma_allocate(real_amatrix,nstates*nstates)
call mma_allocate(imag_amatrix,nstates*nstates)
call mma_allocate(amatrix,nstates*nstates)
call mh5_fetch_dset(restart_fileid,'AMATRIXV-R',real_amatrix)
call mh5_fetch_dset(restart_fileid,'AMATRIXV-I',imag_amatrix)
amatrix(:) = cmplx(real_amatrix,imag_amatrix,kind=wp)
call Put_zarray('AmatrixV',amatrix,nstates*nstates)
call mma_deallocate(amatrix)
call mma_deallocate(real_amatrix)
call mma_deallocate(imag_amatrix)

call mh5_close_file(restart_fileid)
#endif

return

end subroutine Restart_SurfaceHop
