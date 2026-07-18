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
! Copyright (C) 2024, Bhavnesh Jangid                                  *
!***********************************************************************

subroutine LPDFT_TDM()
!> Compute active-space transition 1- and 2-density matrices for LPDFT
!> and store them in ${Project}_tdm_4_lpdft.h5.

use csfbas, only: CONF
use lucia_data, only: CFTP, DStmp, Dtmp, PAtmp, Ptmp
use Lucia_Interface, only: Lucia_Util
use input_ras, only: Key
use rasscf_global, only: IADR15, lRoots, NAC, NACPR2
use general_data, only: JOBIPH, NACTEL, NCONF, STSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
#ifdef _HDF5_
use mh5, only: mh5_close_dset, mh5_close_file, mh5_create_dset_int, mh5_create_dset_real, mh5_create_file, mh5_put_dset
#endif

implicit none
integer(kind=iwp) :: iDisk, iPair, jRoot, kDisk, kRoot, nPairs
integer(kind=iwp), allocatable :: kCnf(:)
real(kind=wp), allocatable :: Tmp(:), VecL(:), VecR(:)
#ifdef _HDF5_
character(len=1024) :: currDir, fileName, project
integer(kind=iwp) :: fileId, oneTdmId, twoTdmId, twoTdmAId, rootPairsId
integer(kind=iwp), allocatable :: rootPairs(:,:)
#endif

! Just the sanity check:
if (.not. Key('LPDF')) return

if (lRoots <= 1) then
  call WarningMessage(1,'LPDF requested, but only one RASSCF root is available; transition densities will not be stored.')
  return
end if

#ifndef _HDF5_
call WarningMessage(1,'LPDF requested, but HDF5 support is disabled; transition densities will not be stored.')
return
#else

nPairs = (lRoots*(lRoots-1))/2

call get_environment_variable('Project',project)

if (len_trim(project) == 0) project = 'Noname'
call get_environment_variable('CurrDir',currDir)
if (len_trim(currDir) > 0) then
  fileName = trim(currDir)//'/'//trim(project)//'_tdm_4_lpdft.h5'
else
  fileName = trim(project)//'_tdm_4_lpdft.h5'
end if

write(u6,*)
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A)') '     Multi-State PDFT Intermediate       '
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A)') 'Computing and storing active 1-TDM and 2-TDM matrices in '//trim(fileName)

call mma_allocate(Tmp,NConf,Label='Tmp')
call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(kCnf,NACTEL,Label='kCnf')
call mma_allocate(Dtmp,NAC*NAC,Label='Dtmp')
call mma_allocate(DStmp,NAC*NAC,Label='DStmp')
call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
call mma_allocate(rootPairs,2,nPairs,Label='rootPairs')

fileId = mh5_create_file(trim(fileName))
oneTdmId = mh5_create_dset_real(fileId,'ACTIVE_1TDM',3,[NAC,NAC,nPairs])
twoTdmId = mh5_create_dset_real(fileId,'ACTIVE_2TDM',2,[NACPR2,nPairs])
twoTdmAId = mh5_create_dset_real(fileId,'ACTIVE_2TDM_ANTISYM',2,[NACPR2,nPairs])
rootPairsId = mh5_create_dset_int(fileId,'ROOT_PAIRS',2,[2,nPairs])

iPair = 0
iDisk = IADR15(4)
call DDafile(JOBIPH,2,Tmp,nConf,iDisk)
do jRoot=2,lRoots
  call DDafile(JOBIPH,2,Tmp,nConf,iDisk)
  call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp,VecL,kCnf)
  kDisk = IADR15(4)
  do kRoot=1,jRoot-1
    call DDafile(JOBIPH,2,Tmp,nConf,kDisk)
    call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp,VecR,kCnf)
    call Lucia_Util('Densi',CI_Vector=VecL(:),RVec=VecR(:),TwoBodyTDM=.true.)

    call mh5_put_dset(oneTdmId,Dtmp(1:NAC*NAC),[NAC,NAC,1],[0,0,iPair])
    call mh5_put_dset(twoTdmId,Ptmp(1:NACPR2),[NACPR2,1],[0,iPair])
    call mh5_put_dset(twoTdmAId,PAtmp(1:NACPR2),[NACPR2,1],[0,iPair])

    iPair = iPair+1
    rootPairs(1,iPair) = jRoot
    rootPairs(2,iPair) = kRoot
  end do
end do

call mh5_put_dset(rootPairsId,rootPairs)
call mh5_close_dset(oneTdmId)
call mh5_close_dset(twoTdmId)
call mh5_close_dset(twoTdmAId)
call mh5_close_dset(rootPairsId)
call mh5_close_file(fileId)

write(u6,'(6X,A,I6,A)') 'Stored transition density matrices in HDF5 for ',iPair,' root pairs.'

call mma_deallocate(Tmp)
call mma_deallocate(VecL)
call mma_deallocate(VecR)
call mma_deallocate(kCnf)
call mma_deallocate(Dtmp)
call mma_deallocate(DStmp)
call mma_deallocate(Ptmp)
call mma_deallocate(PAtmp)
call mma_deallocate(rootPairs)
#endif

end subroutine LPDFT_TDM
