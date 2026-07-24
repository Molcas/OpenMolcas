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

use Index_Functions, only: nTri_Elem
use lucia_data, only: CFTP, DStmp, Dtmp, Ptmp
use Lucia_Interface, only: Lucia_Util
use input_ras, only: Key
use rasscf_global, only: IADR15, lRoots, NAC, NACPAR, NACPR2
use general_data, only: JOBIPH, NACTEL, NCONF, STSYM
use sguga_states, only: SGS, EXS, CIS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
#ifdef _HDF5_
use mh5, only: mh5_close_dset, mh5_close_file, mh5_create_dset_int, mh5_create_dset_real, mh5_create_file, mh5_put_dset
#endif

implicit none
integer(kind=iwp) :: dDisk, iDisk, iPair, jRoot, kDisk, kRoot, nPairs, nSpinTdm2
real(kind=wp), allocatable :: SpinTdm2(:), Tmp(:), VecL(:), VecR(:)
#ifdef _HDF5_
character(len=1024) :: currDir, fileName, project
integer(kind=iwp) :: fileId, oneTdmId, twoTdmId, twoTdmSpinId, rootPairsId
integer(kind=iwp), allocatable :: rootPairs(:,:)
#endif

! Just the sanity check:
if (.not. Key('LPDF')) return

#ifndef _HDF5_
call WarningMessage(1,'LPDF requested, but HDF5 support is disabled; transition densities will not be stored.')
return
#else

nPairs = (lRoots*(lRoots+1))/2
nSpinTdm2 = nTri_Elem(NAC*NAC)

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
call mma_allocate(Dtmp,NAC*NAC,Label='Dtmp')
call mma_allocate(DStmp,NAC*NAC,Label='DStmp')
call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
call mma_allocate(SpinTdm2,nSpinTdm2,Label='SpinTdm2')
call mma_allocate(rootPairs,2,nPairs,Label='rootPairs')

fileId = mh5_create_file(trim(fileName))
oneTdmId = mh5_create_dset_real(fileId,'ACTIVE_1TDM',3,[NAC,NAC,nPairs])
twoTdmId = mh5_create_dset_real(fileId,'ACTIVE_2TDM',2,[NACPR2,nPairs])
twoTdmSpinId = mh5_create_dset_real(fileId,'ACTIVE_2TDM_SPINSUM',2,[nSpinTdm2,nPairs])
rootPairsId = mh5_create_dset_int(fileId,'ROOT_PAIRS',2,[2,nPairs])

iPair = 0
iDisk = IADR15(4)
dDisk = IADR15(3)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,Tmp,nConf,iDisk)
  call SG_ReOrd(SGS(1),EXS(1),STSYM,1,CIS(1)%nCSF(STSYM),Tmp,VecL)
  kDisk = IADR15(4)
  do kRoot=1,jRoot-1
    call DDafile(JOBIPH,2,Tmp,nConf,kDisk)
    call SG_ReOrd(SGS(1),EXS(1),STSYM,1,CIS(1)%nCSF(STSYM),Tmp,VecR)
    SpinTdm2(:) = 0.0_wp
    call Lucia_Util('Densi',CI_Vector=VecL(:),RVec=VecR(:),TwoBodyTDM=.true.,Raw2TDM=SpinTdm2(:))
    call scale_spin_tdm2_for_hamiltonian(SpinTdm2,nac)

    call mh5_put_dset(oneTdmId,Dtmp(1:NAC*NAC),[NAC,NAC,1],[0,0,iPair])
    Ptmp(:) = 0.0_wp
    call mh5_put_dset(twoTdmId,Ptmp(1:NACPR2),[NACPR2,1],[0,iPair])
    call mh5_put_dset(twoTdmSpinId,SpinTdm2(1:nSpinTdm2),[nSpinTdm2,1],[0,iPair])

    iPair = iPair+1
    rootPairs(1,iPair) = jRoot
    rootPairs(2,iPair) = kRoot
  end do

  call DDafile(JOBIPH,2,Dtmp,NACPAR,dDisk)
  call DDafile(JOBIPH,2,DStmp,NACPAR,dDisk)
  call DDafile(JOBIPH,2,Ptmp,NACPR2,dDisk)
  call DDafile(JOBIPH,0,Ptmp,NACPR2,dDisk)

  call Square(Dtmp,DStmp,1,NAC,NAC)
  call mh5_put_dset(oneTdmId,DStmp(1:NAC*NAC),[NAC,NAC,1],[0,0,iPair])
  call mh5_put_dset(twoTdmId,Ptmp(1:NACPR2),[NACPR2,1],[0,iPair])
  SpinTdm2(:) = 0.0_wp
  call mh5_put_dset(twoTdmSpinId,SpinTdm2(1:nSpinTdm2),[nSpinTdm2,1],[0,iPair])

  iPair = iPair+1
  rootPairs(1,iPair) = jRoot
  rootPairs(2,iPair) = jRoot
end do

call mh5_put_dset(rootPairsId,rootPairs)
call mh5_close_dset(oneTdmId)
call mh5_close_dset(twoTdmId)
call mh5_close_dset(twoTdmSpinId)
call mh5_close_dset(rootPairsId)
call mh5_close_file(fileId)

write(u6,'(6X,A,I6,A)') 'Stored transition density matrices in HDF5 for ',iPair,' root pairs.'

call mma_deallocate(Tmp)
call mma_deallocate(VecL)
call mma_deallocate(VecR)
call mma_deallocate(Dtmp)
call mma_deallocate(DStmp)
call mma_deallocate(Ptmp)
call mma_deallocate(SpinTdm2)
call mma_deallocate(rootPairs)
#endif

contains

  subroutine scale_spin_tdm2_for_hamiltonian(tdm2,nAct)
    real(kind=wp), intent(inout) :: tdm2(:)
    integer(kind=iwp), intent(in) :: nAct
    integer(kind=iwp) :: ij, ijij, nPairs

    nPairs = nAct*nAct
    ijij = 0
    do ij=1,nPairs
      ijij = ijij+ij
      tdm2(ijij) = 0.5_wp*tdm2(ijij)
    end do
  end subroutine scale_spin_tdm2_for_hamiltonian

end subroutine LPDFT_TDM
