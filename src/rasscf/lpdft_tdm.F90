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
! Copyright (C) 2026, Bhavnesh Jangid (University of Chicago)          *
!***********************************************************************

subroutine LPDFT_TDM()
!> Compute active-space transition 1-density matrices and symmetrized transition 2-PDMs for LPDFT
!> and store them in ${Project}.tdm_4_lpdft.h5.
!>
!> Cost note:
!> The 2-PDM is obtained only through Mk_pdms (so all CI backends are supported), by polarization:
!> P[L+R] = P[L] + P[R] + (P_LR + P_RL). For nR = lRoots this needs
!>   - nR  Mk_pdms calls, one per root, stored in Pstate,
!>   - nR*(nR-1)/2  Mk_pdms calls, one per root pair, on the unnormalized sum vector,
!> i.e. nR*(nR+1)/2 CI two-body density evaluations in total, plus one Mk_T1DM call per root pair.
!> The extra memory is NACPR2*nR for Pstate and one CI vector (VecSum).
!> Only the symmetrized part (P_LR + P_RL)/2 is kept, which is all that contracts with the
!> symmetric L-PDFT 2e Hamiltonian; the antisymmetric part is neither computed nor needed.

use input_ras, only: Key
use Definitions, only: u6
#ifdef _HDF5_
use ci_interfaces, only: Mk_pdms, Mk_T1DM
use rasscf_global, only: IADR15, lRoots, NAC, NACPAR, NACPR2
use general_data, only: iDoGAS, NCONF, NELEC3, NHOLE1, STSYM
use rasscf_files, only: JOBIPH
use sguga, only: CIS, SG_ReOrd
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp
use mh5, only: mh5_close_dset, mh5_close_file, mh5_create_dset_int, mh5_create_dset_real, mh5_create_file, mh5_put_dset
#endif

implicit none
#ifdef _HDF5_
integer(kind=iwp) :: dDisk, iDisk, iPair, jRoot, kDisk, kRoot, nPairs
real(kind=wp), allocatable :: Dtmp(:), DStmp(:), Ptmp(:), Pstate(:,:), Tmp(:), VecL(:), VecR(:), VecSum(:), Tdm1(:)
character(len=*), parameter :: LPDFTDM = 'LPDFTDM'
character(len=1024) :: fileName
integer(kind=iwp) :: fileId, lFileName, oneTdmId, twoTdmId, rootPairsId
integer(kind=iwp), allocatable :: rootPairs(:,:)
#endif

! Just the sanity check:
if (.not. Key('LPDF')) return

#ifndef _HDF5_
call WarningMessage(1,'LPDF requested, but HDF5 support is disabled; transition densities will not be stored.')
return
#else

! Only CAS-type active spaces are supported (no GAS/RAS1/RAS3).
if (iDoGAS .or. (NHOLE1 /= 0) .or. (NELEC3 /= 0)) then
  call WarningMessage(2,'LPDF is only available for CASSCF-type active spaces, not for GAS or RAS wave functions.')
  call Quit_OnUserError()
end if

nPairs = (lRoots*(lRoots+1))/2

call PrgmTranslate(LPDFTDM,fileName,lFileName)

write(u6,*)
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A)') '     Multi-State PDFT Intermediate       '
write(u6,'(6X,A)') '-----------------------------------------'
write(u6,'(6X,A)') 'Computing and storing active 1-TDM and 2-PDM matrices in '//fileName(:lFileName)

call mma_allocate(Tmp,NConf,Label='Tmp')
call mma_allocate(VecL,NConf,Label='VecL')
call mma_allocate(VecR,NConf,Label='VecR')
call mma_allocate(VecSum,NConf,Label='VecSum')
call mma_allocate(Tdm1,NAC*NAC,Label='Tdm1')
call mma_allocate(Dtmp,NAC*NAC,Label='Dtmp')
call mma_allocate(DStmp,NAC*NAC,Label='DStmp')
call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
call mma_allocate(Pstate,NACPR2,lRoots,Label='Pstate')
call mma_allocate(rootPairs,2,nPairs,Label='rootPairs')

fileId = mh5_create_file(LPDFTDM)
oneTdmId = mh5_create_dset_real(fileId,'ACTIVE_1TDM',3,[NAC,NAC,nPairs])
twoTdmId = mh5_create_dset_real(fileId,'ACTIVE_2TDM',2,[NACPR2,nPairs])
rootPairsId = mh5_create_dset_int(fileId,'ROOT_PAIRS',2,[2,nPairs])

iPair = 0
iDisk = IADR15(4)
dDisk = IADR15(3)
do jRoot=1,lRoots
  call DDafile(JOBIPH,2,Tmp,nConf,iDisk)
  call SG_ReOrd(1,STSYM,1,CIS(1)%nCSF(STSYM),Tmp,VecL)
  ! 2-PDM of the root, needed for the polarization identity below.
  VecSum(:) = VecL(:)
  call Mk_pdms(VecSum,nConf,P=Pstate(:,jRoot),nD=NAC*NAC,nP=NACPR2)
  kDisk = IADR15(4)
  do kRoot=1,jRoot-1
    call DDafile(JOBIPH,2,Tmp,nConf,kDisk)
    call SG_ReOrd(1,STSYM,1,CIS(1)%nCSF(STSYM),Tmp,VecR)
    call Mk_T1DM(VecR,VecL,nConf,Tdm1,NAC*NAC)
    ! Symmetrized transition 2-PDM from polarization of the (bilinear) 2-PDM:
    ! P[L+R] = P[L] + P[R] + (P_LR + P_RL), so (P_LR + P_RL)/2 = (P[L+R] - P[L] - P[R])/2.
    ! This is the only part of the 2-TDM that contracts with the (symmetric) L-PDFT 2e Hamiltonian.
    VecSum(:) = VecL(:)+VecR(:)
    call Mk_pdms(VecSum,nConf,P=Ptmp,nD=NAC*NAC,nP=NACPR2)
    Ptmp(:) = Half*(Ptmp(:)-Pstate(:,jRoot)-Pstate(:,kRoot))

    call mh5_put_dset(oneTdmId,Tdm1(1:NAC*NAC),[NAC,NAC,1],[0,0,iPair])
    call mh5_put_dset(twoTdmId,Ptmp(1:NACPR2),[NACPR2,1],[0,iPair])

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

  iPair = iPair+1
  rootPairs(1,iPair) = jRoot
  rootPairs(2,iPair) = jRoot
end do

call mh5_put_dset(rootPairsId,rootPairs)
call mh5_close_dset(oneTdmId)
call mh5_close_dset(twoTdmId)
call mh5_close_dset(rootPairsId)
call mh5_close_file(fileId)

write(u6,'(6X,A,I6,A)') 'Stored transition density matrices in HDF5 for ',iPair,' root pairs.'

call mma_deallocate(Tmp)
call mma_deallocate(VecL)
call mma_deallocate(VecR)
call mma_deallocate(VecSum)
call mma_deallocate(Tdm1)
call mma_deallocate(Dtmp)
call mma_deallocate(DStmp)
call mma_deallocate(Ptmp)
call mma_deallocate(Pstate)
call mma_deallocate(rootPairs)

#endif

end subroutine LPDFT_TDM
