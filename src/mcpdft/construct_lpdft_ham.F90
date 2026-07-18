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

  subroutine construct_lpdft_ham_from_tdm(CMO)

  use general_data, only: nBas, nFro, nIsh, nOrb, nAsh, nSym, nTot1, nTot2
  use Index_Functions, only: iTri, nTri_Elem
  use rasscf_global, only: ExFac, lRoots, lSquare, NAC, NACPR2, nFint
  use mspdft, only: heff
  use Constants, only: Half, Zero
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp, u6
#ifdef _HDF5_
  use mh5, only: mh5_close_file, mh5_exists_dset, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
#endif

  implicit none
#include "warnings.h"
  real(kind=wp), intent(in) :: CMO(nTot2)
  integer(kind=iwp) :: iPair, iRoot, jRoot, nPairs
  real(kind=wp) :: e1, e2
  real(kind=wp), allocatable :: h1e(:,:), h2e(:,:,:,:), TDM1(:,:,:), TDM2(:,:), TDM2Full(:,:,:,:)
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
  call mma_allocate(rootPairs,2,nPairs,Label='rootPairs')
  call mh5_fetch_dset(mh5id,'ACTIVE_1TDM',TDM1)
  call mh5_fetch_dset(mh5id,'ACTIVE_2TDM',TDM2)

  if (mh5_exists_dset(mh5id,'ROOT_PAIRS')) then
    call mh5_fetch_dset(mh5id,'ROOT_PAIRS',rootPairs)
  else
    call default_root_pairs(rootPairs,nPairs)
  end if

  call mh5_close_file(mh5id)

  call mma_allocate(h1e,NAC,NAC,Label='h1e')
  call mma_allocate(h2e,NAC,NAC,NAC,NAC,Label='h2e')
  call mma_allocate(TDM2Full,NAC,NAC,NAC,NAC,Label='TDM2Full')

  call active_hcore(CMO,h1e)
  call active_2e_integrals(CMO,h2e)

  do iPair=1,nPairs
    iRoot = rootPairs(1,iPair)
    jRoot = rootPairs(2,iPair)
    if ((iRoot < 1) .or. (iRoot > lRoots) .or. (jRoot < 1) .or. (jRoot > lRoots)) then
      call WarningMessage(2,'LPDF transition density HDF5 file contains an invalid ROOT_PAIRS entry.')
      call Quit_OnUserError()
    end if
    call unpack_2e_8fold(TDM2(:,iPair),TDM2Full)
    e1 = contract_1e(h1e,TDM1(:,:,iPair))
    e2 = contract_2e(h2e,TDM2Full)
    heff(iRoot,jRoot) = e1+Half*e2
    heff(jRoot,iRoot) = heff(iRoot,jRoot)
  end do

  write(u6,'(6X,A)') 'LPDF: read transition density matrices from '//trim(fileName)
  write(u6,'(6X,A,I6,A)') 'LPDF: constructed effective Hamiltonian from ',nPairs,' root pairs.'

  call mma_deallocate(TDM2Full)
  call mma_deallocate(h2e)
  call mma_deallocate(h1e)
  call mma_deallocate(rootPairs)
  call mma_deallocate(TDM1)
  call mma_deallocate(TDM2)
#endif

  end subroutine construct_lpdft_ham_from_tdm

  subroutine default_root_pairs(rootPairs,nPairs)

  use rasscf_global, only: lRoots
  use Definitions, only: iwp

  implicit none
  integer(kind=iwp), intent(in) :: nPairs
  integer(kind=iwp), intent(out) :: rootPairs(2,nPairs)
  integer(kind=iwp) :: iPair, iRoot, jRoot

  iPair = 0
  do iRoot=1,lRoots-1
    do jRoot=iRoot+1,lRoots
      iPair = iPair+1
      rootPairs(1,iPair) = iRoot
      rootPairs(2,iPair) = jRoot
    end do
  end do

  end subroutine default_root_pairs

  subroutine active_hcore(CMO,h1e)

  use general_data, only: nBas, nFro, nIsh, nOrb, nAsh, nSym, nTot1, nTot2
  use Index_Functions, only: iTri, nTri_Elem
  use rasscf_global, only: NAC
  use Constants, only: Zero
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(in) :: CMO(nTot2)
  real(kind=wp), intent(out) :: h1e(NAC,NAC)
  integer(kind=iwp) :: iAct, iGlobal, iMO, iOffMO, iSym, jAct, jGlobal, jMO, nActive
  real(kind=wp), allocatable :: hcore(:), hcoreMO(:)

  call mma_allocate(hcore,nTot1,Label='hcore')
  call mma_allocate(hcoreMO,nTot1,Label='hcoreMO')
  call get_hcore(hcore)
  hcoreMO(:) = Zero
  call ao2mo_1e(CMO,hcore,hcoreMO,nSym,nBas,nOrb,nFro)

  h1e(:,:) = Zero
  iOffMO = 1
  iGlobal = 0
  do iSym=1,nSym
    nActive = nAsh(iSym)
    do iAct=1,nActive
      iMO = nIsh(iSym)+iAct
      do jAct=1,nActive
        jMO = nIsh(iSym)+jAct
        jGlobal = iGlobal+jAct
        h1e(iGlobal+iAct,jGlobal) = hcoreMO(iOffMO+iTri(iMO,jMO)-1)
      end do
    end do
    iOffMO = iOffMO+nTri_Elem(nOrb(iSym))
    iGlobal = iGlobal+nActive
  end do

  call mma_deallocate(hcoreMO)
  call mma_deallocate(hcore)

  end subroutine active_hcore

  subroutine active_2e_integrals(CMO,h2e)

  use PrintLevel, only: DEBUG, INSANE
  use mcpdft_output, only: IPRLOC
  use general_data, only: nTot1, nTot2
  use rasscf_global, only: ExFac, lSquare, NAC, NACPR2, nFint
  use Constants, only: Zero
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(in) :: CMO(nTot2)
  real(kind=wp), intent(out) :: h2e(NAC,NAC,NAC,NAC)
  integer(kind=iwp) :: iprLocal
  real(kind=wp), allocatable :: D1A(:), D1I(:), FA(:), FI(:), PUVX(:), TUVX(:)

  call mma_allocate(D1I,nTot2,Label='D1I')
  call mma_allocate(D1A,nTot2,Label='D1A')
  call mma_allocate(FI,nTot1,Label='FI')
  call mma_allocate(FA,nTot1,Label='FA')
  call mma_allocate(PUVX,nFint,Label='PUVX')
  call mma_allocate(TUVX,NACPR2,Label='TUVX')

  D1I(:) = Zero
  D1A(:) = Zero
  FI(:) = Zero
  FA(:) = Zero
  PUVX(:) = Zero
  TUVX(:) = Zero

  iprLocal = 0
  if (IPRLOC(2) == DEBUG) iprLocal = 5
  if (IPRLOC(2) == INSANE) iprLocal = 10
  call TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,iprLocal,lSquare,ExFac)
  call unpack_2e_8fold(TUVX,h2e)

  call mma_deallocate(TUVX)
  call mma_deallocate(PUVX)
  call mma_deallocate(FA)
  call mma_deallocate(FI)
  call mma_deallocate(D1A)
  call mma_deallocate(D1I)

  end subroutine active_2e_integrals

  subroutine unpack_2e_8fold(packed,full)

  use rasscf_global, only: NAC, NACPR2
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(in) :: packed(NACPR2)
  real(kind=wp), intent(out) :: full(NAC,NAC,NAC,NAC)
  integer(kind=iwp) :: i, idx, j, k, l, lMax

  idx = 0
  do i=1,NAC
    do j=1,i
      do k=1,i
        lMax = k
        if (i == k) lMax = j
        do l=1,lMax
          idx = idx+1
          full(i,j,k,l) = packed(idx)
          full(j,i,k,l) = packed(idx)
          full(i,j,l,k) = packed(idx)
          full(j,i,l,k) = packed(idx)
          full(k,l,i,j) = packed(idx)
          full(l,k,i,j) = packed(idx)
          full(k,l,j,i) = packed(idx)
          full(l,k,j,i) = packed(idx)
        end do
      end do
    end do
  end do

  end subroutine unpack_2e_8fold

  function contract_1e(h1e,tdm1) result(value)

  use rasscf_global, only: NAC
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(in) :: h1e(NAC,NAC), tdm1(NAC,NAC)
  real(kind=wp) :: value
  integer(kind=iwp) :: i, j

  value = Zero
  do i=1,NAC
    do j=1,NAC
      value = value+h1e(i,j)*tdm1(i,j)
    end do
  end do

  end function contract_1e

  function contract_2e(h2e,tdm2) result(value)

  use rasscf_global, only: NAC
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  implicit none
  real(kind=wp), intent(in) :: h2e(NAC,NAC,NAC,NAC), tdm2(NAC,NAC,NAC,NAC)
  real(kind=wp) :: value
  integer(kind=iwp) :: i, j, k, l

  value = Zero
  do i=1,NAC
    do j=1,NAC
      do k=1,NAC
        do l=1,NAC
          value = value+h2e(i,j,k,l)*tdm2(i,j,k,l)
        end do
      end do
    end do
  end do

  end function contract_2e

end module construct_lpdft_ham
