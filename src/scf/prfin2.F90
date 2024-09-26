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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

subroutine PrFin2(Ovlp,nDT,OccNo,nEO,CMO,nCMO,note)

use InfSCF, only: iCoCo, jVOut, kIVO, KSDFT, nBas, nBB, nBT, nD, nnB, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
! PAM 2007: Changed dimension of CMO array from nCMO to NBB:
! The larger size is needed here, and the allocated size is nBB.
integer(kind=iwp) :: nDT, nEO, nCMO
real(kind=wp) :: Ovlp(nDT), OccNo(nEO), CMO(nBB)
character(len=80) :: Note
integer(kind=iwp) :: i, iBs, iCMO, i_Or, iSym, iVec, j
real(kind=wp), allocatable :: Scr2(:)

! Write orbitals on the file (the case InVec=3 and nIter=0
! is set up in RdInp)
if (jVOut >= 2) then

  ! Allocate memory for expanded CMO's and occupation numbers
  call mma_allocate(Scr2,nBB,Label='Scr2')

  ! Prepare CMO in symmetry blocks nBas x nBas
  iVec = 0
  iCMO = 0
  !iseed = 153759
  do iSym=1,nSym
    do i=1,nBas(iSym)*nOrb(iSym)
      Scr2(iVec+i) = CMO(iCMO+i)
    end do
    iVec = iVec+nBas(iSym)*nOrb(iSym)
    do i=1,nBas(iSym)*(nBas(iSym)-nOrb(iSym))
      !Scr2(iVec + i) = Random(iseed)-Half
      Scr2(iVec+i) = Zero
    end do
    iVec = iVec+nBas(iSym)*(nBas(iSym)-nOrb(iSym))
    iCMO = iCMO+nOrb(iSym)*nBas(iSym)
  end do
  do i=1,nBB
    CMO(i) = Scr2(i)
  end do

  ! Prepare occupation numbers
  i_Or = 0
  iBs = 0
  do iSym=1,nSym
    do j=1,nOrb(iSym)
      Scr2(iBs+j) = OccNo(i_Or+j)
    end do
    do j=nOrb(iSym)+1,nBas(iSym)
      Scr2(iBs+j) = Zero
    end do
    i_Or = i_Or+nOrb(iSym)
    iBs = iBs+nBas(iSym)
  end do
  do i=1,nnB
    OccNo(i) = Scr2(i)
  end do

  do iSym=1,nSym
    nOrb(iSym) = nBas(iSym)
  end do
  call SetUp_SCF()

  ! Orthogonalize vectors
  call Ortho(CMO,nCMO,Ovlp,nBT)

  ! Write on the file
  if (KSDFT == 'SCF') then
    if (nD == 1) then
      Note = '* SCF orbitals'
      if (kIVO /= 0) Note = '* SCF orbitals + IVO'
      if (iCoCo /= 0) Note = '* SCF orbitals + arbitrary occupations'
    else
      Note = '* UHF orbitals'
      if (kIVO /= 0) Note = '* UHF orbitals + IVO'
      if (iCoCo /= 0) Note = '* UHF orbitals + arbitrary occupations'
    end if
  else
    if (nD == 1) then
      Note = '* RKS-DFT orbitals'
      if (kIVO /= 0) Note = '* RKS-DFT orbitals + IVO'
      if (iCoCo /= 0) Note = '* RKS-DFT orbitals + arbitrary occupations'
    else
      Note = '* UKS-DFT orbitals'
      if (kIVO /= 0) Note = '* UKS-DFT orbitals + IVO'
      if (iCoCo /= 0) Note = '* UKS-DFT orbitals + arbitrary occupations'
    end if
  end if

  ! Deallocate memory for expanded CMO's and occupation numbers
  call mma_deallocate(Scr2)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine PrFin2
