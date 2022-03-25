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

subroutine AllenGinsberg(QMMethod,Eint,Poli,dNuc,Cha,Dip,Qua,iVEC,nDim,lEig,iEig,iQ_Atoms,ip_ExpCento,E_Nuc_Part,lSlater,Eint_Nuc)

use qmstat_global, only: iExtr_Atm
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
character(len=5) :: QMMethod
integer(kind=iwp) :: iVEC, nDim, iEig, iQ_Atoms, ip_ExpCento
real(kind=wp) :: Eint(nTri_Elem(iQ_Atoms),10), Poli(nTri_Elem(iQ_Atoms),10), dNuc(iQ_Atoms), &
                 Cha(nTri_Elem(nDim),nTri_Elem(iQ_Atoms)), Dip(nTri_Elem(nDim),3,nTri_Elem(iQ_Atoms)), &
                 Qua(nTri_Elem(nDim),6,nTri_Elem(iQ_Atoms)), E_Nuc_Part, Eint_Nuc(iQ_Atoms)
logical(kind=iwp) :: lEig, lSlater
integer(kind=iwp) :: i, i1, iAt, iCx, j, jAt, k, kaunt, kaunter, kk, NExpect, NExtrAt, NTotal
real(kind=wp) :: dMp
logical(kind=iwp) :: Check1, Check2
integer(kind=iwp), allocatable :: iCenSet(:)
real(kind=wp), allocatable :: VelP(:), VpoP(:)
#include "warnings.h"

! Set up centre index set. The order of centres are decided by
! the MpProp-program and are hence collected in the get_center
! routine.

! Atom centres

NExtrAt = size(iExtr_Atm)

call mma_allocate(iCenSet,NExtrAt+nTri_Elem(iQ_Atoms-1),label='iCenSet')
iCenSet(1:NExtrAt) = iExtr_Atm(1:NExtrAt)

! Bond centres

kaunter = iQ_Atoms
kaunt = NExtrAt
do iAt=2,iQ_Atoms
  do jAt=1,iAt-1
    kaunter = kaunter+1
    Check1 = .false.
    Check2 = .false.
    do i1=1,NExtrAt
      if (iAt == iCenSet(i1)) Check1 = .true.
      if (jAt == iCenSet(i1)) Check2 = .true.
    end do
    if (Check1 .and. Check2) then
      kaunt = kaunt+1
      iCenSet(kaunt) = kaunter
    end if
  end do
end do

! A minor check.

NExpect = NExtrAt*(nExtrAt+1)/2
NTotal = kaunt
if (NTotal /= NExpect) then
  write(u6,*)
  write(u6,*) ' Error in atom specification for partial perturbation extraction.'
  call Quit(_RC_GENERAL_ERROR_)
end if

! Compute partial nuclear contribution.

E_Nuc_Part = Zero
do iAt=1,NExtrAt
  iCx = iCenSet(iAt)
  if (lSlater) then
    E_Nuc_Part = E_Nuc_Part-(Eint_Nuc(iCx)+Poli(iCx,1))*dNuc(iCx)
  else
    E_Nuc_Part = E_Nuc_Part-(Eint(iCx,1)+Poli(iCx,1))*dNuc(iCx)
  end if
end do

! Set up matrix elements for the partial perturbations.
! Compare with hel, helstate, polink and polins.

call mma_allocate(VelP,nTri_Elem(nDim),label='VelPart')
call mma_allocate(VpoP,nTri_Elem(nDim),label='VpoPart')
VelP(:) = Zero
VpoP(:) = Zero
kk = 0
do i=1,nDim
  do j=1,i
    kk = kk+1
    do k=1,NTotal
      iCx = iCenSet(k)
      dMp = Cha(kk,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,1)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,1)*dMp
      dMp = Dip(kk,1,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,2)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,2)*dMp
      dMp = Dip(kk,2,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,3)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,3)*dMp
      dMp = Dip(kk,3,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,4)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,4)*dMp
      dMp = Qua(kk,1,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,5)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,5)*dMp
      dMp = Qua(kk,3,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,7)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,7)*dMp
      dMp = Qua(kk,6,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,10)*dMp
      VpoP(kk) = VpoP(kk)+Poli(iCx,10)*dMp
      dMp = Qua(kk,2,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,6)*dMp*Two
      VpoP(kk) = VpoP(kk)+Poli(iCx,6)*dMp*Two
      dMp = Qua(kk,4,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,8)*dMp*Two
      VpoP(kk) = VpoP(kk)+Poli(iCx,8)*dMp*Two
      dMp = Qua(kk,5,iCx)
      VelP(kk) = VelP(kk)+Eint(iCx,9)*dMp*Two
      VpoP(kk) = VpoP(kk)+Poli(iCx,9)*dMp*Two
    end do
  end do
end do

call mma_deallocate(iCenSet)

! Collect expectation value for the partial perturbation.

call Expectus(QMMethod,VelP,VelP,VpoP,VpoP,iVEC,nDim,lEig,iEig,ip_ExpCento)

! Deallocate.

call mma_deallocate(VelP)
call mma_deallocate(VpoP)

! Howl

return

end subroutine AllenGinsberg
