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

subroutine Shell_MxSchwz(nSkal,Schwz_Shl)
!***********************************************************************
!                                                                      *
!  Subroutine Shell_MxSchwz:  gets max integral estimates for each     *
!                             shell pair...                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use iSD_data, only: iSD
use Basis_Info, only: DBSC, Shells
use k2_structure, only: IndK2, k2Data
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSkal
real(kind=wp), intent(out) :: Schwz_Shl(nSkal,nSkal)
integer(kind=iwp) :: iCnttp, ijS, ik2, iS, iShell, iShll, jCnttp, jS, jShell, jShll, lDCRR, nDCRR, nij
real(kind=wp) :: Schwz_Tmp
integer(kind=iwp), allocatable :: Pair_Index(:,:)

! loop over shell pair...
Schwz_Shl(:,:) = Zero
do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux .and. (iS /= nSkal)) cycle
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  do jS=1,iS
    jShll = iSD(0,jS)
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle
    if (Shells(jShll)%Aux .and. (jS == nSkal)) cycle
    !write(u6,*) 'Shell_..:iS,jS=',iS,jS
    jShell = iSD(11,jS)
    jCnttp = iSD(13,jS)
    ijS = iTri(iShell,jShell)
    nDCRR = Indk2(2,ijS)
    ik2 = Indk2(3,ijS)
    ! now loop over  R operator...
    if (dbsc(iCnttp)%fMass == dbsc(jCnttp)%fMass) then
      Schwz_tmp = k2data(1,ik2)%EstI
      do lDCRR=2,nDCRR
        Schwz_tmp = max(Schwz_tmp,k2data(lDCRR,ik2)%EstI)
      end do
    else
      Schwz_tmp = Zero
    end if
    Schwz_Shl(jS,iS) = Schwz_tmp
    Schwz_Shl(iS,jS) = Schwz_tmp
  end do
end do
!call RecPrt('Schwz_shl',' ',Schwz_Shl,nSkal,nSkal)
call mma_allocate(Pair_Index,2,nTri_Elem(nSkal),Label='Pair_Index')
nij = 0
do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux .and. (iS /= nSkal)) cycle
  do jS=1,iS
    jShll = iSD(0,jS)
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle
    if (Shells(jShll)%Aux .and. (jS == nSkal)) cycle
    nij = nij+1
    Pair_Index(1,nij) = iS
    Pair_Index(2,nij) = jS
  end do
end do

call Drv2el_ijij(Pair_Index,nij,Schwz_shl,nSkal)

call mma_deallocate(Pair_Index)

end subroutine Shell_MxSchwz
