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

subroutine PickOrb(Nz,Sort,Gref,Sort_ab,Gref_ab,E,Occ,E_ab,Occ_ab,nShowMOs,nShowMOs_ab,isEner,nMOs,myTitle,iType)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

#include "intent.fh"

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use grid_it_globals, only: iMaxDown, iMaxUp, iReq, isAll, iAuMO, isUHF, itRange, NoSort, nReq, Region
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMOs, iType(nMOs)
integer(kind=iwp), intent(_OUT_) :: NZ(*)
integer(kind=iwp), intent(out) :: Sort(nMOs), Sort_ab(merge(nMOs,0,isUHF)), nShowMOs, nShowMOs_ab
integer(kind=iwp), intent(inout) :: Gref(nMOs), Gref_ab(merge(nMOs,0,isUHF))
real(kind=wp), intent(inout) :: E(nMOs), E_ab(merge(nMOs,0,isUHF))
real(kind=wp), intent(in) :: Occ(nMOs), Occ_ab(merge(nMOs,0,isUHF))
logical(kind=iwp), intent(inout) :: isEner
character(len=*), intent(in) :: myTitle
integer(kind=iwp) :: i, iActive, ief, ief_ab, ii, ii_ab, iia, iib, ik, ik_ab, il, il_ab, ishift, ispin, j
real(kind=wp) :: ef, ef_ab, eps, R, s

ispin = index(myTitle,' spin')

ik_ab = 0
ii_ab = 0
il_ab = 0
nShowMOs_ab = 0
Sort(:) = 0
if (isUHF) Sort_ab(:) = 0

ishift = 0
do i=0,nIrrep-1
  if (nBas(i) > 0) then
    do j=1,nBas(i)
      NZ(j+ishift) = i+1
      NZ(j+ishift+nMOs) = j
    end do
  end if
  ishift = ishift+nBas(i)
end do

eps = 1.0e-6_wp
! if no input, but TypeIndex contains 123
if ((iAuMO == -1) .and. (.not. isAll)) then
  iActive = 0
  do i=1,nMOs
    if ((iType(i) >= 3) .and. (iType(i) <= 5)) iActive = iActive+1
  end do
  if (iActive > 0) then
    ii = 1
    do i=1,nMOs
      if ((iType(i) >= 3) .and. (iType(i) <= 5)) then
        Gref(ii) = i
        ii = ii+1
      end if
    end do
    nShowMOs = iActive
    return

  end if
end if

if ((iAuMo == -1) .and. (ispin > 0)) then
  iAuMO = 1
  Region(1) = -Two+eps
  Region(2) = -eps
  itRange = 0
  isEner = .false.
end if

if (isAll) then
  Region(1) = -1000.0_wp
  Region(2) = 1000.0_wp
else
  if ((iAuMO == 1) .and. (itRange == 1)) then
    s = Region(1)
    Region(1) = min(-Region(2),-Region(1))
    Region(2) = max(-Region(2),-s)
    s = Region(1)
    Region(1) = -Region(2)
    Region(2) = -s
  end if
  if ((iAuMO == -1) .and. (itRange == 0)) then
    iAuMO = 1
    Region(1) = -Two+eps
    Region(2) = -eps
  end if
  if ((iAuMO == -1) .and. isEner) then
    Region(1) = -1000.0_wp
    Region(2) = 1000.0_wp
  end if
end if

! 1. user defined number of orbitals. No auto function at all.

if (iAuMO == 0) then
  ishift = 0
  do i=0,nIrrep-1
    if (nBas(i) > 0) then
      Sort(i+1) = ishift
      ! use Sort as temp
    end if
    ishift = ishift+nBas(i)
  end do

  do i=1,nReq
    iia = iReq(i*2-1)
    iib = iReq(i*2)
    if ((iia <= 0) .or. (iia > nIrrep) .or. (iib < 0) .or. (iib > nBas(iia-1))) then
      write(u6,'(a)') 'Requested orbital does not exist'
      call Quit_OnUserError()

    end if
    Gref(i) = Sort(iia)+iib
    if (isUHF) Gref_ab(i) = Gref(i)
  end do

  nShowMOs = nReq
  if (isUHF) nShowMOs_ab = nReq
  return
end if
!***********************************************************************
! Well. The user didn't make an exact request. we need to choose orbitals.
if (itRange == 0) then
  isEner = .false.
  R = Region(2)
  Region(2) = -Region(1)
  Region(1) = -R
end if

Sort(:) = 0
if (.not. isEner) then
  E(:) = -Occ(:)
  if (isUHF) E_ab(:) = -Occ_ab(:)
end if

! Well, now we need to choose rest (nGrid-1) grids.

! Make stupid sorting...

if (NoSort) then
  ik = 0
  do i=1,nMOs
    if ((E(i) > Region(1)) .and. (E(i) < Region(2))) then
      !write(u6,*) 'EE',E(i),Region(1),Region(2)
      ik = ik+1
      Sort(i) = ik
    end if
  end do
  ik_ab = ik
else

  ik = 0
  do i=1,nMOs
    if ((E(i) > Region(1)) .and. (E(i) < Region(2))) then
      do j=1,nMOs
        if ((E(j) >= E(i)) .and. (E(j) >= Region(1)) .and. (E(j) <= Region(2))) then
          if (E(j) == E(i)) then
            if (Occ(j) <= Occ(i)) then
              Sort(i) = Sort(i)+1
              if (ik < Sort(i)) ik = Sort(i)
            end if
          else
            Sort(i) = Sort(i)+1
            if (ik < Sort(i)) ik = Sort(i)
          end if
        end if
      end do
    end if
  end do

  if (isUHF) then
    ik_ab = 0
    do i=1,nMOs
      if ((E_ab(i) > Region(1)) .and. (E_ab(i) < Region(2))) then
        do j=1,nMOs
          if ((E_ab(j) >= E_ab(i)) .and. (E_ab(j) >= Region(1)) .and. (E_ab(j) <= Region(2))) then
            if (E_ab(j) == E_ab(i)) then
              if (Occ_ab(j) <= Occ_ab(i)) then
                Sort_ab(i) = Sort_ab(i)+1
                if (ik_ab < Sort_ab(i)) ik_ab = Sort_ab(i)
              end if
            else
              Sort_ab(i) = Sort_ab(i)+1
              if (ik_ab < Sort_ab(i)) ik_ab = Sort_ab(i)
            end if
          end if
        end do
      end if
    end do
  end if
end if
if ((iAuMO == -1) .and. isEner .and. (.not. isAll)) then
  ef = -1000.0_wp
  ief = 1
  ef_ab = ef
  ief_ab = ief
  do i=1,nMOs
    if ((E(i) > ef) .and. (Occ(i) > eps)) then
      ef = E(i)
      ief = i
    end if
    if (isUHF) then
      if ((E_ab(i) > ef_ab) .and. (Occ_ab(i) > eps)) then
        ef_ab = E(i)
        ief_ab = i
      end if
    end if
  end do
  !write(u6,*) 'ef=',ef
  !write(u6,*) 'ief=',ief,ief_ab
  ii = Sort(ief)
  if (isUHF) ii_ab = Sort_ab(ief_ab)
  do i=1,nMOs
    if ((Sort(i) > ii+iMaxUp) .or. (Sort(i) < ii-iMaxDown)) Sort(i) = 0
    if (isUHF) then
      if ((Sort_ab(i) > ii_ab+iMaxUp) .or. (Sort_ab(i) < ii_ab-iMaxDown)) Sort_ab(i) = 0
    end if
  end do
end if

il = 0
do j=1,ik
  do i=1,nMOs
    if (Sort(i) == j) then
      il = il+1
      GRef(il) = i
    end if
  end do
end do
if (isUHF) then
  il_ab = 0
  do j=1,ik_ab
    do i=1,nMOs
      if (Sort_ab(i) == j) then
        il_ab = il_ab
        GRef_ab(il_ab+1) = i
      end if
    end do
  end do
end if
nShowMOs = il
if (isUHF) nShowMOs_ab = il_ab

return

end subroutine PickOrb
