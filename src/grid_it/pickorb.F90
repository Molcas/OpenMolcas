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

subroutine PickOrb(ipNz,ipSort,ipGref,ipSort_ab,ipGref_ab,ipVol,ipE,ipOcc,ipE_ab,ipOcc_ab,nShowMOs,nShowMOs_ab,isEner,nMOs, &
                   myTitle,ipType)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ipNZ, ipSort, ipGref, ipSort_ab, ipGref_ab, ipVol, ipE, ipOcc, ipE_ab, ipOcc_ab, nMOs, ipType
integer(kind=iwp), intent(out) :: nShowMOs, nShowMOs_ab
integer(kind=iwp), intent(inout) :: isEner
character(len=*), intent(in) :: myTitle
#include "Molcas.fh"
#include "grid.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iActive, ief, ief_ab, ii, ii_ab, iia, iib, ik, ik_ab, il, il_ab, ishift, ispin, j
real(kind=wp) :: ef, ef_ab, eps, R, s

ispin = index(myTitle,' spin')

ik_ab = 0
ii_ab = 0
il_ab = 0
nShowMOs_ab = 0
do i=0,nMOs-1
  iWork(ipSort+i) = 0
  if (isUHF == 1) iWork(ipSort_ab+i) = 0
end do

ishift = 0
do i=0,nIrrep-1
  if (nBas(i) > 0) then
    do j=1,nBas(i)
      iWork(ipNZ+j-1+ishift) = i+1
      iWork(ipNZ+j-1+ishift+nMOs) = j
    end do
  end if
  ishift = ishift+nBas(i)
end do

eps = 1.0e-6_wp
! if no input, but TypeIndex contains 123
if (isAuMO == -1 .and. isAll /= 1) then
  iActive = 0
  do i=0,nMOs-1
    if (iWork(ipType+i) >= 3 .and. iWork(ipType+i) <= 5) iActive = iActive+1
  end do
  if (iActive > 0) then
    ii = 0
    do i=0,nMOs-1
      if (iWork(ipType+i) >= 3 .and. iWork(ipType+i) <= 5) then
        iWork(ipGref+ii) = i+1
        ii = ii+1
      end if
    end do
    nShowMOs = iActive
    goto 555

  end if
end if

if (isAuMo == -1 .and. ispin > 0) then
  isAuMO = 1
  Region(1) = -Two+eps
  Region(2) = -eps
  itRange = 0
  isEner = 0
end if

if (isAll == 0) then
  if (isAuMO == 1 .and. itRange == 1) then
    s = Region(1)
    Region(1) = min(-Region(2),-Region(1))
    Region(2) = max(-Region(2),-s)
    s = Region(1)
    Region(1) = -Region(2)
    Region(2) = -s
  end if
  if (isAuMO == -1 .and. itRange == 0) then
    isAuMO = 1
    Region(1) = -Two+eps
    Region(2) = -eps
  end if
  if (isAuMO == -1 .and. isEner == 1) then
    Region(1) = -1000.0_wp
    Region(2) = 1000.0_wp
  end if
end if
if (isAll == 1) then
  Region(1) = -1000.0_wp
  Region(2) = 1000.0_wp
end if

! 1. user defined number of orbitals. No auto function at all.

if (isAuMO == 0) then
  ishift = 0
  do i=0,nIrrep-1
    if (nBas(i) > 0) then
      iWork(ipSort+i) = ishift
      ! use Sort as temp
    end if
    ishift = ishift+nBas(i)
  end do

  do i=1,nReq
    iia = iReq(i*2-1)
    iib = iReq(i*2)
    if (iia <= 0 .or. iia > nIrrep .or. iib < 0 .or. iib > nBas(iia-1)) then
      write(u6,'(a)') 'Requested orbital does not exist'
      call Quit_OnUserError()

    end if
    iWork(ipGref+i-1) = iWork(ipSort+iia-1)+iib
    if (isUHF == 1) iWork(ipGref_ab+i-1) = iWork(ipGref+i-1)
  end do

  nShowMOs = nReq
  if (isUHF == 1) nShowMOs_ab = nReq
  goto 555
end if
!***********************************************************************
! Well. The user didn't make an exact request. we need to choose orbitals.
if (itRange == 0) then
  isEner = 0
  R = Region(2)
  Region(2) = -Region(1)
  Region(1) = -R
end if

do i=0,nMOs-1
  Work(ipVol+i) = Zero
  iWork(ipSort+i) = 0
  if (isEner == 0) then
    Work(ipE+i) = -Work(ipOcc+i)
    if (isUHF == 1) then
      Work(ipE_ab+i) = -Work(ipOcc_ab+i)
    end if
  end if
end do

! Well, now we need to choose rest (nGrid-1) grids.

! Make stupid sorting...

if (NoSort == 1) then
  ik = 0
  do i=0,nMOs-1
    if (Work(ipE+i) > Region(1) .and. Work(ipE+i) < Region(2)) then
      !write(u6,*) 'EE',Work(ipE+i),Region(1),Region(2)
      ik = ik+1
      iWork(ipSort+i) = ik
    end if
  end do
  ik_ab = ik
else

  ik = 0
  do i=0,nMOs-1
    if (Work(ipE+i) > Region(1) .and. Work(ipE+i) < Region(2)) then
      do j=0,nMOs-1
        if (Work(ipE+j) >= Work(ipE+i) .and. Work(ipE+j) >= Region(1) .and. Work(ipE+j) <= Region(2)) then
          if (Work(ipE+j) == Work(ipE+i)) then
            if (Work(ipOcc+j) <= Work(ipOcc+i)) then
              iWork(ipSort+i) = iWork(ipSort+i)+1
              if (ik < iWork(ipSort+i)) ik = iWork(ipSort+i)
            end if
          else
            iWork(ipSort+i) = iWork(ipSort+i)+1
            if (ik < iWork(ipSort+i)) ik = iWork(ipSort+i)
          end if
        end if
      end do
    end if
  end do

  if (isUHF == 1) then
    ik_ab = 0
    do i=0,nMOs-1
      if (Work(ipE_ab+i) > Region(1) .and. Work(ipE_ab+i) < Region(2)) then
        do j=0,nMOs-1
          if (Work(ipE_ab+j) >= Work(ipE_ab+i) .and. Work(ipE_ab+j) >= Region(1) .and. Work(ipE_ab+j) <= Region(2)) then
            if (Work(ipE_ab+j) == Work(ipE_ab+i)) then
              if (Work(ipOcc_ab+j) <= Work(ipOcc_ab+i)) then
                iWork(ipSort_ab+i) = iWork(ipSort_ab+i)+1
                if (ik_ab < iWork(ipSort_ab+i)) ik_ab = iWork(ipSort_ab+i)
              end if
            else
              iWork(ipSort_ab+i) = iWork(ipSort_ab+i)+1
              if (ik_ab < iWork(ipSort_ab+i)) ik_ab = iWork(ipSort_ab+i)
            end if
          end if
        end do
      end if
    end do
  end if
end if
if (isAuMO == -1 .and. isEner /= 0 .and. isAll == 0) then
  ef = -1000.0_wp
  ief = 0
  ef_ab = ef
  ief_ab = ief
  do i=0,nMOs-1
    if (Work(ipE+i) > ef .and. Work(ipOcc+i) > eps) then
      ef = Work(ipE+i)
      ief = i
    end if
    if (isUHF == 1) then
      if (Work(ipE_ab+i) > ef_ab .and. Work(ipOcc_ab+i) > eps) then
        ef_ab = Work(ipE_ab+i)
        ief_ab = i
      end if
    end if
  end do
  !write(u6,*) 'ef=',ef
  !write(u6,*) 'ief=',ief,ief_ab
  ii = iWork(ipSort+ief)
  if (isUHF == 1) ii_ab = iWork(ipSort_ab+ief_ab)
  do i=0,nMOs-1
    if (iWork(ipSort+i) > ii+iMaxUp .or. iWork(ipSort+i) < ii-iMaxDown) iWork(ipSort+i) = 0
    if (isUHF == 1) then
      if (iWork(ipSort_ab+i) > ii_ab+iMaxUp .or. iWork(ipSort_ab+i) < ii_ab-iMaxDown) iWork(ipSort_ab+i) = 0
    end if
  end do
end if

il = 0
do j=1,ik
  do i=0,nMOs-1
    if (iWork(ipSort+i) == j) then
      iWork(ipGRef+il) = i+1
      il = il+1
    end if
  end do
end do
if (isUHF == 1) then
  il_ab = 0
  do j=1,ik_ab
    do i=0,nMOs-1
      if (iWork(ipSort_ab+i) == j) then
        iWork(ipGRef_ab+il_ab) = i+1
        il_ab = il_ab+1
      end if
    end do
  end do
end if
nShowMOs = il
if (isUHF == 1) nShowMOs_ab = il_ab

555 return

end subroutine PickOrb
