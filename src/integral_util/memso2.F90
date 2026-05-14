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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

function MemSO2(nSD,iSD4)
!***********************************************************************
!  Object: to compile the number of SO block which will be generated   *
!          by the current shell quadruplet.                            *
!                                                                      *
!          Observe that the indices are canonically ordered at the     *
!          time of calling this routine!                               *
!                                                                      *
! Called from: Drv2El                                                  *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             February '90                                             *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MemSO2
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4)
integer(kind=iwp) :: i1, i2, i3, i4, iAO, iCmp, iShell, j1, j12, j2, j2Max, j3, j3Max, j4, jAO, jCmp, jCmpMx, jShell, kAO, kCmp, &
                     kCmpMx, kShell, lAO, lCmp, lCmpMx, lShell
logical(kind=iwp) :: Shij, Shik, Shjl, Shkl

MemSO2 = 0

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

iCmp = iSD4(2,1)
jCmp = iSD4(2,2)
kCmp = iSD4(2,3)
lCmp = iSD4(2,4)
iShell = iSD4(11,1)
jShell = iSD4(11,2)
kShell = iSD4(11,3)
lShell = iSD4(11,4)
iAO = iSD4(7,1)
jAO = iSD4(7,2)
kAO = iSD4(7,3)
lAO = iSD4(7,4)

Shij = (iShell == jShell)
Shkl = (kShell == lShell)
Shik = (iShell == kShell)
Shjl = (jShell == lShell)

if (nIrrep == 1) then

  do i1=1,iCmp
    jCmpMx = jCmp
    if (Shij) jCmpMx = i1
    do i2=1,jCmpMx
      kCmpMx = kCmp
      if (Shik .and. Shjl) kCmpMx = i1
      do i3=1,kCmpMx
        lCmpMx = lCmp
        if (Shkl) lCmpMx = i3
        if (Shik .and. (i1 == i3) .and. Shjl) lCmpMx = i2
        MemSO2 = MemSO2+lCmpMx
      end do
    end do
  end do

else

  do i1=1,iCmp
    jCmpMx = jCmp
    if (Shij) jCmpMx = i1
    do i2=1,jCmpMx
      kCmpMx = kCmp
      if (Shik .and. Shjl) kCmpMx = i1
      do i3=1,kCmpMx
        lCmpMx = lCmp
        if (Shkl) lCmpMx = i3
        if (Shik .and. (i1 == i3) .and. Shjl) lCmpMx = i2
        do i4=1,lCmpMx

          ! Loop over irreps which are spanned by the basis function.
          ! Again, the loop structure is restricted to ensure unique
          ! integrals.

          do j1=0,nIrrep-1
            if (iAOtSO(iAO+i1,j1) < 0) cycle
            j2Max = nIrrep-1
            if (Shij .and. (i1 == i2)) j2Max = j1
            do j2=0,j2Max
              if (iAOtSO(jAO+i2,j2) < 0) cycle
              j12 = ieor(j1,j2)
              j3Max = nIrrep-1
              if (Shik .and. (i1 == i3) .and. Shjl .and. (i2 == i4)) j3Max = j1
              do j3=0,j3Max
                if (iAOtSO(kAO+i3,j3) < 0) cycle
                j4 = ieor(j12,j3)
                if (iAOtSO(lAO+i4,j4) < 0) cycle
                if (Shkl .and. (i3 == i4) .and. (j4 > j3)) cycle
                if (Shik .and. (i1 == i3) .and. Shjl .and. (i2 == i4) .and. (j1 == j3) .and. (j4 > j2)) cycle
                MemSO2 = MemSO2+1

              end do
            end do
          end do

        end do
      end do
    end do
  end do

end if

return

end function MemSO2
