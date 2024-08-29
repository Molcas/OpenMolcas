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

function MemSO2_P(iCmp,jCmp,kCmp,lCmp,iAO,jAO,kAO,lAO)
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
integer(kind=iwp) :: MemSO2_P
integer(kind=iwp), intent(in) :: iCmp, jCmp, kCmp, lCmp, iAO, jAO, kAO, lAO
integer(kind=iwp) :: i1, i2, i3, i4, j1, j12, j2, j3, j4

! Quadruple loop over elements of the basis functions angular
! description. Loops are reduced to just produce unique SO integrals
! Observe that we will walk through the memory in AOInt in a
! sequential way.

if (nIrrep == 1) then

  MemSO2_P = iCmp*jCmp*kCmp*lCmp

else

  MemSO2_P = 0

  do i1=1,iCmp
    do i2=1,jCmp
      do i3=1,kCmp
        do i4=1,lCmp

          ! Loop over irreps which are spanned by the basis function.
          ! Again, the loop structure is restricted to ensure unique
          ! integrals.

          do j1=0,nIrrep-1
            if (iAOtSO(iAO+i1,j1) < 0) cycle
            do j2=0,nIrrep-1
              if (iAOtSO(jAO+i2,j2) < 0) cycle
              j12 = ieor(j1,j2)
              do j3=0,nIrrep-1
                if (iAOtSO(kAO+i3,j3) < 0) cycle
                j4 = ieor(j12,j3)
                if (iAOtSO(lAO+i4,j4) < 0) cycle
                MemSO2_P = MemSO2_P+1

              end do
            end do
          end do

        end do
      end do
    end do
  end do

end if

return

end function MemSO2_P
