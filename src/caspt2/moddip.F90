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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine ModDip()
! modify the Last Dipole Moments array to avoid crash in
! mclr/out_pt2
! Note that Last Dipole Moments are actually not used anywhere

use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: wp, iwp
use caspt2_module, only: LROOTS, NROOTS, NSTATE, ROOT2STATE
use Constants, only: Zero

implicit none
real(kind=wp), allocatable :: DMs1(:,:), DMs2(:,:)
logical(kind=iwp) :: Found
integer(kind=iwp) :: nData, i, j

if (NSTATE == lRoots) return

call qpg_dArray('Last Dipole Moments',Found,nData)
if (nData == 3*NSTATE) return

if ((.not. Found) .or. (nData /= 3*lRoots)) then
  call WarningMessage(2,'Should not happen in ModDip')
  call abend()
end if

call mma_allocate(DMs1,3,nRoots,Label='DMs1')
DMs1(:,:) = Zero
call mma_allocate(DMs2,3,lRoots,Label='DMs2')
DMs2(:,:) = Zero

call Get_dArray('Last Dipole Moments',DMs2,3*LROOTS)
do i=1,lRoots
  j = Root2State(i)
  if (j == 0) cycle
  DMs1(:,j) = DMs2(:,i)
end do
call Put_dArray('Last Dipole Moments',DMs1,3*nRoots)
call mma_deallocate(DMs1)
call mma_deallocate(DMs2)

return

end subroutine ModDip
