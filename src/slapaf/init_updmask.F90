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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************

subroutine Init_UpdMask(nInter)

use Slapaf_Info, only: Coor, Curvilinear, Redundant
use NewH_mod, only: UpdMask
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nInter
integer(kind=iwp) :: iAtom, nAtMM, nsAtom
integer, allocatable :: IsMM(:)

nsAtom = size(Coor,2)

! If redundant Cartesians and no symmetry, use unit matrix for MM atoms

if (Redundant .and. (.not. Curvilinear) .and. (3*nsAtom == nInter)) then
  call mma_allocate(UpdMask,nInter,label='UpdMask')
  call mma_allocate(IsMM,nsAtom,Label='IsMM')
  call MMCount(nsAtom,nAtMM,IsMM)
  do iAtom=1,nsAtom
    UpdMask((iAtom-1)*3+1:(iAtom-1)*3+3) = merge(1,0,IsMM(iAtom) == 1)
  end do
  call mma_deallocate(IsMM)
end if

return

end subroutine Init_UpdMask
