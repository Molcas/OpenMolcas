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

subroutine CoreToPoint(nAt,MP,TP)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt
real(kind=wp), intent(inout) :: MP(nAt*(nAt+1)/2)
real(kind=wp), intent(out) :: TP(nAt)
integer(kind=iwp) :: i, iAt, kAt
real(kind=wp) :: ByggMeraHus, dScaleOff, dScaleOffSave, dToPoint
logical(kind=iwp) :: GoHere
real(kind=wp), allocatable :: NucCh(:)
real(kind=wp), parameter :: ByggareBas(6) = [2.0_wp,8.0_wp,8.0_wp,18.0_wp,18.0_wp,32.0_wp]

! A crude but to the point algorithm to put core electrons and
! nuclei together, and separating the presumably more diffuse
! part of the charge distribution in a separate chunk.

call mma_allocate(NucCh,nAt,label='NucC')
call Get_dArray('Nuclear charge',NucCh,nAt)
kAt = 0
dToPoint = Zero
do iAt=1,nAt
  ByggMeraHus = Zero
  GoHere = .true.
  dScaleOffSave = NucCh(iAt)
  do i=1,6
    dScaleOff = dScaleOffSave-ByggareBas(i)
    if ((dScaleOff <= Zero) .and. GoHere) then
      dToPoint = ByggMeraHus
      GoHere = .false.
    end if
    dScaleOffSave = dScaleOff
    ByggMeraHus = ByggMeraHus+ByggareBas(i)
  end do
  kAt = kAt+iAt
  MP(kAt) = MP(kAt)+dToPoint
  TP(iAt) = NucCh(iAt)-dToPoint
end do
call mma_deallocate(NucCh)

return

end subroutine CoreToPoint
