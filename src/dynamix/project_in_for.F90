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
! Copyright (C) 2020, Morgane Vacher                                   *
!***********************************************************************
!
! *********************************************************************
! *                                                                   *
! * Subroutine to keep in only some nuclear coordinates for the       *
! * forces and do dynamics in reduced dimensionality.                 *
! *                                                                   *
! * 27/05/2020                                                        *
! * Morgane Vacher                                                    *
! *                                                                   *
! *********************************************************************

subroutine project_in_for(force,natom)

use Dynamix_Globals, only: PIN
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(inout) :: force(natom*3)
integer(kind=iwp) :: i, j, p
real(kind=wp) :: pforce
real(kind=wp), allocatable :: Mass(:), pcoo(:,:), pcoo_m(:,:), force_m(:), newforce_m(:)

call mma_allocate(force_m,natom*3)
call mma_allocate(newforce_m,natom*3)
call mma_allocate(pcoo,PIN,natom*3)
call mma_allocate(pcoo_m,PIN,natom*3)
call mma_allocate(Mass,natom)

call Get_dArray('Keep_Coord',pcoo,PIN*natom*3)
call GetMassDx(Mass,natom)

newforce_m(:) = Zero
! Mass-weight the force vector
do i=1,natom
  do j=1,3
    force_m(3*(i-1)+j) = force(3*(i-1)+j)/sqrt(Mass(i))
  end do
end do

do p=1,PIN
  ! Mass-weight the projection vector
  do i=1,natom
    do j=1,3
      pcoo_m(p,3*(i-1)+j) = pcoo(p,3*(i-1)+j)*sqrt(Mass(i))
    end do
  end do
  ! normalise it (needed or not?)
  pcoo_m(p,:) = pcoo_m(p,:)/sqrt(dot_product(pcoo_m(p,:),pcoo_m(p,:)))
  ! Calculate the projection to keep in
  pforce = dot_product(pcoo_m(p,:),force_m)
  newforce_m(:) = newforce_m(:)+pforce*pcoo_m(p,:)
end do

! Un-Mass-weight the force vector
do i=1,natom
  do j=1,3
    force(3*(i-1)+j) = newforce_m(3*(i-1)+j)*sqrt(Mass(i))
  end do
end do

call mma_deallocate(force_m)
call mma_deallocate(newforce_m)
call mma_deallocate(pcoo)
call mma_deallocate(pcoo_m)
call mma_deallocate(Mass)

return

end subroutine project_in_for
