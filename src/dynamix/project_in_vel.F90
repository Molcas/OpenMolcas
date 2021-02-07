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
! * velocities and do dynamics in reduced dimensionality.             *
! *                                                                   *
! * 27/05/2020                                                        *
! * Morgane Vacher                                                    *
! *                                                                   *
! *********************************************************************

subroutine project_in_vel(vel,natom)

implicit real*8(a-h,o-z)
#include "prgm.fh"
#include "warnings.fh"
#include "Molcas.fh"
parameter(ROUTINE='VV_Second')
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
integer :: i, j, p
integer :: natom
real*8, allocatable :: Mass(:)
real*8, dimension(natom*3), intent(inout) :: vel
real*8, dimension(natom*3) :: vel_m, newvel_m
real*8, allocatable :: pcoo(:,:), pcoo_m(:,:)
real*8 :: pvel

call mma_allocate(pcoo,PIN,natom*3)
call mma_allocate(pcoo_m,PIN,natom*3)
call mma_allocate(Mass,natom)

call Get_dArray('Keep_Coord',pcoo,PIN*natom*3)
call GetMassDx(Mass,natom)

newvel_m = 0.0d0
! Mass-weight the velocity vector
do i=1,natom
  do j=1,3
    vel_m(3*(i-1)+j) = vel(3*(i-1)+j)*sqrt(Mass(i))
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
  pvel = dot_product(pcoo_m(p,:),vel_m)
  newvel_m = newvel_m+pvel*pcoo_m(p,:)
end do

! Un-Mass-weight the velocity vector
do i=1,natom
  do j=1,3
    vel(3*(i-1)+j) = newvel_m(3*(i-1)+j)/sqrt(Mass(i))
  end do
end do

call mma_deallocate(pcoo)
call mma_deallocate(pcoo_m)
call mma_deallocate(Mass)

return

end subroutine project_in_vel
