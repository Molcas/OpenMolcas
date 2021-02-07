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
! * Surboutine to project out some nuclear coordinates from the       *
! * forces and do dynamics in reduced dimensionality.                 *
! *                                                                   *
! * 03/04/2020                                                        *
! * Morgane Vacher                                                    *
! *                                                                   *
! *********************************************************************

subroutine project_out_for(force,natom)

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
real*8, dimension(natom*3), intent(inout) :: force
real*8, dimension(natom*3) :: force_m
real*8, allocatable :: pcoo(:,:), pcoo_m(:,:)
real*8 :: pforce

call mma_allocate(pcoo,POUT,natom*3)
call mma_allocate(pcoo_m,POUT,natom*3)
call mma_allocate(Mass,natom)

call Get_dArray('Proj_Coord',pcoo,POUT*natom*3)
call GetMassDx(Mass,natom)

! Mass-weight the force vector
do i=1,natom
  do j=1,3
    force_m(3*(i-1)+j) = force(3*(i-1)+j)/sqrt(Mass(i))
  end do
end do

do p=1,POUT
  ! Mass-weight the projection vector
  do i=1,natom
    do j=1,3
      pcoo_m(p,3*(i-1)+j) = pcoo(p,3*(i-1)+j)*sqrt(Mass(i))
    end do
  end do
  ! normalise it (needed or not?)
  pcoo_m(p,:) = pcoo_m(p,:)/sqrt(dot_product(pcoo_m(p,:),pcoo_m(p,:)))
  ! Project out
  pforce = dot_product(pcoo_m(p,:),force_m)
  force_m = force_m-pforce*pcoo_m(p,:)
end do

! Un-Mass-weight the force vector
do i=1,natom
  do j=1,3
    force(3*(i-1)+j) = force_m(3*(i-1)+j)*sqrt(Mass(i))
  end do
end do

call mma_deallocate(pcoo)
call mma_deallocate(pcoo_m)
call mma_deallocate(Mass)

return

end subroutine project_out_for
