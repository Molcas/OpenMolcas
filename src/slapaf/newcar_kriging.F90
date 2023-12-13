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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine NewCar_Kriging(kIter,SaveBMx,Error)

use Slapaf_Info, only: BMx, BMx_kriging, Cx, Force_dB, mTtAtm, Numerical, PrQ
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: kIter
logical(kind=iwp), intent(in) :: SaveBMx
logical(kind=iwp), intent(inout) :: Error
logical(kind=iwp) :: Numerical_Save, PrQ_Save
real(kind=wp), allocatable :: BMx_Tmp(:,:), Coor(:,:)

call mma_allocate(Coor,3,size(Cx,2),Label='Coor')
Coor(:,:) = Cx(:,:,kIter)

call mma_allocate(BMx_tmp,size(BMx,1),size(BMx,2),Label='BMx_tmp')
BMx_tmp(:,:) = BMx(:,:)

Numerical_Save = Numerical
PrQ_Save = PrQ
Numerical = .false.
PrQ_Save = PrQ
PrQ = .false.

Force_dB = SaveBMx

call NewCar(kIter,size(Coor,2),Coor,mTtAtm,Error)

Numerical = Numerical_Save
PrQ = PrQ_Save
Force_dB = .false.
call mma_deallocate(Coor)

! Stash away this B-matrix for later us with EDiff constraints

if (allocated(BMx_kriging)) call mma_deallocate(BMx_kriging)
call mma_allocate(BMx_kriging,size(BMx,1),size(BMx,2),Label='BMx_kriging')
BMx_kriging(:,:) = BMx(:,:)

if (.not. SaveBMx) BMx(:,:) = BMx_tmp(:,:)

call mma_deallocate(BMx_tmp)

end subroutine NewCar_Kriging
