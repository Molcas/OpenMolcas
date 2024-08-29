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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine SetMltplCenters()
! Thomas Bondo Pedersen, July 2012.
!
! Set multipole centers.

use MpmC, only: Coor_MPM
use Sizes_of_Seward, only: S
use Gateway_Info, only: CoM
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i

! Check
if (S%nMltpl < 0) then
  call WarningMessage(2,'SetMltplCenters: illegal input')
  write(u6,'(A,I10)') 'S%nMltpl=',S%nMltpl
  call Abend()
end if

! Allocate center array
call mma_allocate(Coor_MPM,3,S%nMltpl+1,label='Coor_MPM')

! Set origin as center for overlap (0th order)
Coor_MPM(:,1) = Zero

! Set origin as center for dipole (1st order) and
! center of mass as center for higher-order multipoles
if (S%nMltpl > 0) then
  Coor_MPM(:,2) = Zero
  do i=2,S%nMltpl
    Coor_MPM(:,i+1) = CoM(:)
  end do
end if

end subroutine SetMltplCenters
