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

subroutine PCM_Driver(DMat,V,Q,nTs)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTs
real(kind=wp), intent(inout) :: DMat(nTs,nTs)
real(kind=wp), intent(in) :: V(2,nTs)
real(kind=wp), intent(out) :: Q(2,nTs)

! Computes PCM solvation charges given the nuclear and electronic electrostatic
! potential on each tessera.
! Modifies nuclear repulsion, one-electron and two electron terms.

call dgemm_('N','N',2,nTs,nTs,One,V,2,DMat,nTs,Zero,Q,2)

!do i=1,nTs   ! yma delete later
!  write(u6,*) ' == V(2,iTs) diff == ',i,V(2,i)
!end do

return

end subroutine PCM_Driver
