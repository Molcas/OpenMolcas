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

subroutine OFembed_dmat(Dens,nDens)

use OFembed, only: Do_OFemb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDens
real(kind=wp), intent(inout) :: Dens(nDens)
real(kind=wp), allocatable :: D_Var(:)

if (.not. Do_OFemb) return

call NameRun('AUXRFIL') ! switch RUNFILE name

call mma_allocate(D_var,nDens,Label='D_var')
call get_dArray('D1aoVar',D_var,nDens)
Dens(:) = Dens-D_var
call mma_deallocate(D_Var)

call NameRun('#Pop')    ! switch back to old RUNFILE name

return

end subroutine OFembed_dmat
