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

subroutine DensAB(nBT,nDens,nD,Dens)
!***********************************************************************
!                                                                      *
!     purpose: calculate Density alpha+beta, alpha-beta and dump them  *
!              to runfile                                              *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBT, nDens, nD
real(kind=wp), intent(in) :: Dens(nBT,nD,nDens)
real(kind=wp), allocatable :: Dtemp(:)

if (nD == 1) then
  call Put_dArray('D1ao',Dens(:,1,nDens),nBT)
else
  call mma_allocate(DTemp,nBT,Label='DTemp')

  DTemp(:) = Dens(:,1,nDens)+Dens(:,2,nDens)
  call Put_dArray('D1ao',DTemp,nBT)

  DTemp(:) = Dens(:,1,nDens)-Dens(:,2,nDens)
  call Put_dArray('D1Sao',DTemp,nBT)

  call mma_deallocate(DTemp)
end if

return

end subroutine DensAB
