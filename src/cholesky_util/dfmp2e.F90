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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine DfMP2E(NVar,NOcc,NCore,E,EMin,EMax)
!-----------------------------------------------------------------------
! Function : Define MP2 energy
!-----------------------------------------------------------------------

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NVar, NOcc, NCore
real(kind=wp), intent(in) :: E(NVar)
real(kind=wp), intent(out) :: EMin, EMax
integer(kind=iwp) :: I
integer(kind=iwp), allocatable :: IDX(:)

call mma_allocate(IDX,NVar,Label='IDX')
do I=1,NVar
  IDX(I) = I
end do
!call OrdV(E,IDX,NVar)

EMin = Two*(E(IDX(NCore+NOcc+1))-E(IDX(NCore+NOcc)))
EMax = Two*(E(IDX(NVar))-E(IDX(NCore+1)))
call mma_deallocate(IDX)

return

end subroutine DfMP2E
