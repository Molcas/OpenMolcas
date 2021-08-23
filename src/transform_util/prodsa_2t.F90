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

subroutine ProdsA_2t(AB,iA,iB,CMO,nMO,Y)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iA, iB, nMO
real(kind=wp), intent(in) :: AB(iA,iB), CMO(iA,nMO)
real(kind=wp), intent(out) :: Y(iB,nMO)

call DGEMM_('T','N',iB,nMO,iA,One,AB,iA,CMO,iA,Zero,Y,iB)

return

end subroutine ProdsA_2t
