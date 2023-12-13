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

subroutine O2N(AA,AB,BB,Temp,nA,nB,Error)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nA, nB
real(kind=wp), intent(in) :: AA(nA,nA), AB(nA,nB)
real(kind=wp), intent(out) :: BB(nB,nB), Temp(nA,nB), Error
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
! (1) Project the new basis on to the old space and compute the
!     matrix elements.

call DGEMM_('N','N',nA,nB,nA,One,AA,nA,AB,nA,Zero,Temp,nA)
call DGEMM_('T','N',nB,nB,nA,One,AB,nA,Temp,nA,Zero,BB,nB)
!                                                                      *
!***********************************************************************
!                                                                      *
! (2) Do diagnostics on how poorly or good the
!     new basis spans the old space

Error = DDot_(nA,AA,nA+1,AA,nA+1)-DDot_(nB,BB,nB+1,BB,nB+1)

return

end subroutine O2N
