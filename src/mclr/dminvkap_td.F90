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

subroutine DMInvKap_td(DigPrec,rIn,rOut)
! DigPrec  Diagonal Preconditioner from Prec_td
! rIn      nDensC long orb RHS
! rOut     Trial vector rOut = T^-1B^x

use MCLR_Data, only: nDensC
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: DigPrec(nDensC), rIn(nDensC)
real(kind=wp), intent(out) :: rOut(nDensC)

!-----------------------------------------------------------------------
! Multiply the 1/precond in vector form, rTemp, with RHS in vector form, rIn
! ---> New trial vector.
!-----------------------------------------------------------------------

rOut(:) = rIn(:)/DigPrec(:)

end subroutine DMInvKap_td
