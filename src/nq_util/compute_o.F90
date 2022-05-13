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

subroutine Compute_O(ZA,RA,nAtoms,T,O,Lambda)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: ZA(nAtoms), RA(3,nAtoms), T(3)
real(kind=wp), intent(out) :: O(3,3), Lambda(3)
real(kind=wp) :: EVal(6), M(3,3)

!                                                                      *
!***********************************************************************
!                                                                      *
! Form the nuclear charge moment tensor

call Compute_M(ZA,nAtoms,RA,T,M)
!                                                                      *
!***********************************************************************
!                                                                      *
! Diagonalize the nuclear charge momentum tensor to get
! the principal axis system.

O = reshape([One,Zero,Zero,Zero,One,Zero,Zero,Zero,One],[3,3])
EVal(1) = M(1,1)
EVal(2) = M(2,1)
EVal(3) = M(2,2)
EVal(4) = M(3,1)
EVal(5) = M(3,2)
EVal(6) = M(3,3)
call Jacob(EVal,O,3,3)
!call JacOrd(EVal,O,3,3)
#ifdef _DEBUGPRINT_
call TriPrt('RotGrd: EVal',' ',EVal,3)
call RecPrt('RotGrd: O',' ',O,3,3)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
Lambda(1) = EVal(1)
Lambda(2) = EVal(3)
Lambda(3) = EVal(6)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_O
