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

subroutine cct3_add32(a,b,q,dimp,dimq,dimr,fact)
! this routine does:
! B(p,q,r) <-- fact * A(p,r) for given q

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: q, dimp, dimq, dimr
real(kind=wp), intent(in) :: a(dimp,dimr), fact
real(kind=wp), intent(inout) :: b(dimp,dimq,dimr)

b(:,q,:) = b(:,q,:)+fact*a

return

end subroutine cct3_add32
