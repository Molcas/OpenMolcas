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

subroutine wzero(n,b,idum)

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, idum
real(kind=wp), intent(_OUT_) :: b(*)
integer(kind=iwp) :: i

do i=1,n
  b(i) = 0
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(idum)

end subroutine wzero
