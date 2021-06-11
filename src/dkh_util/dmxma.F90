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

subroutine dmxma(n,transa,transb,a,b,c,alpha)
! Square real matrices multiplication

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
character, intent(in) :: transa, transb
real(kind=wp), intent(in) :: a(n,n), b(n,n), alpha
real(kind=wp), intent(out) :: c(n,n)

call dgemm_(transa,transb,n,n,n,alpha,a,n,b,n,Zero,c,n)

return

end subroutine dmxma
