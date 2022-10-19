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

subroutine cct3_add21(a,b,p,dimp,dimq,fact)
! this routine does:
! B(p,q) <-- fact * A(q) for given p

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: p, dimp, dimq
real(kind=wp) :: a(dimq), b(dimp,dimq), fact
integer(kind=iwp) :: q

do q=1,dimq
  b(p,q) = b(p,q)+fact*a(q)
end do

return

end subroutine cct3_add21
