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

subroutine add32(a,b,q,dimp,dimq,dimr,fact)
! this routine does:
! B(p,q,r) <-- fact * A(p,r) for given q

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: q, dimp, dimq, dimr
real(kind=wp) :: a(dimp,dimr), b(dimp,dimq,dimr), fact
integer(kind=iwp) :: p, r

do r=1,dimr
  do p=1,dimp
    b(p,q,r) = b(p,q,r)+fact*a(p,r)
  end do
end do

return

end subroutine add32
