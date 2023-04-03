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

subroutine pack321(ap,am,b,dimp,dimq,dimr,rc)
! this routine does: B(p,q,r) = A+(p,q,r) - A-(p,r,q) for symq>symr

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimp, dimq, dimr, rc
real(kind=wp) :: ap(dimp,dimq,dimr), am(dimp,dimr,dimq), b(dimp,dimq,dimr)
integer(kind=iwp) :: p, q, r

rc = 0
do r=1,dimr
  do q=1,dimq
    do p=1,dimp
      b(p,q,r) = ap(p,q,r)-am(p,r,q)
    end do
  end do
end do

return

end subroutine pack321
