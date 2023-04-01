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

integer dimp, dimq, dimr, rc
real*8 ap(1:dimp,1:dimq,1:dimr)
real*8 am(1:dimp,1:dimr,1:dimq)
real*8 b(1:dimp,1:dimq,1:dimr)
! help variables
integer p, q, r

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
