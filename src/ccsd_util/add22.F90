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

subroutine add22(a,b,q,dimp,dimq,fact)
! this routine does:
! B(p,q) <-- fact * A(p) for given q

integer dimp, dimq, q
real*8 fact
real*8 b(1:dimp,1:dimq)
real*8 a(1:dimp)
! help variable
integer p

do p=1,dimp
  b(p,q) = b(p,q)+fact*a(p)
end do

return

end subroutine add22
