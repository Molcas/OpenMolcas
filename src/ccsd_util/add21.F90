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

subroutine add21(a,b,p,dimp,dimq,fact)
! this routine does:
! B(p,q) <-- fact * A(q) for given p

integer dimp, dimq, p
real*8 fact
real*8 b(1:dimp,1:dimq)
real*8 a(1:dimq)
! help variable
integer q

do q=1,dimq
  b(p,q) = b(p,q)+fact*a(q)
end do

return

end subroutine add21
