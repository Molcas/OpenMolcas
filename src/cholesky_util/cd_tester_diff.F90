!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_Diff(X,n,Err)

implicit none
integer n
real*8 X(n*n), Err(3)
integer i
real*8 xn2

if (n < 1) then
  Err(1) = 9.876543210d15
  Err(2) = -9.876543210d15
  Err(3) = 9.876543210d15
else
  Err(1) = x(1)
  Err(2) = x(1)
  Err(3) = x(1)*x(1)
  do i=2,n*n
    Err(1) = min(Err(1),x(i))
    Err(2) = max(Err(2),x(i))
    Err(3) = Err(3)+x(i)*x(i)
  end do
  xn2 = dble(n)*dble(n)
  Err(3) = Err(3)/xn2
end if

end subroutine CD_Tester_Diff
