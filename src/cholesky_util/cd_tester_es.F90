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

subroutine CD_Tester_ES(X,n,Err)

implicit none
integer n
real*8 X(n,n), Err(3)
integer i

if (n < 1) then
  Err(1) = 9.876543210d15
  Err(2) = -9.876543210d15
  Err(3) = 9.876543210d15
else
  Err(1) = X(1,1)
  Err(2) = X(1,1)
  Err(3) = X(1,1)*X(1,1)
  do i=1,n
    Err(1) = min(Err(1),X(i,i))
    Err(2) = max(Err(2),X(i,i))
    Err(3) = Err(3)+X(i,i)*X(i,i)
  end do
  Err(3) = sqrt(Err(3)/dble(n))
end if

end subroutine CD_Tester_ES
