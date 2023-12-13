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

subroutine Order_Arrays(mode,A1,N1,N2,P2,SCR1)

use Definitions, only: wp, iwp, u6

implicit none
character(len=4), intent(in) :: mode
integer(kind=iwp), intent(in) :: N1, N2
real(kind=wp), intent(inout) :: A1(N1,N2), P2(N2)
real(kind=wp), intent(out) :: SCR1(N1)
integer(kind=iwp) :: j, k
real(kind=wp) :: tmp

if (mode == 'decr') then
  do j=1,N2-1
    do k=j+1,N2
      if (P2(j) < P2(k)) then
        tmp = P2(j)
        P2(j) = P2(k)
        P2(k) = tmp
        SCR1(:) = A1(:,j)
        A1(:,j) = A1(:,k)
        A1(:,k) = SCR1(:)
      end if
    end do
  end do
else if (mode == 'incr') then
  do j=1,N2-1
    do k=j+1,N2
      if (P2(j) > P2(k)) then
        tmp = P2(j)
        P2(j) = P2(k)
        P2(k) = tmp
        SCR1(:) = A1(:,j)
        A1(:,j) = A1(:,k)
        A1(:,k) = SCR1(:)
      end if
    end do
  end do
else
  write(u6,*) ' In routine Order_Arrays: wrong mode!'
  call Abend()
end if

return

end subroutine Order_Arrays
