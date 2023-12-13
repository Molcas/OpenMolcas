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

subroutine XDR_mkutls(n,TL,TS,Tr,Bk,A,B,R,UL,US,M1,M2,M3,M4)
! Evaluate transform matrices in non-orthogonal basis space

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: TL(n,n), TS(n,n), Tr(n,n), Bk(n,n), A(n), B(n), R(n)
real(kind=wp), intent(out) :: UL(n,n), US(n,n), M1(n,n), M2(n,n), M3(n,n), M4(n,n)
integer(kind=iwp) :: i

do i=1,n
  M1(:,i) = Tr(:,i)*A(i)
  M2(:,i) = Tr(:,i)*A(i)*R(i)
end do
call dmxma(n,'N','N',M1,TL,M3,One)
call dmxma(n,'N','N',M2,TS,M4,One)
M3(:,:) = M3(:,:)-M4(:,:)
call dmxma(n,'N','N',M3,Bk,UL,One)

do i=1,n
  M1(:,i) = Tr(:,i)*B(i)
  M2(:,i) = Tr(:,i)*B(i)/R(i)
end do
call dmxma(n,'N','N',M1,TL,M3,One)
call dmxma(n,'N','N',M2,TS,M4,One)
M3(:,:) = M3(:,:)+M4(:,:)
call dmxma(n,'N','N',M3,Bk,US,One)

return

end subroutine XDR_mkutls
