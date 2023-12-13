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

subroutine XDR_fpFWprop(n,Tr,X,pXp,A,B,R,EL,ES,OL,OS,tmp)
!  Transform property operator (X,pXp) to fpFW picture

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Tr(n,n), A(n), B(n), R(n)
real(kind=wp), intent(inout) :: X(n,n), pXp(n,n)
real(kind=wp), intent(out) :: EL(n,n), ES(n,n), OL(n,n), OS(n,n), tmp(n,n)
integer(kind=iwp) :: i, j
real(kind=wp) :: av, aw

call dmxma(n,'C','N',Tr,X,tmp,One)
call dmxma(n,'N','N',tmp,Tr,X,One)
call dmxma(n,'C','N',Tr,pXp,tmp,One)
call dmxma(n,'N','N',tmp,Tr,pXp,One)

do i=1,n
  do j=1,n
    av = X(j,i)*A(i)*A(j)
    aw = pXp(j,i)*B(i)*B(j)
    EL(j,i) = av+aw
    ES(j,i) = aw/R(i)/R(j)+av*R(i)*R(j)
    OL(j,i) = aw/R(i)-av*R(i)
    OS(j,i) = aw/R(j)-av*R(j)
  end do
end do

return

end subroutine XDR_fpFWprop
