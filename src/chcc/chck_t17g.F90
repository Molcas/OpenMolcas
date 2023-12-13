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

subroutine Chck_T17g(V,dima,adda,dimbe,addbe)
! check V(a',u) - sum(be',b,i)  (b,i|be',a') . [ 2 Ta(b,be',i,u) - Ta(b,be',u,i)]

use chcc_global, only: no, nv, Q3, T1c, T2c
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, adda, dimbe, addbe
real(kind=wp), intent(in) :: V(dima,no)
integer(kind=iwp) :: a, b, bad, be, i, tot, u
real(kind=wp) :: s

tot = 0
bad = 0

do u=1,no
  do a=adda+1,adda+dima

    s = Zero
    do b=1,nv
      do be=1+addbe,addbe+dimbe
        do i=1,no

          s = s+Q3(be,a,b,i)*(Two*(T2c(be,b,u,i)+T1c(be,u)*T1c(b,i))-(T2c(be,b,i,u)+T1c(be,i)*T1c(b,u)))
        end do
      end do
    end do

    if (abs(V(a-adda,u)-s) > 1.0e-10_wp) bad = bad+1
    tot = tot+1

  end do
end do

write(u6,*) ' T17 Chck :',bad,tot

return

end subroutine Chck_T17g
