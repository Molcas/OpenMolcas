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

subroutine Chck_T1(T1,key)
! check T1

use chcc_global, only: Hooc, Hvoc, Hvvc, no, nv, Q1, Q21, Q22, Q3, T1c, T2c
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: T1(nv,no)
integer(kind=iwp), intent(in) :: key
integer(kind=iwp) :: a, b, bad, be, i, j, u
real(kind=wp) :: s

bad = 0

do u=1,no
  do be=1,nv

    s = Zero

    do a=1,nv
      s = s+Hvvc(be,a)*T1c(a,u)
    end do

    do i=1,no
      s = s-Hooc(i,u)*T1c(be,i)
    end do

    do i=1,no
      do a=1,nv
        s = s+Hvoc(a,i)*(Two*T2c(a,be,i,u)-T2c(a,be,u,i)+T1c(a,u)*T1c(be,i))
      end do
    end do

    do i=1,no
      do a=1,nv
        s = s+(Two*Q21(a,i,be,u)-Q22(a,be,i,u))*T1c(a,i)
      end do
    end do

    do i=1,no
      do b=1,nv
        do a=1,nv
          s = s+(Two*Q3(b,be,a,i)-Q3(a,be,b,i))*(T2c(a,b,i,u)+T1c(a,i)*T1c(b,u))
        end do
      end do
    end do

    do i=1,no
      do j=1,no
        do a=1,nv
          s = s-(Two*Q1(a,i,u,j)-Q1(a,j,u,i))*(T2c(a,be,i,j)+T1c(a,i)*T1c(be,j))
        end do
      end do
    end do

    !s = s/(Oeo(u)-Oev(be))

    if (abs(T1(be,u)-s) > 1.0e-10_wp) then
      bad = bad+1
      if (key == 1) T1(be,u) = s
    end if

  end do
end do

write(u6,*) ' T1 test :',bad

return

end subroutine Chck_T1
