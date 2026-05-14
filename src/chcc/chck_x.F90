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

subroutine Chck_X(X,dimbe,addbe,dimga,addga)
! check X(be,u,ga,v)
! do X clenu nieje zahrnuty prispevok od Aex, preto picuje

use chcc_global, only: Ac, Gooc, Gvvc, Jc, Kc, no, nv, Q1, Q21, Q22, Q3, T1c, T2c
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, dimga, addga
real(kind=wp), intent(in) :: X(dimbe,no,dimga,no)
integer(kind=iwp) :: a, bad, be, ga, i, j, u, v
real(kind=wp) :: s, s1

bad = 0
do v=1,no
  do ga=addga+1,addga+dimga
    do u=1,no
      do be=addbe+1,addbe+dimbe

        s = Zero

        ! X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
        s1 = Zero
        do a=1,nv
          s1 = s1+Gvvc(be,a)*T2c(a,ga,u,v)
        end do
        s = s+Two*s1

        ! X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
        s1 = Zero
        do i=1,no
          s1 = s1+Gooc(i,u)*T2c(ga,be,v,i)
        end do
        s = s-Two*s1

        ! T1 cleny
        s1 = Zero

        do a=1,nv
          s1 = s1+Q3(a,ga,be,u)*T1c(a,v)
        end do

        do i=1,no
          s1 = s1-Q1(be,u,i,v)*T1c(ga,i)
        end do

        do i=1,no
          do a=1,nv
            s1 = s1-(Q22(a,ga,i,u)*T1c(be,i))*T1c(a,v)
            s1 = s1-(Q21(a,i,be,u)*T1c(a,v))*T1c(ga,i)
          end do
        end do

        s = s+Two*s1

        ! X4 (T22)   + sum(i,j)   [ Ta(be,ga,i,j) . A(i,j,u,v)  ]
        s1 = Zero
        do i=1,no
          do j=1,no
            s1 = s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
        end do
        s = s+s1

        ! X1   <-  Q(be,u,i,a) . (2 t2(ga,a,v,i) - t2(ga,a,i,v))
        ! calc as (2J(be,i,u,a)-K(i,be,u,a)*(2t2(a,ga,i,v)-t2(ga,a,i,v)
        s1 = Zero
        do i=1,no
          do a=1,nv
            s1 = s1+(Two*Jc(be,i,u,a)-Kc(i,be,u,a))*(Two*T2c(a,ga,i,v)-T2c(ga,a,i,v))
          end do
        end do
        s = s+s1

        ! <- V1(be',u,ga',v)
        s1 = Q21(be,u,ga,v)
        s = s+s1

        if (abs(X(be-addbe,u,ga-addga,v)-s) > 1.0e-10_wp) bad = bad+1

      end do
    end do
  end do
end do

write(u6,*) ' Chck X :',bad

return

end subroutine Chck_X
