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

subroutine Chck_Xred(X,dimbe,addbe,dimga,addga)
! check X(be,u,ga,v)
! do X clenu nieje zahrnuty prispevok od Aex, preto picuje

implicit none
#include "chcc1.fh"
integer dimbe, addbe, dimga, addga
real*8 X(1:dimbe,1:no,1:dimga,1:no)
integer be, u, ga, v, bad
integer a, i, j
real*8 s, s1

bad = 0
do v=1,no
  do ga=addga+1,addga+dimga
    do u=1,no
      do be=addbe+1,addbe+dimbe

        s = 0.0d0

        ! X2 (T24)   + Gvv(a,be)   . t2(ga,a,v,u)
        s1 = 0.0d0
        do a=1,nv
          s1 = s1+Gvvc(be,a)*T2c(a,ga,u,v)
        end do
        s = s+2.0d0*s1

        ! X3 (T25)   - Goo(i,u)    . t2(ga,be,v,i)
        s1 = 0.0d0
        do i=1,no
          s1 = s1+Gooc(i,u)*T2c(ga,be,v,i)
        end do
        s = s-2.0d0*s1

        ! T1 cleny
        s1 = 0.0d0

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

        !red s = s+2.0d0*s1

        ! X4 (T22)   + sum(i,j)   [ Ta(be,ga,i,j) . A(i,j,u,v)  ]
        s1 = 0.0d0
        do i=1,no
          do j=1,no
            s1 = s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
        end do
        s = s+s1

        ! X1   <-  Q(be,u,i,a) . (2 t2(ga,a,v,i) - t2(ga,a,i,v))
        ! calc as (2J(be,i,u,a)-K(i,be,u,a)*(2t2(a,ga,i,v)-t2(ga,a,i,v)
        s1 = 0.d0
        do i=1,no
          do a=1,nv
            s1 = s1+(2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a))*(2.0d0*T2c(a,ga,i,v)-T2c(ga,a,i,v))
          end do
        end do
        s = s+s1

        ! <- V1(be',u,ga',v)
        s1 = Q21(be,u,ga,v)
        s = s+s1

        if (abs(X(be-addbe,u,ga-addga,v)-s) > 1.0d-10) bad = bad+1

      end do
    end do
  end do
end do

write(6,*) ' Chck X :',bad

return

end subroutine Chck_Xred
