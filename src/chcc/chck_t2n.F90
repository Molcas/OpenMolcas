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

subroutine Chck_T2n(T2,dimbe,addbe,dimga,addga,key)
! chek T2n bez menovatelov, nediagonalne

use chcc_global, only: Ac, Bc, Gooc, Gvvc, Jc, Kc, no, nv, OEo, OEv, Q1, Q21, Q22, Q3, T1c, T2c
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, dimga, addga, key
real(kind=wp), intent(in) :: T2(dimbe,dimga,no,no)
integer(kind=iwp) :: a, b, bad, be, bstart, bup, ga, i, j, ntot, u, v
real(kind=wp) :: s, s1

bad = 0
ntot = 0

do v=1,no
  do u=1,no
    do ga=addga+1,addga+dimga
      if (key == 1) then
        bstart = ga
        bup = addbe+dimbe
      else
        bstart = addbe+1
        bup = addbe+dimbe
      end if
      do be=bstart,bup

        ntot = ntot+1
        s = Zero

        !1
        s1 = Q21(be,u,ga,v)
        s = s+s1

        !2
        s1 = Zero
        do j=1,no
          do i=1,no
            s1 = s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
        end do
        s = s+s1

        !3
        s1 = Zero
        do b=1,nv
          do a=1,nv
            s1 = s1+Bc(a,b,be,ga)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
        end do
        s = s+s1

        !4
        s1 = Zero
        do a=1,nv
          s1 = s1+Gvvc(be,a)*T2c(a,ga,u,v)
          s1 = s1+Gvvc(ga,a)*T2c(a,be,v,u)
        end do
        s = s+s1

        !5
        s1 = Zero
        do i=1,no
          s1 = s1+Gooc(i,u)*T2c(be,ga,i,v)
          s1 = s1+Gooc(i,v)*T2c(ga,be,i,u)
        end do
        s = s-s1

        !6
        s1 = Zero
        do a=1,nv
          s1 = s1+Q3(ga,a,be,u)*T1c(a,v)
          s1 = s1+Q3(be,a,ga,v)*T1c(a,u)
        end do
        do i=1,no
          do a=1,nv
            s1 = s1-(Q22(ga,a,i,u)*T1c(be,i))*T1c(a,v)
            s1 = s1-(Q22(be,a,i,v)*T1c(ga,i))*T1c(a,u)
          end do
        end do
        s = s+s1

        !7
        s1 = Zero
        do i=1,no
          s1 = s1+Q1(be,u,i,v)*T1c(ga,i)
          s1 = s1+Q1(ga,v,i,u)*T1c(be,i)
        end do
        do i=1,no
          do a=1,nv
            s1 = s1+(Q21(be,u,a,i)*T1c(a,v))*T1c(ga,i)
            s1 = s1+(Q21(ga,v,a,i)*T1c(a,u))*T1c(be,i)
          end do
        end do
        s = s-s1

        !8
        s1 = Zero
        do i=1,no
          do a=1,nv
            s1 = s1+(Two*Jc(be,i,u,a)-Kc(i,be,u,a))*(Two*T2c(a,ga,i,v)-T2c(ga,a,i,v))
            s1 = s1+(Two*Jc(ga,i,v,a)-Kc(i,ga,v,a))*(Two*T2c(a,be,i,u)-T2c(be,a,i,u))
          end do
        end do
        s = s+Half*s1

        !9
        s1 = Zero
        do i=1,no
          do a=1,nv
            s1 = s1+Kc(i,be,u,a)*T2c(ga,a,i,v)
            s1 = s1+Kc(i,ga,v,a)*T2c(be,a,i,u)
          end do
        end do
        s = s-Half*s1

        !10
        s1 = Zero
        do i=1,no
          do a=1,nv
            s1 = s1+Kc(i,ga,u,a)*T2c(be,a,i,v)
            s1 = s1+Kc(i,be,v,a)*T2c(ga,a,i,u)
          end do
        end do
        s = s-s1

        s = s/(Oeo(u)+Oeo(v)-Oev(be)-Oev(ga))

        if (abs(T2(be-addbe,ga-addga,u,v)-s) > 1.0e-10_wp) then
          bad = bad+1
          !write(u6,*) 'Bad',abs(T2(be-addbe,ga-addga,u,v)
        end if

      end do
    end do
  end do
end do

if (key == 1) then
  write(u6,*) ' Final test T2 dia',bad,ntot
else
  write(u6,*) ' Final test T2 off',bad,ntot
end if

return

end subroutine Chck_T2n
