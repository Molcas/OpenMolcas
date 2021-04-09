************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine Chck_T2n (T2,dimbe,addbe,dimga,addga,key)
c
c        chek T2n bez menovatelov, nediagonalne
c
        implicit none
        integer dimbe,addbe,dimga,addga,key
#include "chcc1.fh"

         real*8 T2(1:dimbe,1:dimga,1:no,1:no)
c
        integer bstart,bup
        integer u,v,be,ga,bad,ntot
        integer i,j,a,b
        real*8 s,s1
c
        bad=0
        ntot=0
c
        do v=1,no
        do u=1,no
        do ga=addga+1,addga+dimga
        if (key.eq.1) then
             bstart=ga
          bup=addbe+dimbe
        else
          bstart=addbe+1
          bup=addbe+dimbe
        end if
        do be=bstart,bup
c
        ntot=ntot+1
        s=0.0d0
c
c1
          s1=Q21(be,u,ga,v)
          s=s+s1
c
c2
          s1=0.0d0
          do j=1,no
          do i=1,no
          s1=s1+Ac(i,j,u,v)*(T2c(be,ga,i,j)+T1c(be,i)*T1c(ga,j))
          end do
          end do
          s=s+s1
c
c3
          s1=0.0d0
          do b=1,nv
          do a=1,nv
           s1=s1+Bc(a,b,be,ga)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
          end do
          s=s+s1
c
c4
          s1=0.0d0
          do a=1,nv
          s1=s1+Gvvc(be,a)*T2c(a,ga,u,v)
          s1=s1+Gvvc(ga,a)*T2c(a,be,v,u)
          end do
          s=s+s1
c
c5
          s1=0.0d0
          do i=1,no
          s1=s1+Gooc(i,u)*T2c(be,ga,i,v)
          s1=s1+Gooc(i,v)*T2c(ga,be,i,u)
          end do
          s=s-s1
c
c6
          s1=0.0d0
          do a=1,nv
          s1=s1+Q3(ga,a,be,u)*T1c(a,v)
          s1=s1+Q3(be,a,ga,v)*T1c(a,u)
          end do
          do i=1,no
          do a=1,nv
          s1=s1-(Q22(ga,a,i,u)*T1c(be,i))*T1c(a,v)
          s1=s1-(Q22(be,a,i,v)*T1c(ga,i))*T1c(a,u)
          end do
          end do
          s=s+s1
c
c7
          s1=0.0d0
          do i=1,no
          s1=s1+Q1(be,u,i,v)*T1c(ga,i)
          s1=s1+Q1(ga,v,i,u)*T1c(be,i)
          end do
          do i=1,no
          do a=1,nv
          s1=s1+(Q21(be,u,a,i)*T1c(a,v))*T1c(ga,i)
          s1=s1+(Q21(ga,v,a,i)*T1c(a,u))*T1c(be,i)
          end do
          end do
          s=s-s1
c
c8
          s1=0.0d0
          do i=1,no
          do a=1,nv
            s1=s1+(2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a))*
     c          (2.0d0*T2c(a,ga,i,v)-T2c(ga,a,i,v))
            s1=s1+(2.0d0*Jc(ga,i,v,a)-Kc(i,ga,v,a))*
     c          (2.0d0*T2c(a,be,i,u)-T2c(be,a,i,u))

          end do
          end do
          s=s+0.5d0*s1
c
c9
          s1=0.0d0
          do i=1,no
          do a=1,nv
          s1=s1+Kc(i,be,u,a)*T2c(ga,a,i,v)
          s1=s1+Kc(i,ga,v,a)*T2c(be,a,i,u)
          end do
          end do
          s=s-0.5d0*s1
c
c10
          s1=0.0d0
          do i=1,no
          do a=1,nv
          s1=s1+Kc(i,ga,u,a)*T2c(be,a,i,v)
          s1=s1+Kc(i,be,v,a)*T2c(ga,a,i,u)
          end do
          end do
          s=s-s1
c
          s=s/(Oeo(u)+Oeo(v)-Oev(be)-Oev(ga))
c
             if  (abs(T2(be-addbe,ga-addga,u,v)-s).gt.1.0d-10) then
c          write (6,*) 'Bad', abs(T2(be-addbe,ga-addga,u,v)
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        if (key.eq.1) then
          write (6,*) ' Final test T2 dia',bad,ntot
        else
          write (6,*) ' Final test T2 off',bad,ntot
        end if
c
        return
        end
