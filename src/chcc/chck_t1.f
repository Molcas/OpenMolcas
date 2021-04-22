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
        subroutine Chck_T1 (T1,key)
c
c        check T1
c
        implicit none
#include "chcc1.fh"
        real*8 T1(1:nv,1:no)
c
c        help var
        integer be,u,i,j,a,b,bad
        integer key
        real*8 s
c
        bad=0
c
        do u=1,no
        do be=1,nv
c
          s=0.0d0
c
          do a=1,nv
          s=s+Hvvc(be,a)*T1c(a,u)
          end do
c
          do i=1,no
          s=s-Hooc(i,u)*T1c(be,i)
          end do
c
          do i=1,no
          do a=1,nv
             s=s+Hvoc(a,i)*
     c      (2.0d0*T2c(a,be,i,u)-T2c(a,be,u,i)+T1c(a,u)*T1c(be,i))
          end do
          end do
c
          do i=1,no
          do a=1,nv
           s=s+(2.0d0*Q21(a,i,be,u)-Q22(a,be,i,u))*T1c(a,i)
          end do
          end do
c
          do i=1,no
          do b=1,nv
          do a=1,nv
             s=s+(2.0d0*Q3(b,be,a,i)-Q3(a,be,b,i))*
     c          (T2c(a,b,i,u)+T1c(a,i)*T1c(b,u))
          end do
          end do
          end do
c
          do i=1,no
          do j=1,no
          do a=1,nv
             s=s-(2.0d0*Q1(a,i,u,j)-Q1(a,j,u,i))*
     c          (T2c(a,be,i,j)+T1c(a,i)*T1c(be,j))
          end do
          end do
          end do
c
c            s=s/(Oeo(u)-Oev(be))
c
          if (abs(T1(be,u)-s).gt.1.0d-10) then
            bad=bad+1
            if (key.eq.1) then
              T1(be,u)=s
            end if
          end if
c
        end do
        end do
c
         write (6,*) ' T1 test :',bad
c
        return
        end
c
