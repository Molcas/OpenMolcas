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
        subroutine Chck_A (AA)
c
c        check AA(ij,u,v)
c
        implicit none
#include "chcc1.fh"
        real*8 AA(1:no*(no+1)/2,1:no,1:no)
c
c        help var
        integer i,j,ij,u,v,a,b,bad
        real*8 s
c
        bad=0
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
c
        do u=1,no
        do v=1,no
c
          s=Q0(u,i,v,j)
c
          do a=1,nv
           s=s+Q1(a,j,u,i)*T1c(a,v)
          end do
c
          do a=1,nv
           s=s+Q1(a,i,v,j)*T1c(a,u)
          end do
c
          do a=1,nv
          do b=1,nv
           s=s+Q21(a,i,b,j)*(T2c(a,b,u,v)+T1c(a,u)*T1c(b,v))
          end do
          end do
c
          Ac(i,j,u,v)=s
c
c          write (6,99) i,j,u,v,AA(ij,u,v),s,AA(ij,u,v)-s
c99          format (4(i2,1x),3(f15.10,1x))
          if (abs(AA(ij,u,v)-s).gt.1.0d-10) then
          bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
c
        do u=1,no
        do v=1,no
        do i=2,no
        do j=1,i-1
          Ac(j,i,v,u)=Ac(i,j,u,v)
        end do
        end do
        end do
        end do
c
        write (6,*) ' A   Chck :',bad
c
        return
        end
