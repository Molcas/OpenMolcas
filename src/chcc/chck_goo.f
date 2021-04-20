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
        subroutine Chck_Goo (Goo)
c
c        check Goo (i,u)
c
        implicit none
#include "chcc1.fh"
        real*8 Goo(1:no,1:no)
c
c        help var
        integer i,u,j,a,bad
c       integer b
        real*8 s
c
        bad=0
c
        do i=1,no
        do u=1,no
c
c          s=0.0d0
c          do j=1,no
c          do a=1,nv
c          do b=1,nv
c          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(a,b,u,j)
c          end do
c          end do
c          end do
c
           s=Hooc(i,u)
c
          do j=1,no
          do a=1,nv
           s=s+(2.0d0*Q1(a,j,i,u)-Q1(a,i,j,u))*T1c(a,j)
          end do
          end do
c
          Gooc(i,u)=s
c
          if (abs(Goo(i,u)-s).gt.1.0d-10) then
          bad=bad+1
          end if
c
        end do
        end do
c
        write (6,*) ' Goo Chck :',bad
c
        return
        end
