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
        subroutine Chck_Hoo (Hoo)
c
c        check Hoo
c
        implicit none
#include "chcc1.fh"
        real*8 Hoo(1:no,1:no)
c
c        help var
        integer i,u,j,a,b,bad
        real*8 s
c
        bad=0
c
        do i=1,no
        do u=1,no
c
          s=0.0d0
          do j=1,no
          do a=1,nv
          do b=1,nv
          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(a,b,u,j)+T1c(a,u)*T1c(b,j))
          end do
          end do
          end do
c
          Hooc(i,u)=s
c
          if (abs(Hoo(i,u)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,*) Hoo(i,u),s
          end if
c
        end do
        end do
c
        write (6,*) ' Hoo Chck :',bad
c
        return
        end
