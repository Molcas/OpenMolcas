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
        subroutine Chck_Hvv (Hvv)
c
c        check Hoo
c
        implicit none
#include "chcc1.fh"
        real*8 Hvv(1:nv,1:nv)
c
c        help var
        integer i,j,a,b,be,bad
        real*8 s
c
        bad=0
c
        do be=1,nv
        do a=1,nv
c
          s=0.0d0
          do i=1,no
          do j=1,no
          do b=1,nv
          s=s+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(be,b,i,j)+T1c(be,i)*T1c(b,j))
          end do
          end do
          end do
          s=-s
c
          Hvvc(be,a)=s
c
          if (abs(Hvv(a,be)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,*) Hoo(i,u),s,
          end if
c
        end do
        end do
c
        write (6,*) ' Hvv Chck :',bad
c
        return
        end
