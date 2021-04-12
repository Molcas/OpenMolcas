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
        subroutine Chck_T17g (V,dima,adda,dimbe,addbe)
c
c        check V(a',u) -
c        sum(be',b,i)  (b,i|be',a') . [ 2 Ta(b,be',i,u) - Ta(b,be',u,i)]
c
        implicit none
#include "chcc1.fh"
        integer dima,adda,dimbe,addbe
        real*8 V(1:dima,1:no)
c
c        help var
        integer a,u,be,b,i,bad,tot
        real*8 s
c
        tot=0
        bad=0
c
        do u=1,no
        do a=adda+1,adda+dima
c
          s=0.0d0
          do b=1,nv
          do be=1+addbe,addbe+dimbe
          do i=1,no
c
           s=s+Q3(be,a,b,i)*(2.0d0*(T2c(be,b,u,i)+T1c(be,u)*T1c(b,i))
     c                         -(T2c(be,b,i,u)+T1c(be,i)*T1c(b,u)))
          end do
          end do
          end do
c
          if (abs(V(a-adda,u)-s).gt.1.0d-10) then
          bad=bad+1
          end if
          tot=tot+1
c
        end do
        end do
c
        write (6,*) ' T17 Chck :',bad,tot
c
        return
        end
