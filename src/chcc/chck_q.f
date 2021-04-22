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
        subroutine Chck_Q (Q,dimbe,addbe,dima,adda)
c
c        check Q(be,u,i,a)
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dima,adda
        real*8 Q(1:dimbe,1:no,1:no,1:dima)
c
        integer be,u,i,a,bad
        real*8 s,sk,sj
c
        bad=0
c
        do a=adda+1,adda+dima
        do i=1,no
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          sj=Jc(be,i,u,a)
          sk=Kc(i,be,u,a)
c
          s=2*sj-sk
c
          if (abs(Q(be-addbe,u,i,a-adda)-s).gt.1.0d-10) then
            bad=bad+1
c           Q(be-addbe,u,i,a-adda)=s
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Q :',bad
c
        return
        end
