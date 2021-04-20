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
        subroutine Chck_Vx (V)
c
c        check V
c
        implicit none
#include "chcc1.fh"
         real*8 V(1:nv,1:no,1:nv,1:no)
c
        integer be,u,i,a,bad
        real*8 s
c
        bad=0
        do a=1,nv
        do i=1,no
        do u=1,no
        do be=1,nv
c
           s=2.0d0*Jc(be,i,u,a)-Kc(i,be,u,a)
           s=Kc(i,be,u,a)
c
          if (abs(V(be,u,a,i)-s).gt.1.0d-10) then
            bad=bad+1
            V(be,u,a,i)=s
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Vx :',bad
c
        return
        end
