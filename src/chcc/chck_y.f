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
        subroutine Chck_Y (Y,dimbe,addbe,dimga,addga)
c
c        check Y(be,u,ga,v)
c
        implicit none
#include "chcc1.fh"
        integer dimbe,addbe,dimga,addga
        real*8 Y(1:dimbe,1:no,1:dimga,1:no)
c
        integer be,u,ga,v,bad
        integer a,i
        real*8 s
c
        bad=0
        do v=1,no
        do ga=addga+1,addga+dimga
        do u=1,no
        do be=addbe+1,addbe+dimbe
c
          s=0.0d0
          do i=1,no
          do a=1,nv
           s=s+Kc(i,be,u,a)*T2c(ga,a,i,v)
          end do
          end do
c
          if (abs(Y(be-addbe,u,ga-addga,v)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck Y :',bad
c
        return
        end
