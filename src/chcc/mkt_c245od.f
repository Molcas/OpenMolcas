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
        subroutine MkT_C245od (T2,X,Y,dimbe,dimga,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <<-
c        C2                + 1/2 X(ga',v,be',u)
c        C4                - 1/2 Y(ga',v,be',u)
c        C5                - 1   Y(ga',u,be',v)
c        for beGrp>gaGrp
c
        implicit none
        integer dimbe,dimga,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 X(1:dimga,1:no,1:dimbe,1:no)
        real*8 Y(1:dimga,1:no,1:dimbe,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do be=1,dimbe
              do ga=1,dimga
                 T2(be,ga,u,v)=T2(be,ga,u,v)
     c                       +(X(ga,v,be,u)-Y(ga,v,be,u))/2
     c                       -Y(ga,u,be,v)
              end do
            end do
          end do
        end do
c
        return
        end
