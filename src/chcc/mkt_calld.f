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
        subroutine MkT_CAlld (T2,X,Y,dimbe,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <-
c           C1                + 1/2 X(be',u,ga',v)
c           C2                + 1/2 X(ga',v,be',u)
c           C3                - 1/2 Y(be',u,ga',v)
c           C4                - 1/2 Y(ga',v,be',u)
c           C5                - 1   Y(ga',u,be',v)
c           C6                - 1   Y(be',v,ga',u)
c        for beGrp=gaGrp
c
        implicit none
        integer dimbe,no
        real*8 T2(1:dimbe,1:dimbe,1:no,1:no)
        real*8 X(1:dimbe,1:no,1:dimbe,1:no)
        real*8 Y(1:dimbe,1:no,1:dimbe,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do ga=1,dimbe
c              do be=1,dimbe  - povodne, stacia iba cleny be>=ga
              do be=ga,dimbe
        T2(be,ga,u,v)=(X(be,u,ga,v)-Y(be,u,ga,v))/2-Y(be,v,ga,u)
              end do
            end do
          end do
        end do
c
        do v=1,no
          do u=1,no
            do be=1,dimbe
c              do ga=1,dimbe  - povodne, stacia iba cleny be>=ga
              do ga=1,be
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
