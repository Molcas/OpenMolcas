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
        subroutine MkT_C136od (T2,X,Y,dimbe,dimga,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <-
c        C1                + 1/2 X(be',u,ga',v)
c        C3                - 1/2 Y(be',u,ga',v)
c        C6                - 1   Y(be',v,ga',u)
c        for beGrp>gaGrp
c
        implicit none
        integer dimbe,dimga,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 X(1:dimbe,1:no,1:dimga,1:no)
        real*8 Y(1:dimbe,1:no,1:dimga,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do ga=1,dimga
              do be=1,dimbe
        T2(be,ga,u,v)=(X(be,u,ga,v)-Y(be,u,ga,v))/2-Y(be,v,ga,u)
              end do
            end do
          end do
        end do
c
        return
        end
