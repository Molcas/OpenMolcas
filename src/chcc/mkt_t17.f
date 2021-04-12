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
        subroutine MkT_T17 (Ta,T,dimb,dimbe,no)
c
c       this routine do:
c       Ta(i,b',be',u) <- 2 T(be',b',u,i) - T(be',b',i,u)
c
c       N.B. Qvajt odflaknute
c
        implicit none
        integer dimb,dimbe,no
        real*8 T(1:dimbe,1:dimb,1:no,1:no)
        real*8 Ta(1:no,1:dimb,1:dimbe,1:no)
c
c       help variables
        integer i,u,be,b
c
        do u=1,no
        do i=1,no
        do b=1,dimb
        do be=1,dimbe
          Ta(i,b,be,u)=2.0d0*T(be,b,u,i)-T(be,b,i,u)
        end do
        end do
        end do
        end do
c
c
        return
        end
