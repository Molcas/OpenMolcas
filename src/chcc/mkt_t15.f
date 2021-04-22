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
        subroutine MkT_T15 (Tp,T2,T11,T12,dimbe,dima,no)
c
c       this routine do:
c       Tp(be',u,a',i) <- 2 t2(be,a,u,i)
c                       - t2(be,a,i,u)+t12(a,u).t11(be,i)/2
c        N.B. mozno sa to da este vylepsit
c
        implicit none
        integer dimbe,dima,no
        real*8 Tp(1:dimbe,1:no,1:dima,1:no)
        real*8 T2(1:dimbe,1:dima,1:no,1:no)
        real*8 T11(1:dimbe,1:no)
        real*8 T12(1:dima,1:no)
c
c       help variables
        integer be,a,u,i
        real*8 c1

        do i=1,no
          do a=1,dima
            do u=1,no
            c1=T12(a,u)
              do be=1,dimbe
                 Tp(be,u,a,i)=2.0d0*T2(be,a,u,i)
     c                       -T2(be,a,i,u)+c1*T11(be,i)
              end do
            end do
          end do
        end do
c
c
        return
        end
