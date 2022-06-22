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
        subroutine ExpT2 (T2p,T2u,dima,dimab,no)
c
c        this routine do:
c        Make T2u(a',b',i,j) from T2p(a'b',i,j)
c
c        N.B. Kvajt odflaknute, vypocet ab sa da dat zefektivnit
c
        implicit none
        integer dima,dimab,no
        real*8 T2p(1:dimab,1:no,1:no)
        real*8 T2u(1:dima,1:dima,1:no,1:no)
c
c        help variables
        integer i,j,a,b,ab,ba0
c
c
        do j=1,no
        do i=1,no
          do b=1,dima
            ba0=b*(b-1)/2
            do a=1,b
              T2u(a,b,i,j)=T2p(ba0+a,j,i)
            end do
            do a=1+b,dima
              ab=a*(a-1)/2+b
              T2u(a,b,i,j)=T2p(ab,i,j)
            end do
          end do
        end do
        end do
c
c
        return
        end
