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
        subroutine MkTau_chcc (T2,T11,T12,dima,dimb,no,f1,f2)
c
c        this routine do:
c        T2(a',b',i,j) = f1 . T2(a',b',i,j) + f2 . T11(a',i) . T12(b',j)
c
c        N.B. Kvajt odflaknute
c
        implicit none
        integer dima,dimb,no
        real*8 f1,f2
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 T11(1:dima,1:no)
        real*8 T12(1:dimb,1:no)
c
c        help variables
        integer i,j,a,b
        real*8 c
c
c
        do j=1,no
          do i=1,no
            do b=1,dimb
              c=t12(b,j)*f2
              do a=1,dima
                t2(a,b,i,j)=f1*t2(a,b,i,j)+c*t11(a,i)
              end do
            end do
          end do
        end do
c
c
        return
        end
