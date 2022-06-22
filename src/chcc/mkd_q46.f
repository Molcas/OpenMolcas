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
        subroutine MkD_Q46 (D,T2,T1a,T1b,dima,dimb,no)
c
c       this routine do:
c       Create D(a',j,b',i) [NB. D(a,b,i,j) =  t2(a,b,j,i)-T2(a,b,i,j)]
c       from following available arrays (permuted as given):
c       T2(a',j,b',i) [NB. T2abij = 1/2 t2abij + tai . tjb ]
c       T1a(a',p)
c       T1b(b',p)
c
c       N.B. Kvajt odflaknute

c
        implicit none
        integer dima,dimb,no
        real*8 D(1:dima,1:no,1:dimb,1:no)
        real*8 T2(1:dima,1:no,1:dimb,1:no)
        real*8 T1a(1:dima,1:no)
        real*8 T1b(1:dimb,1:no)
c
c       help variables
        integer a,b,i,j
c
        do i=1,no
        do b=1,dimb
        do j=1,no
        do a=1,dima
        D(a,j,b,i)=2.0d0*(T2(a,i,b,j)-T1a(a,j)*T1b(b,i))-T2(a,j,b,i)
        end do
        end do
        end do
        end do
c
        return
        end
