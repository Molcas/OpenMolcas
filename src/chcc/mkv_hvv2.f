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
        subroutine MkV_Hvv2 (Va,V,dima,dimb,no)
c
c        this routine do:
c        Va(b',i,j,a') <- 2(a'i|b'j)-(a'j|b'i)  V(b'I|a'J)
c
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dimb,1:no,1:no,1:dima)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c        help variables
        integer i,j,a,b
c
c
        do a=1,dima
        do j=1,no
        do i=1,no
        do b=1,dimb
          Va(b,i,j,a)=2.0d0*V(b,j,a,i)-V(b,i,a,j)
        end do
        end do
        end do
        end do
c
        return
        end
