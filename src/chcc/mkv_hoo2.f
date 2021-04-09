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
        subroutine MkV_Hoo2 (V2,V,dima,dimb,no)
c
c       this routine do:
c       Make AntiSymetric integrals
c        V2(i,a',b',j) <- 2 V(b'i|a'j) - V(b'j|a'i)
c
        implicit none
        integer dimb,dima,no
        real*8 V2(1:no,1:dima,1:dimb,1:no)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c       help variables
        integer a,b,i,j
c
        do j=1,no
        do b=1,dimb
        do a=1,dima
        do i=1,no
          V2(i,a,b,j)=2.0d0*V(b,j,a,i)-V(b,i,a,j)
        end do
        end do
        end do
        end do
c
c
        return
        end
