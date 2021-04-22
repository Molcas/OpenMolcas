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
        subroutine MkI_Q47 (Va,V,dimb,dima,no)
c
c       this routine do:
c       Create Va(B',o_A,o_B,A')  needed in step Q48
c       from following available array (permuted as given):
c       V1(B',o_A,o_B,A') = (A',o_A|B',o_B)
c       N.B.
c
c       N.B. Kvajt odflaknute, aj koment k rutine odflaknuty

c
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dimb,1:no,1:no,1:dima)
        real*8 V(1:dimb,1:no,1:no,1:dima)
c
c       help variables
        integer a,b,i,j
c
        do a=1,dima
        do j=1,no
        do i=1,no
        do b=1,dimb
          Va(b,i,j,a)=2.0d0*V(b,j,i,a)-V(b,i,j,a)
        end do
        end do
        end do
        end do
c
        return
        end
