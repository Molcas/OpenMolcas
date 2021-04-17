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
        subroutine MkV_A4 (Vp,V,dimb,dima,no,dimij)
c
c       this routine do:
c       Vp(a,b,ij) <- (ai|bj) from V(b,j,a,i)
c
        implicit none
        integer dima,dimb,no,dimij
        real*8 Vp(1:dima,1:dimb,1:dimij)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c       help variables
        integer i,j,ij,a,b
c
        ij=0
        do i=1,no
        do j=1,i
          ij=ij+1
          do b=1,dimb
          do a=1,dima
            Vp(a,b,ij)=V(b,j,a,i)
          end do
          end do
        end do
        end do
c
        return
        end
