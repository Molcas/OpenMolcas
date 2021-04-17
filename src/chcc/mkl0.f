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
        subroutine MkL0 (V)
c
c        L0(m,i,j) <- V(m,ij)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nc,1:no*(no+1)/2)
c
c        help var
        integer i,j,ij,m
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
          do m=1,nc
             L0k(m,i,j)=V(m,ij)
             L0k(m,j,i)=V(m,ij)
          end do
        end do
        end do
c
        return
        end
