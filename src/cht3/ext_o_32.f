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
        subroutine ext_o_32 (A,B,nc,no,dima,occ_ind)
c
c this routine do :
c
c extract B (m,a')_i <- A (m,i,a')
c
        implicit none
        integer i1,i2,occ_ind,dima,nc,no
        real*8 A(1:nc,1:no,1:dima),B(1:nc,1:dima)
c
        do i2=1,dima
        do i1=1,nc
c
        B(i1,i2)=A(i1,occ_ind,i2)
c
        end do
        end do
c
        return
        end
