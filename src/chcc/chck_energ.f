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
        subroutine Chck_Energ
c
c        calc energy from T1c and t2c
c
        implicit none
#include "chcc1.fh"
c
c        help var
        integer i,j,a,b
        real*8 e
c
        e=0.0d0
c
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
          e=e+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(a,b,i,j)+T1c(a,i)*T1c(b,j))
        end do
        end do
        end do
        end do
c
        write (6,*) ' Energia Checkeroo',e
c
        return
        end
