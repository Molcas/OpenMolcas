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
        subroutine UpG_T2od (T2,dima,adda,dimb,addb)
c
c        upgrade T2 (off-diagonal)
c
        implicit none
#include "chcc1.fh"
        integer dima,adda,dimb,addb
        real*8 T2(1:dima,1:dimb,1:no,1:no)
c
c        help var
        integer i,j,a,b
c
        do j=1,no
        do i=1,no
          do b=1,dimb
          do a=1,dima
            T2c(a+adda,b+addb,i,j)=T2(a,b,i,j)
            T2c(b+addb,a+adda,j,i)=T2(a,b,i,j)
          end do
          end do
        end do
        end do
c
        return
        end
