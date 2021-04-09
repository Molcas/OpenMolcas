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
        subroutine Chck_t2sym
c
c        chek T2c symmetry abij = baji
c
        implicit none
#include "chcc1.fh"
c
        integer i,j,a,b,bad
c
        bad=0
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
          if (abs(T2c(a,b,i,j)-T2c(b,a,j,i)).gt.1.0d-10) then
          bad=bad+1
          end if
        end do
        end do
        end do
        end do
c
        write (6,*) ' T2 Symm Check: ',bad
c
        return
        end
