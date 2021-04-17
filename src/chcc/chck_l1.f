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
        subroutine Chck_L1 (L1,dima,adda)
c
c        this routine test L1
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
         real*8 L1(1:nc,1:dima,1:no)
c
c        help var
        integer m,a,bad,i,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do i=1,no
        do a=adda+1,adda+dima
        do m=1,nc
c
           s=L1k(m,i,a)
c
           if (abs(L1(m,a-adda,i)-s).gt.1.0d-10) then
            bad=bad+1
            L1(m,a-adda,i)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
c
        write (6,*) ' L1   ',bad,ntot
c
        return
        end
